/* CAT Bridge web app: inputs, checks, worker, progress, results. All data stay in the browser. */
(function () {
  'use strict';

  const WHEEL = 'catbridge-2.0.0-py3-none-any.whl';
  const ROLES = ['gene', 'metabolite', 'design', 'annotation'];
  const ROLE_LABEL = { gene: 'Gene table', metabolite: 'Metabolite table', design: 'Design table', annotation: 'Annotation table' };
  const EXAMPLE = { gene: 'gene.tsv', metabolite: 'metabolite.tsv', design: 'design.tsv', annotation: 'annotation.tsv' };
  const EXAMPLE_TARGET = 'Capsaicin';
  const METRIC_HINT = {
    Granger: 'Granger causality test; the p-value is used (lower is stronger).',
    CCM: 'Convergent cross mapping (higher is stronger).',
    CCA: 'Canonical correlation (higher is stronger).',
    DTW: 'Dynamic time warping distance (lower is stronger).',
    CCF: 'Cross-correlation within the lag (higher is stronger).',
    Spearman: 'Rank correlation (higher is stronger).',
    Pearson: 'Linear correlation (higher is stronger).',
  };
  const STEP_LABEL = {
    runtime: 'Load Python', packages: 'Load numpy, pandas, SciPy, scikit-learn, statsmodels', install: 'Install CAT Bridge',
    deseq: 'Load DESeq2', read: 'Read and check tables', fc: 'Compute fold change', metric: 'Compute metric', cluster: 'Cluster time series',
    annot: 'Add annotation', score: 'Score and rank genes',
  };
  const $ = (id) => document.getElementById(id);
  const fmtInt = (n) => Number(n).toLocaleString('en-US');
  const fmtBytes = (b) => (b < 1024 ? b + ' B' : b < 1048576 ? (b / 1024).toFixed(0) + ' KB' : (b / 1048576).toFixed(1) + ' MB');
  const fmtSec = (s) => (s < 60 ? s.toFixed(s < 10 ? 1 : 0) + ' s' : Math.floor(s / 60) + ' min ' + String(Math.round(s % 60)).padStart(2, '0') + ' s');
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  const S = {
    files: {}, infos: {}, target: '', metabolites: [],
    worker: null, baseReady: false, deseqReady: false, prepError: null, prepDetail: '',
    steps: {}, running: false, runId: 0, runStart: 0, ticker: null, runPlan: null,
    check: null, post: {}, result: null,
  };

  /* ---- messages ---------------------------------------------------------------------------------------- */
  function renderMsgs(ul, items) {
    ul.innerHTML = '';
    items.slice().sort((a, b) => (a.level === b.level ? 0 : a.level === 'error' ? -1 : 1)).forEach((it) => {
      const li = document.createElement('li');
      li.className = 'msg msg-' + (it.level === 'error' ? 'error' : it.level === 'ok' ? 'ok' : 'warn');
      li.innerHTML = '<span class="kind">' + (it.level === 'error' ? 'Error' : it.level === 'ok' ? 'OK' : 'Warning') + '</span> <span class="body">' + esc(it.text) + '</span>';
      ul.appendChild(li);
    });
  }

  /* ---- worker ------------------------------------------------------------------------------------------ */
  function newWorker() {
    if (S.worker) S.worker.terminate();
    const w = new Worker('worker.js?v=7', { type: 'module' });
    S.worker = w;
    S.baseReady = false; S.deseqReady = false; S.prepError = null; S.prepDetail = '';
    S.steps = {};
    w.onmessage = (ev) => onWorker(ev.data, w);
    w.onerror = (e) => onWorker({ type: 'fatal', message: e.message || 'The worker failed to start.' }, w);
    return w;
  }
  function warmUp() {
    if (!S.worker) newWorker();
    S.worker.postMessage({ type: 'init', wheelUrl: WHEEL, deseq: !!(S.files.design && S.infos.design && !S.infos.design.fatal) });
    renderRuntime();
  }
  function setStep(key, detail) {
    const order = ['runtime', 'packages', 'install', 'deseq', 'read', 'fc', 'metric', 'cluster', 'annot', 'score'];
    const idx = order.indexOf(key);
    if (idx < 0) return;
    order.forEach((k, i) => {
      const st = S.steps[k] || (S.steps[k] = { state: 'pending', detail: '' });
      if (i < idx && st.state !== 'done') { st.state = 'done'; st.t1 = performance.now(); }
    });
    const st = S.steps[key] || (S.steps[key] = { state: 'pending', detail: '' });
    if (st.state !== 'active') { st.state = 'active'; st.t0 = performance.now(); }
    if (detail && ['read', 'cluster'].includes(key) === false && key !== 'metric' && key !== 'fc') detail = '';
    if (detail) st.detail = key === 'metric' ? '' : detail;
  }
  function finishStep(key) {
    const st = S.steps[key] || (S.steps[key] = { state: 'pending', detail: '' });
    if (st.state !== 'done') { st.state = 'done'; st.t1 = performance.now(); }
  }
  function onWorker(m, w) {
    if (w !== S.worker) return;
    if (m.type === 'progress') {
      if (m.step === 'fc-detail') { if (S.steps.fc) S.steps.fc.detail = m.detail; }
      else if (m.step === 'packages') { setStep('packages'); }
      else if (m.step === 'metric') { setStep('metric', m.detail); }
      else if (m.step === 'read') { setStep('read'); }
      else setStep(m.step, m.detail);
      renderRuntime(); renderProgress();
    } else if (m.type === 'ready') {
      if (m.stage === 'base') { ['runtime', 'packages', 'install'].forEach(finishStep); S.baseReady = true; }
      if (m.stage === 'deseq') { finishStep('deseq'); S.deseqReady = true; S.versions = m.versions; }
      renderRuntime(); renderProgress();
    } else if (m.type === 'notice') {
      S.notice = m.text;
    } else if (m.type === 'result') {
      onResult(m);
    } else if (m.type === 'fatal') {
      S.prepError = m.message;
      if (S.running) fail('Could not load or run the Python environment: ' + m.message + ' Check your internet connection and try again.', m.stack || m.message, true);
      renderRuntime();
    }
  }
  function renderRuntime() {
    const n = $('runtime-note');
    if (!n) return;
    if (S.prepError) n.textContent = 'Python could not be loaded. Check your internet connection and reload the page.';
    else if (!S.worker) n.textContent = '';
    else if (S.baseReady && (S.deseqReady || !S.files.design)) n.textContent = 'Python is ready.';
    else if (S.baseReady) n.textContent = 'Loading DESeq2 in the background.';
    else {
      const cur = ['runtime', 'packages', 'install'].find((k) => S.steps[k] && S.steps[k].state === 'active');
      n.textContent = 'Loading Python in the background' + (cur ? ': ' + STEP_LABEL[cur].toLowerCase() : '') + '. The first visit downloads about 45 MB, later visits use the cache.';
    }
  }

  /* ---- progress ---------------------------------------------------------------------------------------- */
  function plannedSteps() {
    const p = S.runPlan; if (!p) return [];
    const list = [];
    if (!p.baseWasReady) list.push('runtime', 'packages', 'install');
    if (p.design && !p.deseqWasReady) list.push('deseq');
    list.push('read', 'fc', 'metric');
    if (p.clusters) list.push('cluster');
    if (p.annotation) list.push('annot');
    list.push('score');
    return list;
  }
  function stepLabel(k) {
    if (k === 'metric') return 'Compute ' + S.runPlan.metric;
    if (k === 'fc') return S.runPlan.design ? 'Compute fold change with DESeq2' : 'Compute fold change';
    return STEP_LABEL[k];
  }
  function renderProgress() {
    if (!S.running) return;
    const list = plannedSteps(); const ol = $('steps');
    let done = 0; let html = '';
    list.forEach((k) => {
      const st = S.steps[k] || { state: 'pending', detail: '' };
      if (st.state === 'done') done++;
      const dur = st.state === 'done' && st.t0 ? ' ' + fmtSec((st.t1 - st.t0) / 1000) : '';
      let detail = st.detail || '';
      if (k === 'fc' && S.runPlan.contrast) detail = 'log2(' + S.runPlan.contrast + ')' + (st.state === 'active' && st.detail && !/^log2/.test(st.detail) ? ': ' + st.detail : '');
      if (k === 'fc' && st.state === 'active' && S.runPlan.design && !S.steps.fc.detail) detail = detail || '';
      html += '<li class="' + st.state + '"><span class="dot"></span><span>' + esc(stepLabel(k)) + (detail && st.state !== 'pending' ? ' <span class="detail">' + esc(detail) + '</span>' : '') + '</span><span class="time">' + esc(dur) + '</span></li>';
    });
    ol.innerHTML = html;
    const active = list.find((k) => (S.steps[k] || {}).state === 'active');
    const pct = Math.min(99, Math.round(((done + (active ? 0.5 : 0)) / list.length) * 100));
    const bar = $('progress-bar'); bar.setAttribute('aria-valuenow', String(pct)); bar.firstElementChild.style.width = pct + '%';
    $('progress-title').textContent = active ? stepLabel(active) + '...' : 'Starting...';
    let note = '';
    if (list.includes('runtime')) note += 'First visit: Python and its packages (about 45 MB) are downloaded once and cached by your browser. ';
    if (S.runPlan.design) note += (S.runPlan.genes > 10000 ? 'DESeq2 on ' + fmtInt(S.runPlan.genes) + ' genes can take several minutes. ' : 'DESeq2 can take a minute or more. ') + 'The page stays responsive and you can cancel at any time.';
    $('progress-note').textContent = note.trim();
  }
  function tick() { $('progress-elapsed').textContent = 'Elapsed ' + fmtSec((performance.now() - S.runStart) / 1000); }

  /* ---- failure ----------------------------------------------------------------------------------------- */
  function slotOf(text) {
    for (const r of ROLES) if (S.files[r] && text.indexOf(S.files[r].name) >= 0) return r;
    const t = text.toLowerCase();
    if (/annotation/.test(t)) return 'annotation';
    if (/design/.test(t) || /\bgroup\b/.test(t)) return 'design';
    if (/target|metabolite/.test(t)) return /target/.test(t) ? 'target' : 'metabolite';
    if (/gene/.test(t)) return 'gene';
    return null;
  }
  function fail(message, trace, fatal) {
    stopRun();
    $('progress').hidden = true; $('results').hidden = true;
    $('failure').hidden = false;
    renderMsgs($('failure-msgs'), [{ level: 'error', text: message }]);
    $('failure-trace').textContent = trace || '';
    $('failure').querySelector('details').hidden = !trace;
    if (!fatal) {
      const slot = slotOf(message);
      if (slot) { S.post[slot] = (S.post[slot] || []).concat([{ level: 'error', text: message }]); renderSlotMsgs(); }
    }
    $('failure').scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    updateRunState();
  }

  /* ---- slots and checks -------------------------------------------------------------------------------- */
  function slotEl(role) { return document.querySelector('.slot[data-role="' + role + '"]'); }

  async function setFile(role, file) {
    clearOutputs();
    S.files[role] = file; S.infos[role] = null;
    const slot = slotEl(role); const drop = slot.querySelector('.drop');
    drop.classList.add('has-file'); drop.classList.remove('is-error');
    slot.querySelector('.drop-main').textContent = file.name;
    slot.querySelector('.drop-meta').textContent = fmtBytes(file.size) + ' · checking...';
    slot.querySelector('[data-remove]').hidden = false; slot.querySelector('.tag').hidden = true;
    const my = file;
    const info = await CBChecks.analyzeFile(role, file);
    if (S.files[role] !== my) return;
    S.infos[role] = info;
    if (role === 'metabolite') setMetabolites(info);
    if (!info.fatal && !S.worker) warmUp();
    else if (!info.fatal && role === 'design') warmUp();
    recheck();
  }
  function removeFile(role) {
    clearOutputs();
    delete S.files[role]; delete S.infos[role];
    const slot = slotEl(role); const drop = slot.querySelector('.drop');
    drop.classList.remove('has-file', 'is-error');
    slot.querySelector('.drop-main').innerHTML = '<strong>Choose a file</strong> or drop it here';
    slot.querySelector('.drop-meta').textContent = '.tsv, .csv or .txt';
    slot.querySelector('[data-remove]').hidden = true; slot.querySelector('.tag').hidden = false;
    slot.querySelector('input[type=file]').value = '';
    if (role === 'metabolite') setMetabolites(null);
    recheck();
  }
  function describe(info) {
    if (!info || info.fatal && !info.header.length) return '';
    const h = info.header; const cols = h.length - 1;
    if (info.role === 'gene') return fmtInt(info.nRows) + ' genes · ' + cols + ' columns';
    if (info.role === 'metabolite') return fmtInt(info.nRows) + ' metabolites · ' + cols + ' columns';
    if (info.role === 'design') return fmtInt(info.nRows) + ' samples' + (info.groups ? ' · ' + new Set(info.groups).size + ' groups' : '');
    return fmtInt(info.nRows) + ' rows';
  }

  function state() {
    return {
      infos: S.infos, target: S.target, metric: $('metric').value,
      params: { clusters: $('p-clusters').value, lag: $('p-lag').value, E: $('p-E').value, tau: $('p-tau').value },
    };
  }
  function recheck() {
    S.check = CBChecks.crossChecks(state());
    renderSlotMsgs();
    updateRunState();
  }
  function renderSlotMsgs() {
    const c = S.check || { slots: {} };
    ROLES.concat(['target', 'params']).forEach((r) => {
      const info = S.infos[r];
      const items = (info ? info.issues : []).concat(c.slots[r] || [], S.post[r] || []);
      renderMsgs($('msgs-' + r), items);
      const hasErr = items.some((i) => i.level === 'error');
      if (ROLES.includes(r)) {
        const slot = slotEl(r); const drop = slot.querySelector('.drop');
        drop.classList.toggle('is-error', hasErr && !!S.files[r]);
        if (S.files[r]) {
          slot.querySelector('.drop-meta').textContent = fmtBytes(S.files[r].size) + (info ? (describe(info) ? ' · ' + describe(info) : '') : ' · checking...');
        }
      }
      if (r === 'target') $('target').setAttribute('aria-invalid', hasErr ? 'true' : 'false');
      if (r === 'params' && hasErr) $('adv').open = true;
    });
    const hint = $('target-hint');
    hint.textContent = c.contrast
      ? 'Fold change = log2(' + c.contrast.end + ' / ' + c.contrast.start + '), the interval in which ' + S.target + ' rises fastest.'
      : 'Search the rows of the metabolite table.';
  }
  function updateRunState() {
    const c = S.check || { slots: {}, todo: [] };
    let nErr = 0; const where = [];
    ROLES.concat(['target', 'params']).forEach((r) => {
      const items = ((S.infos[r] ? S.infos[r].issues : []).concat(c.slots[r] || [], S.post[r] || [])).filter((i) => i.level === 'error');
      if (items.length) { nErr += items.length; where.push(r === 'params' ? 'Advanced' : r === 'target' ? 'Target metabolite' : ROLE_LABEL[r]); }
    });
    const todo = (c.todo || []).slice();
    const btn = $('btn-run'); const st = $('run-status');
    const busy = S.running;
    btn.disabled = busy || nErr > 0 || todo.length > 0;
    btn.setAttribute('aria-disabled', btn.disabled ? 'true' : 'false');
    $('btn-cancel').hidden = !busy;
    btn.hidden = busy;
    st.classList.toggle('has-error', nErr > 0);
    let text = '';
    if (busy) text = 'Running. You can cancel at any time.';
    else if (nErr) text = 'Fix ' + nErr + (nErr === 1 ? ' error' : ' errors') + ' before running: ' + where.join(', ') + '.';
    else if (todo.length) text = 'To run: ' + todo.join(', ') + '.';
    else text = 'Ready.' + ((c.slots.gene || []).concat(c.slots.metabolite || [], c.slots.design || [], c.slots.annotation || [], c.slots.target || []).some((i) => i.level === 'warn') ? ' Warnings do not block the run.' : '');
    st.textContent = text;
  }
  function clearOutputs() {
    S.post = {};
    if (!S.running) { $('failure').hidden = true; }
  }

  /* ---- target combobox --------------------------------------------------------------------------------- */
  function setMetabolites(info) {
    const input = $('target');
    S.metabolites = info && !info.fatal ? Array.from(new Set(info.names)) : [];
    input.disabled = !S.metabolites.length;
    input.placeholder = S.metabolites.length ? 'Search ' + fmtInt(S.metabolites.length) + ' metabolites' : 'Add the metabolite table first';
    if (!S.metabolites.includes(S.target)) { S.target = ''; input.value = ''; }
  }
  const combo = { open: false, items: [], active: -1 };
  function comboFilter(q) {
    q = q.trim().toLowerCase();
    if (!q) return S.metabolites.slice();
    const starts = []; const has = [];
    S.metabolites.forEach((n) => { const l = n.toLowerCase(); if (l.startsWith(q)) starts.push(n); else if (l.includes(q)) has.push(n); });
    return starts.concat(has);
  }
  function comboRender() {
    const ul = $('target-list'); const input = $('target');
    const all = comboFilter(input.value === S.target ? '' : input.value);
    combo.items = all.slice(0, 50);
    ul.innerHTML = '';
    if (!all.length) { ul.innerHTML = '<li class="more" role="presentation">No metabolite matches</li>'; }
    combo.items.forEach((n, i) => {
      const li = document.createElement('li'); li.id = 'opt-' + i; li.setAttribute('role', 'option'); li.textContent = n;
      li.setAttribute('aria-selected', i === combo.active ? 'true' : 'false');
      li.addEventListener('mousedown', (e) => { e.preventDefault(); comboPick(n); });
      ul.appendChild(li);
    });
    if (all.length > 50) { const li = document.createElement('li'); li.className = 'more'; li.setAttribute('role', 'presentation'); li.textContent = fmtInt(all.length - 50) + ' more: keep typing to narrow the list'; ul.appendChild(li); }
    input.setAttribute('aria-activedescendant', combo.active >= 0 ? 'opt-' + combo.active : '');
  }
  function comboOpen(open) {
    combo.open = open; $('target-list').hidden = !open; $('target').setAttribute('aria-expanded', String(open));
    if (open) comboRender();
  }
  function comboPick(name) {
    S.target = name; $('target').value = name; comboOpen(false); clearOutputs(); recheck();
  }
  function initCombo() {
    const input = $('target');
    input.addEventListener('focus', () => { if (!input.disabled) { input.select(); combo.active = -1; comboOpen(true); } });
    input.addEventListener('input', () => {
      combo.active = 0; comboOpen(true);
      const exact = S.metabolites.find((n) => n === input.value.trim());
      S.target = exact || ''; clearOutputs(); recheck();
      if (!exact && input.value.trim()) { /* the list is the way to choose */ }
    });
    input.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown') { e.preventDefault(); if (!combo.open) comboOpen(true); combo.active = Math.min(combo.items.length - 1, combo.active + 1); comboRender(); scrollActive(); }
      else if (e.key === 'ArrowUp') { e.preventDefault(); combo.active = Math.max(0, combo.active - 1); comboRender(); scrollActive(); }
      else if (e.key === 'Enter') { if (combo.open && combo.active >= 0 && combo.items[combo.active]) { e.preventDefault(); comboPick(combo.items[combo.active]); } }
      else if (e.key === 'Escape') { comboOpen(false); }
    });
    input.addEventListener('blur', () => {
      comboOpen(false);
      const typed = input.value.trim();
      const hit = S.metabolites.find((n) => n === typed) || S.metabolites.find((n) => n.toLowerCase() === typed.toLowerCase());
      if (hit) { S.target = hit; input.value = hit; }
      else if (typed) { S.target = ''; S.check = CBChecks.crossChecks(state()); S.check.slots.target.push({ level: 'error', text: "'" + typed + "' is not a row of the metabolite table. Choose a metabolite from the list." }); renderSlotMsgs(); updateRunState(); return; }
      recheck();
    });
  }
  function scrollActive() { const li = $('opt-' + combo.active); if (li) li.scrollIntoView({ block: 'nearest' }); }

  /* ---- run --------------------------------------------------------------------------------------------- */
  async function run() {
    if ($('btn-run').disabled) return;
    S.post = {}; renderSlotMsgs();
    $('failure').hidden = true; $('results').hidden = true;
    const c = S.check;
    S.running = true; S.runId++;
    const id = 'job' + S.runId;
    const metric = $('metric').value;
    S.runPlan = {
      metric, design: !!S.files.design, annotation: !!S.files.annotation, clusters: Number($('p-clusters').value) > 0,
      baseWasReady: S.baseReady, deseqWasReady: S.deseqReady, genes: S.infos.gene.nRows,
      contrast: c.contrast ? c.contrast.end + ' / ' + c.contrast.start : null,
    };
    ['read', 'fc', 'metric', 'cluster', 'annot', 'score', 'deseq'].forEach((k) => { if (!(k === 'deseq' && S.steps.deseq && S.steps.deseq.state !== 'pending')) delete S.steps[k]; });
    S.runStart = performance.now(); tick(); clearInterval(S.ticker); S.ticker = setInterval(tick, 500);
    $('progress').hidden = false; updateRunState(); renderProgress();
    $('progress').scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    if (!S.worker) newWorker();
    const files = []; const transfer = [];
    for (const r of ROLES) {
      if (!S.files[r]) continue;
      const buffer = await S.files[r].arrayBuffer();
      files.push({ role: r, name: S.files[r].name, buffer }); transfer.push(buffer);
    }
    if (!S.running) return;
    const params = {
      target: S.target, metric, cluster_count: Number($('p-clusters').value), lag: Number($('p-lag').value), E: Number($('p-E').value), tau: Number($('p-tau').value),
    };
    S.worker.postMessage({ type: 'run', id, files, params, wheelUrl: WHEEL }, transfer);
  }
  function stopRun() {
    S.running = false; clearInterval(S.ticker);
    $('progress').hidden = true;
    updateRunState();
  }
  function cancelRun() {
    if (!S.running) return;
    stopRun();
    newWorker();
    renderRuntime();
    warmUp();
    $('run-status').textContent = 'Run cancelled. You can run again.';
  }

  /* ---- results ----------------------------------------------------------------------------------------- */
  const R = { rows: [], cols: [], head: [], order: [], shown: 100, query: '', sort: { key: null, dir: 1 }, sel: -1, hay: null, nameIdx: new Map(), curveRow: new Map(), scatter: null };

  function onResult(m) {
    const secs = (performance.now() - S.runStart) / 1000;
    stopRun();
    if (!m.ok) { fail(m.message, m.errorType + ': ' + m.message + '\n\n' + (m.traceback || ''), false); attachWarnings(m.warnings, false); return; }
    S.result = m; S.result.secs = secs;
    $('failure').hidden = true;
    $('results').hidden = false;
    buildResults(m);
    attachWarnings(m.warnings, true);
    $('results').scrollIntoView({ block: 'start', behavior: 'smooth' });
    window.__lastRun = { secs, runMs: m.runMs };
  }
  function attachWarnings(ws, toResults) {
    const general = [];
    (ws || []).forEach((w) => {
      const item = { level: 'warn', text: w.text };
      const slot = slotOf(w.text);
      if (slot && slot !== 'target') S.post[slot] = (S.post[slot] || []).concat([item]);
      else general.push(item);
    });
    if (S.notice) { general.push({ level: 'warn', text: S.notice }); S.notice = null; }
    renderSlotMsgs();
    if (toResults) renderMsgs($('result-warnings'), general);
    else if (general.length) { const ul = $('failure-msgs'); general.forEach((g) => { const li = document.createElement('li'); li.className = 'msg msg-warn'; li.innerHTML = '<span class="kind">Warning</span><span class="body">' + esc(g.text) + '</span>'; ul.appendChild(li); }); }
  }
  function buildResults(m) {
    const metric = $('metric').value; const target = S.target;
    R.cols = m.columns; R.rows = m.rows; R.metric = metric;
    R.ci = {}; m.columns.forEach((c, i) => { R.ci[c] = i; });
    R.nameIdx = new Map(m.rows.map((r, i) => [r[0], i]));
    const nT = m.times.length;
    R.curveRow = new Map(m.seriesNames.map((n, i) => [n, i]));
    R.series = m.series; R.nT = nT; R.times = m.times; R.targetSeries = m.target;
    R.hasDesign = !!S.files.design;
    R.hay = null; R.query = ''; $('table-search').value = ''; R.sort = { key: null, dir: 1 }; R.shown = 100;
    R.order = m.rows.map((_, i) => i);

    const fcc = m.fcContrast;
    const sum = $('result-sum');
    sum.innerHTML = '';
    const a = document.createElement('p'); a.className = 'main';
    a.textContent = fmtInt(m.rows.length) + ' genes ranked by score. Metric: ' + metric + '. Target metabolite: ' + target + '.';
    const b = document.createElement('p'); b.className = 'sub';
    if (fcc) b.innerHTML = 'Fold change = log2(<code>' + esc(fcc[0]) + '</code> / <code>' + esc(fcc[1]) + '</code>), the interval in which ' + esc(target) + ' rises fastest.' + (m.lfcShrunk === false ? ' DESeq2 shrinkage was not applied.' : '');
    else b.textContent = target + ' never rises between adjacent time points, so the fold change is empty for all genes.';
    const c2 = document.createElement('p'); c2.className = 'sub'; c2.textContent = 'Finished in ' + fmtSec(S.result.secs) + (S.runPlan && !S.runPlan.baseWasReady ? ' (including loading Python).' : '.');
    sum.append(a, b, c2);

    // scatter
    const xi = R.ci[metric]; const yi = R.ci.log2FoldChange;
    const N = m.rows.length; const xs = new Float64Array(N); const ys = new Float64Array(N); const top = new Uint8Array(N);
    const gr = metric === 'Granger';
    for (let i = 0; i < N; i++) {
      const xv = m.rows[i][xi]; const yv = m.rows[i][yi];
      xs[i] = xv === null ? NaN : (gr ? -Math.log10(Math.max(xv, 1e-300)) : xv);
      ys[i] = yv === null ? NaN : yv;
      top[i] = i < 100 ? 1 : 0;
    }
    R.xs = xs; R.ys = ys;
    if (R.scatter) R.scatter.destroy();
    const xLabel = gr ? 'Granger (−log10 p)' : metric === 'DTW' ? 'DTW distance' : metric;
    R.scatter = CBCharts.scatter($('chart-scatter'), {
      x: xs, y: ys, top, xLabel, yLabel: 'log2 fold change', selected: 0, label: 'Scatter plot of ' + xLabel + ' against log2 fold change for every gene',
      topLabel: 'Top ' + Math.min(100, N),
      onPick: (i) => selectRow(i, { reveal: true, flash: true }),
      tip: (i) => m.rows[i][0] + ' (rank ' + (i + 1) + ')  ' + metric + ' ' + (m.rows[i][xi] === null ? 'n/a' : +Number(m.rows[i][xi]).toPrecision(3)) + ', log2 fold change ' + (m.rows[i][yi] === null ? 'n/a' : (+m.rows[i][yi]).toFixed(2)),
    });
    $('cap-scatter').textContent = 'Each point is a gene. Select a point to find it in the table.' + (metric === 'DTW' ? ' For DTW a smaller distance is stronger.' : '') + (fcc ? '' : ' The fold change is empty, so no points are drawn.');

    buildTableHead(); renderTable();
    selectRow(0, { reveal: false });
  }

  function buildTableHead() {
    const thead = $('table').querySelector('thead');
    const defs = [{ key: 'rank', label: 'Rank', num: true }];
    R.cols.forEach((c) => {
      const label = c === 'Name' ? 'Gene' : c === 'log2FoldChange' ? 'log2 fold change' : c;
      defs.push({ key: c, label, num: !['Name', 'Description', 'Cluster'].includes(c), idx: R.ci[c] });
    });
    R.head = defs;
    thead.innerHTML = '';
    const tr = document.createElement('tr');
    defs.forEach((d) => {
      const th = document.createElement('th'); th.scope = 'col'; if (d.num) th.className = 'num';
      th.setAttribute('aria-sort', 'none');
      const b = document.createElement('button'); b.type = 'button';
      b.innerHTML = esc(d.label) + '<span class="arrow" aria-hidden="true"></span>';
      b.addEventListener('click', () => sortBy(d.key));
      th.appendChild(b); tr.appendChild(th);
    });
    thead.appendChild(tr);
  }
  function sortBy(key) {
    if (R.sort.key === key) R.sort.dir = -R.sort.dir; else { R.sort.key = key; R.sort.dir = (key === 'Name' || key === 'Description') ? 1 : (key === 'rank' ? 1 : (key === R.metric && (R.metric === 'Granger' || R.metric === 'DTW') ? 1 : -1)); }
    applyView(); renderTable(); updateSortUI();
  }
  function updateSortUI() {
    const ths = $('table').querySelectorAll('thead th');
    R.head.forEach((d, i) => {
      const on = R.sort.key === d.key;
      ths[i].setAttribute('aria-sort', on ? (R.sort.dir === 1 ? 'ascending' : 'descending') : 'none');
      ths[i].querySelector('.arrow').textContent = on ? (R.sort.dir === 1 ? '▲' : '▼') : '';
    });
  }
  function applyView() {
    const q = R.query.trim().toLowerCase();
    let idx = R.rows.map((_, i) => i);
    if (q) {
      if (!R.hay) { const di = R.ci.Description; R.hay = R.rows.map((r) => (r[0] + ' ' + (di !== undefined && r[di] ? r[di] : '')).toLowerCase()); }
      idx = idx.filter((i) => R.hay[i].includes(q));
    }
    const k = R.sort.key;
    if (k && k !== 'rank') {
      const ci = R.ci[k]; const dir = R.sort.dir;
      idx.sort((a, b) => {
        const x = R.rows[a][ci]; const y = R.rows[b][ci];
        if (x === null && y === null) return a - b; if (x === null) return 1; if (y === null) return -1;
        const c = typeof x === 'string' ? x.localeCompare(y, 'en', { numeric: true }) : x - y;
        return c ? c * dir : a - b;
      });
    } else if (k === 'rank' && R.sort.dir === -1) idx.reverse();
    R.order = idx;
  }
  function fmtCell(d, v) {
    if (v === null || v === undefined) return '';
    if (d.key === 'Description' || d.key === 'Name' || d.key === 'Cluster') return String(v);
    if (d.key === 'log2FoldChange') return (+v).toFixed(2);
    if (d.key === 'Score') return (+v).toFixed(3);
    const a = Math.abs(v);
    return a !== 0 && a < 0.001 ? (+v).toExponential(2) : (+v).toFixed(3);
  }
  function renderTable() {
    const tb = $('table').querySelector('tbody'); const n = R.order.length; const upto = Math.min(R.shown, n);
    let html = '';
    for (let k = 0; k < upto; k++) {
      const i = R.order[k]; const r = R.rows[i];
      html += '<tr data-i="' + i + '" tabindex="' + (i === R.sel ? 0 : -1) + '" aria-selected="' + (i === R.sel) + '">';
      R.head.forEach((d) => {
        if (d.key === 'rank') html += '<td class="num">' + (i + 1) + '</td>';
        else {
          const v = fmtCell(d, r[d.idx]);
          html += '<td class="' + (d.num ? 'num' : d.key === 'Description' ? 'desc' : d.key === 'Name' ? 'name' : '') + '">' + esc(v) + '</td>';
        }
      });
      html += '</tr>';
    }
    tb.innerHTML = n ? html : '<tr><td colspan="' + R.head.length + '" class="empty">No gene matches the search.</td></tr>';
    const q = R.query.trim();
    $('table-count').textContent = q ? fmtInt(n) + ' of ' + fmtInt(R.rows.length) + ' genes match' : fmtInt(R.rows.length) + ' genes';
    $('table-foot-count').textContent = n ? 'Showing ' + fmtInt(upto) + ' of ' + fmtInt(n) : '';
    $('btn-more').hidden = upto >= n;
    $('btn-more').textContent = 'Show ' + Math.min(100, n - upto) + ' more';
  }
  function selectRow(i, opt) {
    opt = opt || {};
    R.sel = i;
    if (opt.reveal) {
      let pos = R.order.indexOf(i);
      if (pos < 0) { R.query = ''; $('table-search').value = ''; applyView(); pos = R.order.indexOf(i); }
      if (pos >= R.shown) R.shown = Math.ceil((pos + 1) / 100) * 100;
      renderTable();
    } else {
      $('table').querySelectorAll('tbody tr').forEach((tr) => { const on = Number(tr.dataset.i) === i; tr.setAttribute('aria-selected', String(on)); tr.tabIndex = on ? 0 : -1; });
    }
    if (opt.reveal) {
      const tr = $('table').querySelector('tbody tr[data-i="' + i + '"]');
      if (tr) { tr.scrollIntoView({ block: 'center', behavior: 'smooth' }); if (opt.flash) { tr.classList.remove('flash'); void tr.offsetWidth; tr.classList.add('flash'); } }
    }
    drawCurve(i);
    if (R.scatter) R.scatter.select(i);
  }
  function drawCurve(i) {
    const name = R.rows[i][0]; const row = R.curveRow.get(name);
    const vals = row === undefined ? R.times.map(() => NaN) : Array.from(R.series.subarray(row * R.nT, (row + 1) * R.nT));
    CBCharts.lineChart($('chart-curve'), {
      times: R.times, xLabel: 'Time point', yLabel: 'Scaled level (0–1)', label: 'Time course of ' + name + ' and ' + S.target,
      series: [
        { name: S.target, values: R.targetSeries.map((v) => (v === null ? NaN : v)), color: getComputedStyle(document.documentElement).getPropertyValue('--c-target').trim(), shape: 'square' },
        { name, values: vals, color: getComputedStyle(document.documentElement).getPropertyValue('--c-gene').trim(), shape: 'circle' },
      ],
    });
    $('cap-curve').textContent = name + ' (rank ' + (i + 1) + ') and ' + S.target + ', each scaled to 0–1.' + (R.hasDesign ? ' Points are group means.' : '') + ' Select a table row or a point in the scatter plot to change the gene.';
  }
  function initTable() {
    const tb = $('table').querySelector('tbody');
    tb.addEventListener('click', (e) => { const tr = e.target.closest('tr[data-i]'); if (tr) selectRow(Number(tr.dataset.i), {}); });
    tb.addEventListener('keydown', (e) => {
      const tr = e.target.closest('tr[data-i]'); if (!tr) return;
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); selectRow(Number(tr.dataset.i), {}); }
      else if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
        e.preventDefault();
        const nx = e.key === 'ArrowDown' ? tr.nextElementSibling : tr.previousElementSibling;
        if (nx && nx.dataset.i) { nx.focus(); selectRow(Number(nx.dataset.i), {}); }
      }
    });
    let t = null;
    $('table-search').addEventListener('input', (e) => {
      clearTimeout(t);
      t = setTimeout(() => { R.query = e.target.value; R.shown = 100; applyView(); renderTable(); }, 150);
    });
    $('btn-more').addEventListener('click', () => { R.shown += 100; renderTable(); });
    $('btn-download').addEventListener('click', () => {
      const blob = new Blob([S.result.csv], { type: 'text/csv;charset=utf-8' });
      const a = document.createElement('a'); a.href = URL.createObjectURL(blob);
      a.download = 'catbridge_' + S.target.replace(/[^A-Za-z0-9._-]+/g, '_') + '_' + R.metric + '.csv';
      document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    });
  }

  /* ---- example, reset, wiring ---------------------------------------------------------------------------- */
  async function loadExample() {
    const btn = $('btn-example'); btn.disabled = true; btn.textContent = 'Loading example...';
    try {
      const got = {};
      await Promise.all(ROLES.map(async (r) => {
        const resp = await fetch('example/' + EXAMPLE[r]);
        if (!resp.ok) throw new Error('example/' + EXAMPLE[r] + ' (HTTP ' + resp.status + ')');
        got[r] = new File([await resp.blob()], EXAMPLE[r], { type: 'text/tab-separated-values' });
      }));
      ROLES.forEach((r) => removeFile(r));
      for (const r of ROLES) await setFile(r, got[r]);
      if (S.metabolites.includes(EXAMPLE_TARGET)) comboPick(EXAMPLE_TARGET);
    } catch (e) {
      fail('The example data could not be loaded: ' + e.message + '. Check your internet connection and try again.', '', true);
    } finally { btn.disabled = false; btn.textContent = 'Load example data'; }
  }
  function resetAll() {
    if (S.running) cancelRun();
    ROLES.forEach((r) => removeFile(r));
    S.target = ''; $('target').value = '';
    $('p-clusters').value = 8; $('p-lag').value = 1; $('p-E').value = 3; $('p-tau').value = 1; $('metric').value = 'CCM'; metricHint();
    $('results').hidden = true; $('failure').hidden = true; S.post = {};
    recheck();
  }
  function metricHint() { $('metric-hint').textContent = METRIC_HINT[$('metric').value]; }

  function init() {
    if (!window.Worker || !window.WebAssembly) {
      $('run-status').textContent = 'This browser does not support WebAssembly or web workers. Use a current version of Chrome, Edge, Firefox or Safari.';
    }
    ROLES.forEach((r) => {
      const input = $('file-' + r);
      input.addEventListener('change', () => { if (input.files[0]) setFile(r, input.files[0]); });
      const drop = slotEl(r).querySelector('.drop');
      ['dragenter', 'dragover'].forEach((ev) => drop.addEventListener(ev, () => drop.classList.add('over')));
      ['dragleave', 'drop'].forEach((ev) => drop.addEventListener(ev, () => drop.classList.remove('over')));
      slotEl(r).querySelector('[data-remove]').addEventListener('click', () => { removeFile(r); input.focus(); });
    });
    $('form').addEventListener('submit', (e) => { e.preventDefault(); run(); });
    $('btn-cancel').addEventListener('click', cancelRun);
    $('btn-example').addEventListener('click', loadExample);
    $('btn-reset').addEventListener('click', resetAll);
    $('metric').addEventListener('change', () => { metricHint(); clearOutputs(); recheck(); });
    ['p-clusters', 'p-lag', 'p-E', 'p-tau'].forEach((id) => $(id).addEventListener('input', () => { clearOutputs(); recheck(); }));
    initCombo(); initTable(); metricHint(); recheck(); renderRuntime();
    window.__cb = { S, R, setFile, removeFile, resetAll, comboPick, run, cancelRun };
  }
  init();
})();
