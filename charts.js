/* Two small charts without a library: a time-course line chart (SVG) and a scatter plot (canvas, tens of thousands of points). */
(function (root) {
  'use strict';
  const NS = 'http://www.w3.org/2000/svg';
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  function el(tag, attrs, text) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function niceTicks(lo, hi, n) {
    if (!(hi > lo)) { hi = lo + 1; }
    const raw = (hi - lo) / Math.max(1, n);
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    const f = raw / p;
    const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
    const start = Math.ceil(lo / step - 1e-9) * step;
    const ticks = [];
    for (let v = start; v <= hi + step * 1e-6; v += step) ticks.push(+v.toFixed(12));
    return { ticks, step };
  }
  function fmt(v, step) {
    if (v === 0) return '0';
    const d = step >= 1 ? 0 : Math.min(6, Math.ceil(-Math.log10(step)));
    if (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0)) return v.toExponential(0);
    return v.toFixed(d);
  }
  function short(s, n) { s = String(s); return s.length > n ? s.slice(0, n - 1) + '…' : s; }
  function scale01(values) {
    const f = values.filter(Number.isFinite);
    if (!f.length) return values.map(() => NaN);
    const lo = Math.min.apply(null, f); const hi = Math.max.apply(null, f);
    return values.map((v) => (Number.isFinite(v) ? (hi > lo ? (v - lo) / (hi - lo) : 0.5) : NaN));
  }
  const num = (v) => (Number.isFinite(v) ? (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0) ? v.toPrecision(3) : String(+v.toPrecision(4))) : 'missing');

  /* ---- line chart --------------------------------------------------------------------------------------- */
  // spec: {times:[], series:[{name, values:[], color, shape:'circle'|'square'}], yLabel, xLabel}
  function lineChart(box, spec) {
    box.__spec = spec;
    if (!box.__ro) { box.__ro = new ResizeObserver(() => { if (box.__w && Math.abs(box.clientWidth - box.__w) > 1) lineChart(box, box.__spec); }); box.__ro.observe(box); }
    box.__w = box.clientWidth;
    box.innerHTML = '';
    const W = Math.max(280, box.clientWidth); const H = W < 420 ? 300 : 320;
    const m = { l: 52, r: 30, t: 44, b: 56 };
    const iw = W - m.l - m.r; const ih = H - m.t - m.b;
    const svg = el('svg', { viewBox: '0 0 ' + W + ' ' + H, role: 'img', class: 'chart-svg', 'aria-label': spec.label || 'Time course' });
    const n = spec.times.length;
    const x = (i) => m.l + (n === 1 ? iw / 2 : (iw * i) / (n - 1));
    const y = (v) => m.t + ih * (1 - v);
    // grid + y axis
    [0, 0.25, 0.5, 0.75, 1].forEach((t) => {
      svg.appendChild(el('line', { x1: m.l, x2: m.l + iw, y1: y(t), y2: y(t), stroke: css('--line'), 'stroke-width': 1 }));
      svg.appendChild(el('text', { x: m.l - 8, y: y(t) + 4, 'text-anchor': 'end' }, String(t)));
    });
    const yt = el('text', { class: 'axis-title', transform: 'translate(14,' + (m.t + ih / 2) + ') rotate(-90)', 'text-anchor': 'middle' }, spec.yLabel);
    svg.appendChild(yt);
    // x labels (thin out when crowded, rotate when long)
    const maxLen = Math.max.apply(null, spec.times.map((t) => String(t).length));
    const slot = n > 1 ? iw / (n - 1) : iw;
    const rotate = maxLen * 7.5 > slot * 0.9;
    const every = rotate ? Math.max(1, Math.ceil(14 / slot)) : Math.max(1, Math.ceil((maxLen * 6.5) / slot));
    spec.times.forEach((t, i) => {
      svg.appendChild(el('line', { x1: x(i), x2: x(i), y1: m.t + ih, y2: m.t + ih + 4, stroke: css('--ink-3') }));
      if (i % every) return;
      if (rotate) svg.appendChild(el('text', { x: x(i), y: m.t + ih + 16, transform: 'rotate(-35 ' + x(i) + ' ' + (m.t + ih + 16) + ')', 'text-anchor': 'end' }, short(t, 14)));
      else svg.appendChild(el('text', { x: x(i), y: m.t + ih + 18, 'text-anchor': 'middle' }, short(t, 14)));
    });
    svg.appendChild(el('line', { x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih, stroke: css('--ink-3') }));
    svg.appendChild(el('text', { class: 'axis-title', x: m.l + iw / 2, y: H - 6, 'text-anchor': 'middle' }, spec.xLabel));
    // series
    spec.series.forEach((s, si) => {
      const sc = scale01(s.values);
      let d = ''; let pen = false;
      sc.forEach((v, i) => { if (!Number.isFinite(v)) { pen = false; return; } d += (pen ? 'L' : 'M') + x(i).toFixed(1) + ' ' + y(v).toFixed(1); pen = true; });
      svg.appendChild(el('path', { d, fill: 'none', stroke: s.color, 'stroke-width': 2.5, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
      sc.forEach((v, i) => {
        if (!Number.isFinite(v)) return;
        const g = el('g', {});
        const shape = s.shape === 'square'
          ? el('rect', { x: x(i) - 4.5, y: y(v) - 4.5, width: 9, height: 9, fill: s.color, stroke: '#fff', 'stroke-width': 1.5 })
          : el('circle', { cx: x(i), cy: y(v), r: 5, fill: s.color, stroke: '#fff', 'stroke-width': 1.5 });
        g.appendChild(shape);
        g.appendChild(el('title', {}, s.name + ' at ' + spec.times[i] + ': ' + num(s.values[i])));
        svg.appendChild(g);
      });
      // legend
      const lx = m.l + si * Math.min(iw / 2, 220); const ly = 16;
      if (s.shape === 'square') svg.appendChild(el('rect', { x: lx, y: ly - 5, width: 10, height: 10, fill: s.color }));
      else svg.appendChild(el('circle', { cx: lx + 5, cy: ly, r: 5, fill: s.color }));
      svg.appendChild(el('line', { x1: lx - 4, x2: lx + 14, y1: ly, y2: ly, stroke: s.color, 'stroke-width': 2.5 }));
      svg.appendChild(el('text', { x: lx + 20, y: ly + 4 }, short(s.name, Math.max(10, Math.floor(Math.min(iw / 2, 220) / 7.5)))));
    });
    box.appendChild(svg);
  }

  /* ---- scatter (canvas) ---------------------------------------------------------------------------------- */
  // spec: {x:Float64Array|[], y:[], top:Uint8Array (1 = highlighted group), xLabel, yLabel, onPick(index), selected, tip(i)}
  function scatter(box, spec) {
    box.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', spec.label || 'Scatter plot');
    const tip = document.createElement('div'); tip.className = 'chart-tip'; tip.hidden = true;
    box.appendChild(canvas); box.appendChild(tip);
    const state = { selected: spec.selected == null ? -1 : spec.selected, W: 0, H: 0, m: null, sx: null, sy: null, hover: -1 };
    const N = spec.x.length;

    let xlo = Infinity; let xhi = -Infinity; let ylo = Infinity; let yhi = -Infinity;
    for (let i = 0; i < N; i++) {
      const a = spec.x[i]; const b = spec.y[i];
      if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
      if (a < xlo) xlo = a; if (a > xhi) xhi = a; if (b < ylo) ylo = b; if (b > yhi) yhi = b;
    }
    if (!(xhi >= xlo)) { xlo = 0; xhi = 1; ylo = 0; yhi = 1; }
    const padx = (xhi - xlo || 1) * 0.04; const pady = (yhi - ylo || 1) * 0.06;
    xlo -= padx; xhi += padx; ylo = Math.min(ylo - pady, 0 - pady * 0.2); yhi += pady;

    function draw() {
      const dpr = window.devicePixelRatio || 1;
      const W = Math.max(280, box.clientWidth); const H = W < 420 ? 300 : 320;
      state.W = W; state.H = H;
      canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
      canvas.style.height = H + 'px';
      const c = canvas.getContext('2d'); c.setTransform(dpr, 0, 0, dpr, 0, 0);
      c.clearRect(0, 0, W, H);
      const m = { l: 52, r: 16, t: 36, b: 48 }; state.m = m;
      const iw = W - m.l - m.r; const ih = H - m.t - m.b;
      const sx = (v) => m.l + ((v - xlo) / (xhi - xlo)) * iw;
      const sy = (v) => m.t + ih - ((v - ylo) / (yhi - ylo)) * ih;
      state.sx = sx; state.sy = sy;
      c.font = '12px ' + css('--font'); c.textBaseline = 'middle';
      // grid
      const xt = niceTicks(xlo, xhi, Math.max(3, Math.floor(iw / 90)));
      const yt = niceTicks(ylo, yhi, 5);
      c.strokeStyle = css('--line'); c.lineWidth = 1; c.fillStyle = css('--ink-2');
      c.textAlign = 'right';
      yt.ticks.forEach((t) => { const py = Math.round(sy(t)) + 0.5; c.beginPath(); c.moveTo(m.l, py); c.lineTo(m.l + iw, py); c.stroke(); c.fillText(fmt(t, yt.step), m.l - 8, py); });
      c.textAlign = 'center';
      xt.ticks.forEach((t) => { const px = Math.round(sx(t)) + 0.5; c.beginPath(); c.moveTo(px, m.t); c.lineTo(px, m.t + ih); c.stroke(); c.fillText(fmt(t, xt.step), px, m.t + ih + 14); });
      // zero line of the fold change
      if (ylo < 0 && yhi > 0) { c.strokeStyle = css('--ink-3'); c.beginPath(); c.moveTo(m.l, Math.round(sy(0)) + 0.5); c.lineTo(m.l + iw, Math.round(sy(0)) + 0.5); c.stroke(); }
      c.strokeStyle = css('--ink-3'); c.strokeRect(m.l + 0.5, m.t + 0.5, iw, ih);
      // points: others first, then top group, then selection
      c.fillStyle = 'rgba(139,145,168,0.38)';
      for (let i = 0; i < N; i++) {
        if (spec.top[i]) continue;
        const a = spec.x[i]; const b = spec.y[i];
        if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
        c.beginPath(); c.arc(sx(a), sy(b), 2.2, 0, 6.2832); c.fill();
      }
      c.fillStyle = 'rgba(74,84,178,0.85)'; c.strokeStyle = '#fff'; c.lineWidth = 0.8;
      for (let i = 0; i < N; i++) {
        if (!spec.top[i]) continue;
        const a = spec.x[i]; const b = spec.y[i];
        if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
        c.beginPath(); c.arc(sx(a), sy(b), 3.4, 0, 6.2832); c.fill(); c.stroke();
      }
      const s = state.selected;
      if (s >= 0 && Number.isFinite(spec.x[s]) && Number.isFinite(spec.y[s])) {
        c.beginPath(); c.arc(sx(spec.x[s]), sy(spec.y[s]), 6.5, 0, 6.2832);
        c.fillStyle = css('--c-target'); c.fill(); c.lineWidth = 2; c.strokeStyle = '#fff'; c.stroke();
      }
      // axis titles
      c.fillStyle = css('--ink'); c.font = '500 12px ' + css('--font'); c.textAlign = 'center';
      c.fillText(spec.xLabel, m.l + iw / 2, H - 10);
      c.save(); c.translate(14, m.t + ih / 2); c.rotate(-Math.PI / 2); c.fillText(spec.yLabel, 0, 0); c.restore();
      // legend
      c.font = '12px ' + css('--font'); c.textAlign = 'left'; let lx = m.l;
      [[css('--primary'), spec.topLabel || 'Top 100'], ['rgba(139,145,168,0.7)', 'Other genes'], [css('--c-target'), 'Selected']].forEach(([col, t]) => {
        c.fillStyle = col; c.beginPath(); c.arc(lx + 5, 16, 5, 0, 6.2832); c.fill();
        c.fillStyle = css('--ink-2'); c.fillText(t, lx + 14, 16); lx += 14 + c.measureText(t).width + 16;
      });
    }

    function nearest(px, py) {
      let best = -1; let bd = 12 * 12;
      for (let i = 0; i < N; i++) {
        const a = spec.x[i]; const b = spec.y[i];
        if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
        const dx = state.sx(a) - px; const dy = state.sy(b) - py;
        const d = dx * dx + dy * dy;
        if (d < bd || (d === bd && best >= 0 && spec.top[i] && !spec.top[best])) { bd = d; best = i; }
      }
      return best;
    }
    const pos = (ev) => { const r = canvas.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    canvas.addEventListener('pointermove', (ev) => {
      const [px, py] = pos(ev); const i = nearest(px, py);
      state.hover = i;
      canvas.style.cursor = i >= 0 ? 'pointer' : 'crosshair';
      if (i < 0) { tip.hidden = true; return; }
      tip.hidden = false; tip.textContent = spec.tip(i);
      const bw = box.clientWidth; tip.style.left = Math.min(bw - 60, Math.max(60, state.sx(spec.x[i]))) + 'px'; tip.style.top = state.sy(spec.y[i]) + 'px';
    });
    canvas.addEventListener('pointerleave', () => { tip.hidden = true; });
    canvas.addEventListener('click', (ev) => {
      const [px, py] = pos(ev); const i = nearest(px, py);
      if (i >= 0) spec.onPick(i);
    });
    const ro = new ResizeObserver(() => { if (Math.abs(box.clientWidth - state.W) > 1) draw(); });
    ro.observe(box);
    draw();
    return {
      select(i) { state.selected = i; draw(); },
      destroy() { ro.disconnect(); },
    };
  }

  root.CBCharts = { lineChart, scatter };
})(self);
