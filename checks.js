/* CAT Bridge input checks (run in the browser before the analysis).
 * Per file: extension, empty file, binary content, delimiter, header, row length, first column, duplicate names,
 *           non-numeric / infinite / negative values (all rows, first 40 MB at most).
 * Across files: matching columns, design, annotation, target metabolite, advanced parameters.
 * A problem is {level: 'error' | 'warn', text}. Errors block the run, warnings do not.
 */
(function (root) {
  'use strict';

  const MAX_BYTES = 40 * 1024 * 1024;
  const NA = new Set(['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN',
    '<NA>', 'N/A', 'NA', 'NULL', 'NaN', 'None', 'n/a', 'nan', 'null']);
  const NUM = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;
  const INF = /^[+-]?inf(inity)?$/i;
  const DELIMS = { '\t': 'tab', ',': 'comma', ';': 'semicolon' };
  const ROLE_NAME = { gene: 'gene', metabolite: 'metabolite', design: 'design', annotation: 'annotation' };

  const err = (text) => ({ level: 'error', text });
  const warn = (text) => ({ level: 'warn', text });
  const list = (a, n) => a.slice(0, n || 5).map((x) => "'" + x + "'").join(', ') + (a.length > (n || 5) ? ', ...' : '');
  const plural = (n, one, many) => n + ' ' + (n === 1 ? one : many || one + 's');
  const fmtNum = (v) => (Number.isInteger(v) ? String(v) : String(+v.toPrecision(6)));

  function splitLine(line, delim) {
    if (line.indexOf('"') < 0) return line.split(delim);
    const out = []; let cur = ''; let q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
      else if (c === '"') q = true;
      else if (c === delim) { out.push(cur); cur = ''; }
      else cur += c;
    }
    out.push(cur);
    return out;
  }

  function editDistance(a, b) {
    if (Math.abs(a.length - b.length) > 3) return 9;
    let prev = Array.from({ length: b.length + 1 }, (_, i) => i);
    for (let i = 1; i <= a.length; i++) {
      const cur = [i];
      for (let j = 1; j <= b.length; j++) cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
      prev = cur;
    }
    return prev[b.length];
  }
  function similar(name, pool) {
    let best = null; let bd = 3;
    for (const p of pool) { const d = editDistance(name.toLowerCase(), p.toLowerCase()); if (d < bd) { bd = d; best = p; } }
    return best;
  }

  function extensionOf(name) {
    const m = /\.([A-Za-z0-9]+)$/.exec(name);
    return m ? m[1].toLowerCase() : '';
  }

  /* Check one file. Returns {role, name, size, ext, delim, header, names, rows, issues, fatal}.
   * `rows` (parsed numeric rows) is kept for the metabolite table only. */
  async function analyzeFile(role, file) {
    const info = { role, name: file.name, size: file.size, ext: extensionOf(file.name), delim: null, header: [], names: [],
      nRows: 0, issues: [], fatal: false, rowMap: null, valuesNegative: 0, nonInteger: false, desc: null };
    const bad = (t) => { info.issues.push(err(t)); info.fatal = true; return info; };
    const ext = info.ext;

    if (ext === 'xlsx' || ext === 'xls') {
      return bad('Excel files are not supported. In Excel choose File > Save As > CSV (Comma delimited), then select the .csv file.');
    }
    if (['gz', 'zip', 'bz2', 'xz', 'zst'].includes(ext)) {
      return bad('Compressed files are not supported here. Unzip the file and select the .tsv, .csv or .txt file.');
    }
    if (!['tsv', 'csv', 'txt'].includes(ext)) {
      return bad(ext ? "'." + ext + "' files are not supported. Use a .tsv, .csv or .txt file." : 'The file has no extension. Use a .tsv, .csv or .txt file.');
    }
    if (file.size === 0) return bad('The file is empty (0 bytes). Select a file that contains a header line and data rows.');

    const truncated = file.size > MAX_BYTES;
    const buf = await file.slice(0, MAX_BYTES).arrayBuffer();
    const bytes = new Uint8Array(buf);
    if (bytes[0] === 0x50 && bytes[1] === 0x4b && (bytes[2] === 3 || bytes[2] === 5)) {
      return bad('The content is a zip archive (an Excel file renamed to .' + ext + '?). Save it as CSV or TSV text and select that file.');
    }
    if (bytes.subarray(0, 4096).indexOf(0) >= 0) {
      return bad('The file contains binary data, not text. Save it as a text table (.tsv or .csv) and select that file.');
    }
    let text = new TextDecoder('utf-8').decode(bytes);
    if (text.charCodeAt(0) === 0xfeff) text = text.slice(1);
    let lines = text.split(/\r\n|\n|\r/);
    if (truncated) lines.pop();                       // last line is cut
    lines = lines.filter((l) => l.trim() !== '');
    if (!lines.length) return bad('The file contains no text. Select a file with a header line and data rows.');
    if (truncated) info.issues.push(warn('The file is larger than 40 MB; only the first 40 MB were checked.'));

    // delimiter
    const head = lines[0];
    const expected = ext === 'csv' ? ',' : '\t';
    const counts = {};
    Object.keys(DELIMS).forEach((d) => { counts[d] = head.split(d).length - 1; });
    let delim = expected;
    if (!counts[expected]) {
      const alt = Object.keys(DELIMS).filter((d) => d !== expected).sort((a, b) => counts[b] - counts[a])[0];
      if (counts[alt] > 0) {
        delim = alt;
        info.issues.push(warn("The extension '." + ext + "' means " + DELIMS[expected] + '-separated, but the file is ' + DELIMS[alt] +
          '-separated. It will be read as ' + DELIMS[alt] + '-separated.'));
      } else {
        return bad('The header has only one column, so no separator was found. Separate columns with tabs (.tsv, .txt) or commas (.csv).');
      }
    }
    info.delim = delim;
    const header = splitLine(head, delim).map((s) => s.trim());
    info.header = header;
    if (header.length < 2) return bad('The header has only one column. Separate columns with tabs (.tsv, .txt) or commas (.csv).');
    if (lines.length < 2) return bad('The file has a header but no data rows. Add at least one row below the header.');

    // duplicate column names
    const seen = new Map(); const dupCols = [];
    header.slice(1).forEach((h) => { if (h === '') return; if (seen.has(h)) { if (!dupCols.includes(h)) dupCols.push(h); } seen.set(h, 1); });
    if (dupCols.length) {
      info.issues.push(err('Duplicated column names: ' + list(dupCols, 5) + '. Give every ' + (role === 'design' ? 'column' : 'sample or time point') + ' a unique name.'));
      info.fatal = true;
    }

    if (role === 'design') return analyzeDesign(info, lines, delim);

    const nCol = header.length;
    const numeric = role !== 'annotation';
    const keepRows = role === 'metabolite';
    const names = []; const rows = []; const nameSet = new Map(); const dupNames = [];
    const nonNum = []; const infCells = []; let nNonNum = 0; let nInf = 0; let neg = 0; let negEx = null; let nonInt = false;
    const longRows = []; const shortRows = []; let blankNames = []; let numericNames = 0; let comment = 0;
    for (let i = 1; i < lines.length; i++) {
      const cells = splitLine(lines[i], delim);
      const nm = (cells[0] || '').trim();
      if (role === 'annotation' && nm.startsWith('#')) { comment++; continue; }
      if (cells.length > nCol) { if (longRows.length < 3) longRows.push([i, cells.length]); else longRows.push(null); }
      else if (cells.length < nCol) shortRows.push([i, cells.length]);
      names.push(nm);
      if (nm === '') blankNames.push(i);
      else if (NUM.test(nm)) numericNames++;
      if (nm !== '') { if (nameSet.has(nm)) { if (dupNames.length < 6 && !dupNames.includes(nm)) dupNames.push(nm); } else nameSet.set(nm, i); }
      if (numeric) {
        const vals = keepRows ? new Array(nCol - 1) : null;
        for (let c = 1; c < nCol; c++) {
          const raw = c < cells.length ? cells[c].trim() : '';
          let v = NaN;
          if (!NA.has(raw)) {
            if (NUM.test(raw)) {
              v = parseFloat(raw);
              if (v < 0) { neg++; if (!negEx) negEx = "row '" + nm + "', column '" + header[c] + "': " + raw; }
              if (!Number.isInteger(v)) nonInt = true;
            } else if (INF.test(raw)) { nInf++; if (infCells.length < 3) infCells.push("row '" + nm + "', column '" + header[c] + "'"); }
            else { nNonNum++; if (nonNum.length < 3) nonNum.push("row '" + nm + "', column '" + header[c] + "': '" + raw + "'"); }
          }
          if (keepRows) vals[c - 1] = v;
        }
        if (keepRows) rows.push(vals);
      }
    }
    info.nRows = names.length; info.names = names; info.nonInteger = nonInt; info.valuesNegative = neg; info.negExample = negEx;
    if (keepRows) { info.rowMap = new Map(); names.forEach((n, i) => { if (!info.rowMap.has(n)) info.rowMap.set(n, rows[i]); }); }
    if (!names.length) return bad('The file has no data rows. Add at least one row below the header.');

    const what = role === 'gene' ? 'gene' : role === 'metabolite' ? 'metabolite' : 'gene';
    if (longRows.length) {
      const ex = longRows.filter(Boolean).map((r) => 'data row ' + r[0] + ' has ' + r[1]).join('; ');
      info.issues.push(err('Some rows have more values than the header has columns (' + nCol + '): ' + ex + '. Check separators, quotes and decimal commas (use a point: 1.5).'));
      info.fatal = true;
    }
    if (shortRows.length) {
      info.issues.push(warn(plural(shortRows.length, 'row has', 'rows have') + ' fewer values than the header (first: data row ' + shortRows[0][0] +
        '). Missing cells are read as empty.'));
    }
    if (blankNames.length) {
      info.issues.push(err(plural(blankNames.length, 'row has', 'rows have') + ' no name in the first column (data row ' + blankNames[0] + '). Fill in the ' + what + ' name or delete the row.'));
      info.fatal = true;
    }
    if (numericNames / Math.max(1, names.length) > 0.5) {
      info.issues.push(err('The first column holds numbers (e.g. ' + names.slice(0, 3).join(', ') + '), not names. Add a first column with the ' + what + ' names.'));
      info.fatal = true;
    }
    if (role === 'gene' && dupNames.length) {
      info.issues.push(err('Duplicated gene names: ' + list(dupNames, 5) + '. Every gene name must appear once.'));
      info.fatal = true;
    }
    if (numeric) {
      if (nNonNum) {
        info.issues.push(err(plural(nNonNum, 'cell is', 'cells are') + ' not a number, e.g. ' + nonNum.join('; ') + '. Use plain numbers, leave missing values empty (or NA), and use a point as the decimal mark.'));
        info.fatal = true;
      }
      if (nInf) {
        info.issues.push(err(plural(nInf, 'cell is', 'cells are') + ' infinite, e.g. ' + infCells.join('; ') + '. Replace with a number or leave the cell empty.'));
        info.fatal = true;
      }
      if (neg && role === 'metabolite') {
        info.issues.push(warn(plural(neg, 'value is', 'values are') + ' negative (e.g. ' + negEx + '). Check that the table holds abundances or log-ratios you intend to use.'));
      }
      // negative gene values are judged together with the design table (see crossChecks)
    }
    if (comment) info.commentRows = comment;
    return info;
  }

  function analyzeDesign(info, lines, delim) {
    const header = info.header;
    const gi = header.findIndex((h, i) => i >= 0 && h.toLowerCase() === 'group');
    if (gi < 1) {
      info.issues.push(err("The design table needs a column named 'group' (the time point of each sample); its columns are " + list(header, 6) +
        ". The first column must hold the sample names."));
      info.fatal = true;
    }
    const samples = []; const groups = []; const dup = []; const noGroup = [];
    const seen = new Set();
    for (let i = 1; i < lines.length; i++) {
      const cells = splitLine(lines[i], delim).map((s) => s.trim());
      if (cells[0] === '') continue;
      samples.push(cells[0]);
      const g = gi >= 1 ? (cells[gi] || '') : '';
      groups.push(g);
      if (seen.has(cells[0])) { if (!dup.includes(cells[0])) dup.push(cells[0]); }
      seen.add(cells[0]);
      if (gi >= 1 && g === '') noGroup.push(cells[0]);
    }
    info.samples = samples; info.groups = groups; info.nRows = samples.length;
    if (!samples.length) { info.issues.push(err('The design table has no sample rows.')); info.fatal = true; }
    if (dup.length) { info.issues.push(err('The design table lists the same sample more than once: ' + list(dup, 5) + '. Keep one row per sample.')); info.fatal = true; }
    if (noGroup.length) { info.issues.push(err('Samples without a group: ' + list(noGroup, 5) + '. Fill in the group (time point) of every sample.')); info.fatal = true; }
    return info;
  }

  /* ---- target series -------------------------------------------------------------------------------------- */
  function groupMeans(columns, values, design) {
    // values: numbers per column; returns {labels, means}. Groups in order of first appearance in `columns`.
    const gOf = new Map(); design.samples.forEach((s, i) => gOf.set(s, design.groups[i]));
    const order = []; const acc = new Map();
    columns.forEach((c, i) => {
      const g = gOf.get(c); if (g === undefined) return;
      if (!acc.has(g)) { acc.set(g, { s: 0, n: 0 }); order.push(g); }
      if (Number.isFinite(values[i])) { const a = acc.get(g); a.s += values[i]; a.n++; }
    });
    return { labels: order, means: order.map((g) => { const a = acc.get(g); return a.n ? a.s / a.n : NaN; }) };
  }

  function steepestRise(labels, series) {
    let best = -Infinity; let bi = -1;
    for (let i = 0; i + 1 < series.length; i++) {
      const d = series[i + 1] - series[i];
      if (Number.isFinite(d) && d > best) { best = d; bi = i; }
    }
    return bi >= 0 && best > 0 ? { end: labels[bi + 1], start: labels[bi] } : null;
  }

  /* ---- cross-file checks ---------------------------------------------------------------------------------- */
  // state: {infos:{gene,metabolite,design,annotation}, target, metric, params:{lag,E,tau,clusters}}
  // returns {slots:{gene:[], metabolite:[], design:[], annotation:[], target:[], params:[]}, todo:[], contrast, timePoints}
  function crossChecks(state) {
    const out = { slots: { gene: [], metabolite: [], design: [], annotation: [], target: [], params: [] }, todo: [], contrast: null, timePoints: null, targetSeries: null };
    const I = state.infos; const g = I.gene; const m = I.metabolite; const d = I.design; const a = I.annotation;
    const ok = (x) => x && !x.fatal;
    if (!g) out.todo.push('add the gene table');
    if (!m) out.todo.push('add the metabolite table');
    const gOK = ok(g); const mOK = ok(m); const dOK = ok(d);

    // first column looks like a sample name
    [['gene', g, m], ['metabolite', m, g]].forEach(([role, t, other]) => {
      if (!ok(t)) return;
      const known = new Set(dOK ? d.samples : (ok(other) ? other.header.slice(1) : []));
      if (t.header[0] && known.has(t.header[0])) {
        out.slots[role].push(err("The first column is named '" + t.header[0] + "', which is a sample name: the table seems to have no name column. Add a first column with the " + role + ' names.'));
        t.fatal = true;
      }
    });

    // time points / samples
    let columnsG = gOK ? g.header.slice(1) : null; let columnsM = mOK ? m.header.slice(1) : null;
    if (gOK && mOK && !dOK) {
      if (d && !dOK) { /* the design table has its own error */ }
      const cg = columnsG; const cm = columnsM;
      if (cg.length !== cm.length || cg.some((x, i) => x !== cm[i])) {
        let k = 0; while (k < cg.length && k < cm.length && cg[k] === cm[k]) k++;
        const onlyG = cg.filter((x) => !cm.includes(x)); const onlyM = cm.filter((x) => !cg.includes(x));
        let t;
        if (!onlyG.length && !onlyM.length && cg.length === cm.length) {
          t = 'The gene and metabolite tables have the same columns in a different order (first difference at column ' + (k + 1) + ": '" + cg[k] + "' vs '" + cm[k] + "'). Reorder the columns so that both tables match.";
        } else {
          t = 'The columns of the gene and metabolite tables do not match (first difference at column ' + (k + 1) + ": '" + (cg[k] === undefined ? '(none)' : cg[k]) + "' in the gene table, '" + (cm[k] === undefined ? '(none)' : cm[k]) +
            "' in the metabolite table" + (cg.length === cm.length ? '' : '; ' + cg.length + ' vs ' + cm.length + ' columns') + '). ';
          if (onlyG.length) t += 'Only in gene: ' + list(onlyG, 4) + '. ';
          if (onlyM.length) t += 'Only in metabolite: ' + list(onlyM, 4) + '. ';
          t += 'Both tables need the same samples or time points in the same order.';
        }
        out.slots.metabolite.push(err(t));
        out.slots.gene.push(err('Columns differ from the metabolite table (first difference at column ' + (k + 1) + '). See the metabolite table.'));
      }
    }
    if (dOK) {
      const known = new Set(d.samples);
      [['gene', g, columnsG], ['metabolite', m, columnsM]].forEach(([role, t, cols]) => {
        if (!ok(t)) return;
        const miss = cols.filter((c) => !known.has(c));
        if (miss.length) {
          const hints = miss.slice(0, 3).map((x) => { const c = similar(x, d.samples.filter((y) => !cols.includes(y))); return c ? "'" + x + "' vs '" + c + "'?" : null; }).filter(Boolean);
          out.slots.design.push(err(plural(miss.length, 'column', 'columns') + ' of the ' + role + " table " + (miss.length === 1 ? 'is' : 'are') + ' not in the first column of the design table: ' + list(miss, 4) +
            (hints.length ? ' (similar design samples: ' + hints.join(', ') + ')' : '') + '. The design table needs one row for every sample, with the same name as the column header.'));
        }
      });
      const used = new Set([].concat(columnsG || [], columnsM || []));
      const unused = d.samples.filter((s) => !used.has(s));
      if (unused.length && (gOK || mOK) && !out.slots.design.length) {
        out.slots.design.push(warn(plural(unused.length, 'design sample is', 'design samples are') + ' not a column of the gene or metabolite table and will be ignored: ' + list(unused, 4) + '.'));
      }
      // order of the groups and replicates
      if (gOK && mOK && !out.slots.design.some((x) => x.level === 'error')) {
        const gm = groupMeans(columnsG, new Array(columnsG.length).fill(1), d).labels;
        const mm = groupMeans(columnsM, new Array(columnsM.length).fill(1), d).labels;
        if (gm.length !== mm.length || gm.some((x, i) => x !== mm[i])) {
          out.slots.metabolite.push(err('The time points (groups) of the gene and metabolite tables differ: ' + list(gm, 5) + ' vs ' + list(mm, 5) + '. Both tables need the same groups in the same order.'));
        }
      }
    }

    // time points
    if (gOK) {
      let nT = columnsG.length; let labels = columnsG;
      if (dOK && !out.slots.design.some((x) => x.level === 'error')) {
        labels = groupMeans(columnsG, new Array(columnsG.length).fill(1), d).labels; nT = labels.length;
      }
      out.timePoints = nT;
      if (nT < 4) out.slots.gene.push(warn('Only ' + nT + ' time point' + (nT === 1 ? '' : 's') + '. Granger and CCM cannot be computed and the other metrics are unreliable; use at least 4 (more is better).'));
      out.timeLabels = labels;
    }

    // negative / non-integer counts
    if (gOK) {
      if (g.valuesNegative && dOK) out.slots.gene.push(err(plural(g.valuesNegative, 'value is', 'values are') + ' negative (e.g. ' + g.negExample + '). DESeq2 needs raw counts, which cannot be negative. Use raw counts, or remove the design table.'));
      else if (g.valuesNegative) out.slots.gene.push(warn(plural(g.valuesNegative, 'value is', 'values are') + ' negative (e.g. ' + g.negExample + '). The fold change takes log2 of the expression and will be empty for those values.'));
      if (dOK && g.nonInteger) out.slots.gene.push(warn('Values are not integers. DESeq2 needs raw counts (not TPM/FPKM); values are rounded and the fold change is not valid for normalised data.'));
    }

    // annotation
    if (ok(a)) {
      const di = a.header.findIndex((h, i) => i >= 1 && h.toLowerCase() === 'description');
      if (di < 1) {
        out.slots.annotation.push(err("The annotation table needs a column named 'Description' (the gene description); its columns are " + list(a.header, 6) + '. The first column must hold the gene names.'));
      } else if (gOK) {
        const genes = new Set(g.names);
        let hit = 0; let first = null;
        a.names.forEach((n) => { if (genes.has(n)) hit++; else if (first === null) first = n; });
        const nA = a.names.length;
        if (hit === 0) {
          out.slots.annotation.push(err("None of the " + nA + " annotation names matches a gene name (e.g. annotation '" + a.names[0] + "' vs gene '" + g.names[0] + "'). The first column of the annotation table must hold the gene names of the gene table."));
        } else if (hit / g.names.length < 0.5) {
          out.slots.annotation.push(warn(Math.round(100 * hit / g.names.length) + '% of the genes (' + hit + ' of ' + g.names.length + ') have an annotation. Check that the first column holds the gene names used in the gene table.'));
        }
        a.matchInfo = { hit, genes: g.names.length };
      }
    }

    // target metabolite
    const tname = state.target;
    if (mOK && !tname) out.todo.push('choose a target metabolite');
    if (mOK && tname) {
      const row = m.rowMap.get(tname);
      const cntT = m.names.filter((n) => n === tname).length;
      const designErr = d && !dOK ? true : out.slots.design.some((x) => x.level === 'error') || out.slots.metabolite.some((x) => x.level === 'error');
      if (!row) out.slots.target.push(err("'" + tname + "' is not a row of the metabolite table. Choose a metabolite from the list."));
      else if (cntT > 1) out.slots.target.push(err("'" + tname + "' occurs " + cntT + ' times in the metabolite table. Metabolite names must be unique; rename or delete the duplicates.'));
      else if (!designErr) {
        let labels = columnsM; let series = row;
        if (dOK) { const gm = groupMeans(columnsM, row, d); labels = gm.labels; series = gm.means; }
        out.targetSeries = { labels, series };
        const missing = labels.filter((_, i) => !Number.isFinite(series[i]));
        if (missing.length) {
          out.slots.target.push(err("'" + tname + "' has no value at " + list(missing, 4) + '. Fill in the missing value(s) or choose another metabolite.'));
        } else if (series.every((v) => v === series[0])) {
          out.slots.target.push(err("'" + tname + "' is constant over time (" + fmtNum(series[0]) + ' at every time point), so no association can be computed. Choose another metabolite.'));
        } else {
          const rise = steepestRise(labels, series);
          if (rise) {
            out.contrast = rise;
            if (dOK && gOK) {
              const cnt = new Map(); const gOf = new Map(); d.samples.forEach((s, i) => gOf.set(s, d.groups[i]));
              columnsG.forEach((c) => { const gr = gOf.get(c); cnt.set(gr, (cnt.get(gr) || 0) + 1); });
              const few = [rise.end, rise.start].filter((x) => (cnt.get(x) || 0) < 2);
              if (few.length) {
                out.slots.design.push(err('The fold change compares ' + rise.end + ' with ' + rise.start + ' (the largest rise of ' + tname + '), but ' + list(few, 2) + ' ' + (few.length === 1 ? 'has' : 'have') +
                  ' only 1 sample. DESeq2 needs at least 2 samples per group. Add replicates to the tables and the design table, or remove the design table.'));
              }
            }
          } else {
            out.slots.target.push(warn("'" + tname + "' never rises between adjacent time points. The fold change will be empty for all genes and counts as 0 in the score."));
          }
        }
      }
    }
    // groups with a single sample that are not part of the comparison
    if (dOK && gOK && out.slots.design.every((x) => x.level !== 'error')) {
      const cnt = new Map(); const gOf = new Map(); d.samples.forEach((s, i) => gOf.set(s, d.groups[i]));
      columnsG.forEach((c) => { const gr = gOf.get(c); if (gr !== undefined) cnt.set(gr, (cnt.get(gr) || 0) + 1); });
      const single = [...cnt.entries()].filter(([, n]) => n < 2).map(([k]) => k);
      if (single.length === cnt.size && cnt.size) {
        out.slots.design.push(err('Every group has a single sample, so there are no replicates. DESeq2 needs at least 2 samples per group. Remove the design table, or add replicates.'));
      } else if (single.length) {
        out.slots.design.push(warn('Group ' + list(single, 4) + ' has only 1 sample. Groups with replicates are averaged; DESeq2 needs at least 2 samples in the two groups of the fold change.'));
      }
    }

    // parameters
    const P = state.params || {};
    const nT = out.timePoints;
    const intIn = (v, lo, name, extra) => {
      if (v === '' || v === null || v === undefined || !Number.isInteger(Number(v)) || Number(v) < lo) {
        out.slots.params.push(err(name + ' must be a whole number of at least ' + lo + (extra || '') + '.')); return false;
      }
      return true;
    };
    const clusters = P.clusters === '' ? NaN : Number(P.clusters);
    if (!(Number.isInteger(clusters) && clusters >= 0)) out.slots.params.push(err('Clusters must be a whole number: 0 (no clustering) or at least 2.'));
    else if (clusters === 1) out.slots.params.push(err('Clusters must be 0 (no clustering) or at least 2.'));
    else if (gOK && clusters > g.names.length) out.slots.params.push(err('Clusters (' + clusters + ') cannot exceed the number of genes (' + g.names.length + ').'));
    const lagOK = intIn(P.lag, 1, 'Lag');
    const eOK = intIn(P.E, 1, 'E');
    const tauOK = intIn(P.tau, 1, 'Tau');
    if (nT && state.metric === 'Granger' && lagOK && nT <= 3 * Number(P.lag) + 1) {
      out.slots.params.push(warn('Granger needs more than 3 x lag + 1 = ' + (3 * Number(P.lag) + 1) + ' time points; there are ' + nT + '. All Granger values will be empty. Lower the lag or choose another metric.'));
    }
    if (nT && state.metric === 'CCM' && eOK && tauOK) {
      const emb = nT - (Number(P.E) - 1) * Number(P.tau);
      if (emb < Number(P.E) + 2) out.slots.params.push(warn('CCM needs at least E + 2 = ' + (Number(P.E) + 2) + ' points after embedding; there are ' + emb + ' (' + nT + ' time points, E = ' + P.E + ', tau = ' + P.tau + '). All CCM values will be empty. Lower E or tau, or choose another metric.'));
    }
    return out;
  }

  const api = { analyzeFile, crossChecks, steepestRise, groupMeans, splitLine, extensionOf };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.CBChecks = api;
})(typeof self !== 'undefined' ? self : this);
