/* CAT Bridge web worker: runs Python (Pyodide) off the main thread.
 * Must be started as a MODULE worker: new Worker('worker.js', { type: 'module' }).
 *
 * main -> worker
 *   {type:'init', wheelUrl}                         load Python, packages and CAT Bridge
 *   {type:'run', id, files:[{role, name, buffer}], params, wheelUrl}
 * worker -> main
 *   {type:'progress', step, detail}                  step: runtime | packages | install | deseq | read | fc | metric | cluster | annot | score
 *   {type:'ready', versions, loadMs}
 *   {type:'notice', text}
 *   {type:'result', id, ok, ...}                      see PY_BRIDGE
 *   {type:'fatal', message}
 */
'use strict';

const PYODIDE_VERSION = '314.0.8';
const PYODIDE_BASE = 'https://cdn.jsdelivr.net/pyodide/v' + PYODIDE_VERSION + '/full/';
const PYODIDE_PACKAGES = ['numpy', 'pandas', 'scipy', 'scikit-learn', 'statsmodels', 'micropip'];
// DESeq2 (only with a design table). Versions tested in the browser; pydantic-settings is an undeclared requirement of anndata.
const DESEQ_PINNED = ['pydeseq2==0.6.0', 'anndata==0.13.4', 'pydantic-settings'];
const DESEQ_LOOSE = ['pydeseq2', 'pydantic-settings'];

const post = (m, transfer) => self.postMessage(m, transfer || []);
const progress = (step, detail) => post({ type: 'progress', step, detail: detail || '' });

let pyodide = null;
let basePromise = null;
let deseqPromise = null;

const PY_BRIDGE = `
import json, time, warnings, traceback
import numpy as np
import pandas as pd

_series = None

def _clean(s, names):
    for path, shown in names.items():
        s = s.replace(path, shown)
    return s

def run_corr(params_json):
    global _series
    p = json.loads(params_json)
    names = p.get('path_names', {})
    out = {}
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        try:
            from catbridge import catbridge as cat
            from catbridge.data_io import read_upload, _check_design
            from catbridge.preprocess import repeat_aggregation_mean
            col = p['metric']
            data = cat.compute_corr(
                p['gene'], p['metabolite'], p.get('design'), p.get('annotation'),
                p['target'], p.get('cluster_count') or None,
                lag=p.get('lag', 1), E=p.get('E', 3), tau=p.get('tau', 1),
                metrics=[col], n_cpus=1, random_state=0, verbose=True)
            t1 = time.perf_counter()
            print('[catbridge-step] score', flush=True)
            score = cat.compute_score(data, col)
            fc = data.attrs.get('fc_contrast')

            full = score.copy()
            d2 = data.set_index('Name')
            extra = [c for c in d2.columns if c not in full.columns]
            full = full.join(d2[extra], on='Name')
            full.index.name = 'Rank'
            cols = ['Name']
            if 'Description' in full.columns:
                cols.append('Description')
            cols += [col, 'log2FoldChange']
            if 'Cluster' in full.columns:
                cols.append('Cluster')
            cols.append('Score')
            view = full[cols].copy()
            for c in cols:
                if c not in ('Name', 'Description', 'Cluster') and pd.api.types.is_numeric_dtype(view[c]):
                    view[c] = view[c].astype(float).round(6)
            view = view.astype(object).where(view.notna(), None)

            # series for the time-course plot (replicates averaged), re-read quietly
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                g = read_upload(p['gene']); g.columns = g.columns.map(str)
                m = read_upload(p['metabolite']); m.columns = m.columns.map(str)
                if p.get('design'):
                    dz = _check_design(read_upload(p['design']), p['design'])
                    g = repeat_aggregation_mean(g, dz)
                    m = repeat_aggregation_mean(m, dz)
            _series = np.ascontiguousarray(g.to_numpy(dtype=np.float32)).ravel()
            tvals = m.loc[p['target']].to_numpy(dtype=float)

            out = {'ok': True,
                   'columns': cols,
                   'rows': view.values.tolist(),
                   'csv': full.to_csv(),
                   'fcContrast': list(fc) if fc is not None else None,
                   'lfcShrunk': bool(data.attrs['lfc_shrunk']) if 'lfc_shrunk' in data.attrs else None,
                   'seriesNames': [str(x) for x in g.index],
                   'times': [str(x) for x in g.columns],
                   'target': [None if not np.isfinite(v) else float(v) for v in tvals],
                   'timings': {'compute_s': t1 - t0, 'total_s': time.perf_counter() - t0}}
        except BaseException as e:
            msg = e.args[0] if isinstance(e, KeyError) and e.args and isinstance(e.args[0], str) else str(e)
            out = {'ok': False, 'errorType': type(e).__name__, 'message': _clean(msg, names),
                   'traceback': _clean(traceback.format_exc(), names)}
    seen, ws = set(), []
    for w in wl:
        text = _clean(f"{w.message}", names)
        if text not in seen:
            seen.add(text)
            ws.append({'category': w.category.__name__, 'text': text})
    out['warnings'] = ws
    return json.dumps(out)
`;

// pydeseq2 0.6 in Pyodide: no processes/threads in joblib, and joblib rejects inner_max_num_threads there.
const PY_PATCH_DESEQ = `
import joblib
from pydeseq2 import default_inference as di
_init = di.DefaultInference.__init__
_d = list(_init.__defaults__); _d[-1] = 'sequential'; _init.__defaults__ = tuple(_d)
def _parallel_backend(backend, n_jobs=-1, inner_max_num_threads=None, **kw):
    return joblib.parallel_backend(backend, n_jobs=n_jobs, **kw)
di.parallel_backend = _parallel_backend
`;

// ---- output lines -> progress ------------------------------------------------------------------------------
function onLine(line) {
  let m;
  if (/^\[catbridge-step\] score/.test(line)) return progress('score', 'Scoring and ranking');
  if (!/^\[catbridge\s/.test(line)) {
    if (/Fitting|Running|Replacing|Calculating|Refitting|LFC/i.test(line)) progress('fc-detail', line.trim());
    return;
  }
  const msg = line.replace(/^\[catbridge\s+[\d.]+s\]\s*/, '');
  if (/^reading the/.test(msg)) progress('read', msg.replace(/\s\S*[\\/]\S+$/, ''));
  else if ((m = msg.match(/^fold change: log2\((\S+) \/ (\S+)\)/))) progress('fc', m[1] + ' / ' + m[2]);
  else if (/^fold change:/.test(msg)) progress('fc', '');
  else if (/^time series clustering/.test(msg)) progress('cluster', msg);
  else if (/^annotation/.test(msg)) progress('annot', '');
  else if (/^(Granger|CCM|CCA|DTW|CCF|Spearman|Pearson)$/.test(msg)) progress('metric', msg);
}

// ---- loading -----------------------------------------------------------------------------------------------
function loadBase(wheelUrl) {
  if (basePromise) return basePromise;
  basePromise = (async () => {
    const t0 = performance.now();
    progress('runtime', 'Python ' + PYODIDE_VERSION.split('.')[0] + ' (WebAssembly)');
    const { loadPyodide } = await import(PYODIDE_BASE + 'pyodide.mjs');
    pyodide = await loadPyodide({
      indexURL: PYODIDE_BASE,
      stdout: (s) => onLine(s),
      stderr: (s) => onLine(s),
    });
    progress('packages', 'numpy, pandas, scipy, scikit-learn, statsmodels');
    await pyodide.loadPackage(PYODIDE_PACKAGES, {
      messageCallback: (m) => progress('packages', m),
      errorCallback: () => {},
    });
    progress('install', 'CAT Bridge');
    const url = new URL(wheelUrl, self.location.href).href;
    const resp = await fetch(url);
    if (!resp.ok) throw new Error('Cannot load ' + url + ' (HTTP ' + resp.status + ').');
    const wheelName = url.split('/').pop().split('?')[0];
    pyodide.FS.writeFile('/tmp/' + wheelName, new Uint8Array(await resp.arrayBuffer()));
    await pyodide.runPythonAsync('import micropip\nawait micropip.install("emfs:/tmp/' + wheelName + '")');
    pyodide.runPython('import catbridge');
    pyodide.runPython(PY_BRIDGE);
    post({ type: 'ready', stage: 'base', loadMs: Math.round(performance.now() - t0) });
  })().catch((e) => { basePromise = null; throw e; });
  return basePromise;
}

function loadDeseq() {
  if (deseqPromise) return deseqPromise;
  deseqPromise = (async () => {
    progress('deseq', 'pydeseq2 (from PyPI)');
    const install = (list) => pyodide.runPythonAsync('import micropip\nawait micropip.install(' + JSON.stringify(list) + ')');
    try {
      await install(DESEQ_PINNED);
    } catch (e) {
      post({ type: 'notice', text: 'Pinned DESeq2 versions could not be installed (' + String(e.message || e).split('\n').pop() + '); using the latest versions.' });
      await install(DESEQ_LOOSE);
    }
    pyodide.runPython(PY_PATCH_DESEQ);
    const versions = JSON.parse(pyodide.runPython(`
import json, importlib.metadata as md
def v(n):
    try: return md.version(n)
    except Exception: return None
json.dumps({n: v(n) for n in ['catbridge', 'numpy', 'pandas', 'scipy', 'scikit-learn', 'statsmodels', 'pydeseq2', 'anndata']})
`));
    post({ type: 'ready', stage: 'deseq', versions });
  })().catch((e) => { deseqPromise = null; throw e; });
  return deseqPromise;
}

// ---- running -----------------------------------------------------------------------------------------------
function safeName(s) { return String(s).replace(/[^A-Za-z0-9._-]/g, '_'); }

async function run(msg) {
  const { id, files, params, wheelUrl } = msg;
  await loadBase(wheelUrl);
  if (files.some((f) => f.role === 'design')) await loadDeseq();
  const dir = '/work/' + id;
  pyodide.FS.mkdirTree(dir);
  const p = Object.assign({ cluster_count: 8, metric: 'CCM' }, params, { path_names: {} });
  for (const f of files) {
    // the extension is kept: catbridge picks the reader from it
    const path = dir + '/' + f.role + '_' + safeName(f.name);
    pyodide.FS.writeFile(path, new Uint8Array(f.buffer));
    p[f.role] = path;
    p.path_names[path] = f.name;
  }
  progress('read', 'Reading tables');
  const t0 = performance.now();
  const fn = pyodide.globals.get('run_corr');
  let out;
  try { out = JSON.parse(fn(JSON.stringify(p))); } finally { fn.destroy(); }
  for (const f of files) { try { pyodide.FS.unlink(p[f.role]); } catch (e) { /* ignore */ } }
  out.runMs = Math.round(performance.now() - t0);
  const transfer = [];
  if (out.ok) {
    const px = pyodide.globals.get('_series');
    out.series = px.toJs();           // Float32Array, row-major (genes x time points)
    px.destroy();
    transfer.push(out.series.buffer);
  }
  post(Object.assign({ type: 'result', id }, out), transfer);
}

self.onmessage = async (ev) => {
  const msg = ev.data;
  try {
    if (msg.type === 'init') {
      await loadBase(msg.wheelUrl);
      if (msg.deseq) await loadDeseq();
      post({ type: 'ready', stage: 'init' });
    } else if (msg.type === 'run') {
      await run(msg);
    }
  } catch (e) {
    post({ type: 'fatal', message: String((e && e.message) || e), stack: e && e.stack });
  }
};
