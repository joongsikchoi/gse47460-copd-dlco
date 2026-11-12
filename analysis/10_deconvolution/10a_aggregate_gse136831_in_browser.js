// Step 10a: aggregate GSE136831 single-cell counts into a pseudobulk reference, run in a browser page on
// https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136831/suppl/ (same-origin fetch; full-body GETs of large files
// failed in the pane, so files are read with 256 KB / 8 MB Range requests). Run the blocks in order in the page console.
const base = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136831/suppl/';
async function getRange(u, a, b) { let last; for (let k = 0; k < 6; k++) { try {
  const r = await fetch(u, { headers: { Range: `bytes=${a}-${b}` }, cache: 'no-store' }); if (r.status !== 206) throw new Error('status ' + r.status);
  return new Uint8Array(await r.arrayBuffer()); } catch (e) { last = e; await new Promise(z => setTimeout(z, 700)); } } throw new Error('range fail ' + a + ' ' + last); }
async function getAll(u, len, cs = 262144) { const out = new Uint8Array(len); for (let a = 0; a < len; a += cs) out.set(await getRange(u, a, Math.min(len - 1, a + cs - 1)), a); return out; }
async function gunzipText(buf) { const ds = new DecompressionStream('gzip'); const w = ds.writable.getWriter(); w.write(buf); w.close();
  const rd = ds.readable.getReader(); const dec = new TextDecoder(); const parts = [];
  while (true) { const { done, value } = await rd.read(); if (done) break; parts.push(dec.decode(value, { stream: true })); } parts.push(dec.decode()); return parts.join(''); }

// 1) metadata, genes, barcodes (sizes from the GEO listing)
const meta = await gunzipText(await getAll(base + 'GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz', 4308752));
const genesTxt = await gunzipText(await getAll(base + 'GSE136831_AllCells.GeneIDs.txt.gz', 323848));
const L = meta.trim().split('\n').slice(1).map(l => l.replace(/"/g, '').split('\t'));   // barcode order == matrix columns (checked)
const G = { AT1: ['ATI'], AT2: ['ATII'], CapA: ['VE_Capillary_A'], CapB: ['VE_Capillary_B'], VE_other: ['VE_Arterial', 'VE_Venous', 'VE_Peribronchial'],
  Lymphatic: ['Lymphatic'], Fibroblast: ['Fibroblast', 'Myofibroblast'], SMC_Pericyte: ['SMC', 'Pericyte'], Mesothelial: ['Mesothelial'],
  Airway_epi: ['Basal', 'Club', 'Goblet', 'Ionocyte', 'PNEC', 'Aberrant_Basaloid'], Ciliated: ['Ciliated'], Mac_Alveolar: ['Macrophage_Alveolar'],
  Mac_Interstitial: ['Macrophage'], Monocyte: ['cMonocyte', 'ncMonocyte'], DC: ['cDC1', 'cDC2', 'pDC', 'DC_Mature', 'DC_Langerhans'], Mast: ['Mast'],
  T_NK: ['T', 'T_Cytotoxic', 'T_Regulatory', 'NK', 'ILC_A', 'ILC_B'], B_Plasma: ['B', 'B_Plasma'] };
const names = Object.keys(G), map = {}; names.forEach((g, i) => G[g].forEach(t => map[t] = i));
const grp = new Int8Array(L.length).fill(-1), inv = new Float64Array(L.length), nPer = new Array(names.length).fill(0);
L.forEach((r, j) => { if (r[6] === 'IPF') return; const g = map[r[4]]; if (g === undefined) return; grp[j] = g; inv[j] = 1 / Number(r[1]); nPer[g]++; });

// 2) stream the 2.0 GB MatrixMarket file: rows = genes (45,947), cols = cells (312,928); S[gene, group] += count / nUMI
const U = base + 'GSE136831_RawCounts_Sparse.mtx.gz', total = 2150633813, K = names.length;
const ds = new DecompressionStream('gzip'); const w = ds.writable.getWriter(); const rd = ds.readable.getReader();
(async () => { for (let a = 0; a < total; a += 8 * 1024 * 1024) await w.write(await getRange(U, a, Math.min(total - 1, a + 8 * 1024 * 1024 - 1))); await w.close(); })();
let S = null, dims = null; const nums = [0, 0, 0]; let ni = 0, cur = 0, inNum = false, comment = false, lineStart = true;
while (true) { const { done, value } = await rd.read(); if (done) break;
  for (let p = 0; p < value.length; p++) { const c = value[p];
    if (lineStart) { comment = (c === 37); lineStart = false; }
    if (c === 10) { if (!comment) { if (inNum) nums[ni++] = cur;
        if (ni === 3) { if (dims === null) { dims = [...nums]; S = new Float64Array(dims[0] * K); }
          else { const g = grp[nums[1] - 1]; if (g >= 0) S[(nums[0] - 1) * K + g] += nums[2] * inv[nums[1] - 1]; } } }
      ni = 0; cur = 0; inNum = false; lineStart = true; comment = false; continue; }
    if (comment) continue;
    if (c >= 48 && c <= 57) { cur = cur * 10 + (c - 48); inNum = true; } else if (inNum) { nums[ni++] = cur; cur = 0; inNum = false; } } }

// 3) merge the two macrophage groups, mean CP10K per group, marker export (see README)
const sym = genesTxt.trim().split('\n').slice(1).map(l => l.replace(/"/g, '').split('\t')[1]); const n = sym.length;
const ia = names.indexOf('Mac_Alveolar'), ii = names.indexOf('Mac_Interstitial');
const nn = names.filter((x, i) => i !== ia && i !== ii).concat(['Macrophage']); const K2 = nn.length; const M2 = new Float64Array(n * K2);
for (let i = 0; i < n; i++) { let k = 0; for (let g = 0; g < K; g++) { if (g === ia || g === ii) continue; M2[i * K2 + k++] = S[i * K + g] / nPer[g] * 1e4; }
  M2[i * K2 + K2 - 1] = (S[i * K + ia] + S[i * K + ii]) / (nPer[ia] + nPer[ii]) * 1e4; }
const sel = new Set();
for (let g = 0; g < K2; g++) { const sc = []; for (let i = 0; i < n; i++) { const v = M2[i * K2 + g]; if (v < 0.5) continue; let mo = 0;
    for (let h = 0; h < K2; h++) if (h !== g && M2[i * K2 + h] > mo) mo = M2[i * K2 + h]; const r = Math.log2(v + 1) - Math.log2(mo + 1); if (r > 0.5) sc.push([r, i]); }
  sc.sort((a, b) => b[0] - a[0]); sc.slice(0, 100).forEach(x => sel.add(x[1])); }
JSON.stringify({ groups: nn, genes: [...sel].map(i => sym[i]), rows: [...sel].map(i => Array.from(M2.slice(i * K2, i * K2 + K2)).map(v => +v.toPrecision(4))) });
