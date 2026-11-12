"""Re-render the manuscript figures at submission quality: TIFF 600 dpi (LZW) + EPS + PDF.
Patches Figure.savefig so every PNG the original scripts write is also written in the other formats.
Usage: python make_hires.py <12_figures_english.py> <07 script> <11 script> <analysis dir> <qc dir> <deconv dir> <out dir>
"""
import sys, runpy, warnings
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure

s12, s07, s11, A, Q, D, out = sys.argv[1:8]
Path(out).mkdir(parents=True, exist_ok=True)
_orig = Figure.savefig
made = []

def savefig(self, fname, *a, **kw):
    p = Path(str(fname))
    _orig(self, p, *a, **kw)
    if p.suffix.lower() != ".png":
        return
    kw2 = {k: v for k, v in kw.items() if k not in ("dpi", "pil_kwargs", "format")}
    _orig(self, p.with_suffix(".tif"), *a, dpi=600,
          pil_kwargs={"compression": "tiff_lzw"}, **kw2)
    _orig(self, p.with_suffix(".pdf"), *a, **kw2)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            _orig(self, p.with_suffix(".eps"), *a, **kw2)
            eps = "ok" + ("" if not w else " (%d warning)" % len(w))
        except Exception as e:
            eps = "FAILED: %s" % e
    made.append((p.name, eps))

Figure.savefig = savefig
sys.argv = [s12, s07, s11, A, Q, D, out]
runpy.run_path(s12, run_name="__main__")
print("\n=== written ===")
for n, e in made:
    print(" ", n, "-> tif/pdf/eps:", e)
