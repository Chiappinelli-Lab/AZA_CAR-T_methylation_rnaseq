#!/usr/bin/env python3
"""Compare the rerun figures with the reference (current-project) figures listed in manifest.tsv.

PDF: each page rendered at 300 dpi with pdftoppm, then pixel-diffed (PIL). PNG: pixel-diffed directly.
SVG: text compared after canonicalizing what matplotlib randomizes per run (element ids drawn from a uuid salt,
and the creation date); if the text still differs, the result defers to the same figure's PDF render.
Writes verification/figure_comparison.tsv. Env data_analysis; pdftoppm from PATH or $PDFTOPPM.

usage: compare_figures.py [manifest.tsv]
"""
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MANIFEST = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "manifest.tsv"
PDFTOPPM = os.environ.get("PDFTOPPM") or shutil.which("pdftoppm")
if not PDFTOPPM:
    sys.exit("pdftoppm not found (set PDFTOPPM)")


def render_pdf(pdf, outdir):
    subprocess.run([PDFTOPPM, "-r", "300", "-png", str(pdf), str(outdir / "p")], check=True)
    return sorted(outdir.glob("p*.png"))


def pixel_diff(a, b):
    """(identical?, detail) for two images. Compares every channel with numpy: Image.getbbox() on an RGBA
    difference looks at the alpha channel only, which hid all RGB differences in an earlier version."""
    ia, ib = Image.open(a).convert("RGBA"), Image.open(b).convert("RGBA")
    if ia.size != ib.size:
        return False, f"size {ia.size} vs {ib.size}"
    differs = (np.asarray(ia) != np.asarray(ib)).any(axis=2)
    if not differs.any():
        return True, "pixel-identical"
    ys, xs = np.nonzero(differs)
    bbox = (int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1)
    return False, f"{int(differs.sum())} pixels differ ({100 * differs.mean():.4f}%), bbox {bbox}"


def compare_pdf(a, b):
    with tempfile.TemporaryDirectory() as ta, tempfile.TemporaryDirectory() as tb:
        pa, pb = render_pdf(a, Path(ta)), render_pdf(b, Path(tb))
        if len(pa) != len(pb):
            return False, f"{len(pa)} vs {len(pb)} pages"
        out = [pixel_diff(x, y) for x, y in zip(pa, pb)]
        ok = all(o for o, _ in out)
        return ok, "; ".join(d for _, d in out)


def canon_svg(path):
    s = Path(path).read_text()
    s = re.sub(r"<dc:date>.*?</dc:date>", "<dc:date/>", s, flags=re.S)
    ids = {}
    for m in re.finditer(r'id="([^"]+)"', s):
        ids.setdefault(m.group(1), f"id{len(ids)}")
    if ids:
        pat = re.compile("|".join(sorted(map(re.escape, ids), key=len, reverse=True)))
        s = pat.sub(lambda m: ids[m.group(0)], s)
    return s


rows, pdf_result = [], {}
with open(MANIFEST) as fh:
    entries = [r for r in csv.DictReader(fh, delimiter="\t") if r["type"] in ("pdf", "png", "svg")]
entries.sort(key=lambda r: {"pdf": 0, "png": 1, "svg": 2}[r["type"]])  # PDFs first: SVGs defer to them
for r in entries:
    ref, new = Path(r["current_path"]), REPO / r["new_path"]
    if not new.exists():
        rows.append((r["new_path"], r["type"], "MISSING", "not produced by the rerun"))
        continue
    if r["type"] == "pdf":
        ok, detail = compare_pdf(ref, new)
        pdf_result[str(new.with_suffix(""))] = ok
    elif r["type"] == "png":
        ok, detail = pixel_diff(ref, new)
    else:
        ok = canon_svg(ref) == canon_svg(new)
        detail = "text identical after id/date canonicalization"
        if not ok:
            twin = pdf_result.get(str(new.with_suffix("")))
            if twin is None:
                detail = "text differs; no PDF twin to check"
            else:
                ok, detail = twin, f"text differs; PDF twin {'pixel-identical' if twin else 'differs'}"
    rows.append((r["new_path"], r["type"], "PASS" if ok else "FAIL", detail))

out = HERE / "figure_comparison.tsv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["new_path", "type", "status", "detail"])
    w.writerows(sorted(rows))
n = {s: sum(1 for x in rows if x[2] == s) for s in ("PASS", "FAIL", "MISSING")}
print(f"figures: {n} -> {out}")
