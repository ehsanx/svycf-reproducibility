#!/usr/bin/env python3
"""Rebuild every display into a temporary folder and compare it, byte for byte, with results/displays/.

    python tests/check_displays.py

Exit status 0 when every file is identical (no file missing, none extra). Building Figures 1, 2 and S1 needs pdflatex
on the PATH; config/python-requirements.txt records the TeX version they were made with.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TIKZ = {"main/fig_targets.pdf", "main/fig_folds_variances.pdf", "supplement/figS_proof_map.pdf"}


def files(d):
    return {p.relative_to(d).as_posix(): p.read_bytes() for p in Path(d).rglob("*") if p.is_file()}


with tempfile.TemporaryDirectory() as tmp:
    subprocess.run([sys.executable, "-B", str(ROOT / "src" / "displays" / "make_sim_displays.py")],
                   env=dict(os.environ, DISPLAY_OUT=tmp, PYTHONDONTWRITEBYTECODE="1"), check=True)
    new, ref = files(tmp), files(ROOT / "results" / "displays")
same = sorted(k for k in ref if new.get(k) == ref[k])
differ = sorted(k for k in ref if k in new and new[k] != ref[k])
missing, extra = sorted(set(ref) - set(new)), sorted(set(new) - set(ref))
print("%d of %d display files identical" % (len(same), len(ref)))
for label, xs in (("differ", differ), ("not rebuilt", missing), ("extra", extra)):
    for x in xs:
        print("  %s: %s" % (label, x))
if TIKZ & set(differ):
    print("  (the TikZ figures are PDF output of pdflatex: another TeX version can change their bytes, not their content)")
sys.exit(0 if len(same) == len(ref) and not extra else 1)
