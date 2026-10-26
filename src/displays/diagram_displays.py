"""The paper's conceptual figures: TikZ diagrams compiled to vector PDFs. Reads no data.

Typesets diagrams.tex (the diagrams of diagram_macros.tex in the paper's fonts and notation) once per figure:
    fig_targets.pdf           Figure 1, populations and targets                       (main-text folder)
    fig_folds_variances.pdf   Figure 2, whole-PSU folds and the two variance scales   (main-text folder)
    figS_proof_map.pdf        Figure S1, the four terms of the expansion              (supplement folder)
The supplement typesets Figure S1 itself, from the same macro, so that its lemma numbers are cross-references; the
standalone copy prints the numbers, and when the supplement's .aux file is present beside this folder
(../appendix/appendix.aux) they are checked against it. The PDFs carry no dates, so a rerun with the same TeX
installation reproduces the same bytes.

Called by make_sim_displays.py, which passes the two output folders; by hand:
    python diagram_displays.py <main-text folder> <supplement folder>
Requires pdflatex (TeX Live or MiKTeX, with standalone, TikZ and fix-cm) on the PATH.
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGURES = (("fig_targets", "main"), ("fig_folds_variances", "main"), ("figS_proof_map", "supp"))
# The lemma numbers printed by the standalone Figure S1 (diagrams.tex), in the order S_des, T_1, T_2, S_mod.
PROOF_MAP_LEMMAS = (("lem:designclt", "1"), ("lem:crossfit", "4"), ("lem:prodrate", "3"), ("lem:modterm", "2"))


def check_lemma_numbers():
    aux = HERE.parent / "appendix" / "appendix.aux"
    if not aux.is_file():
        return
    text = aux.read_text(encoding="utf-8", errors="replace")
    for label, number in PROOF_MAP_LEMMAS:
        m = re.search(r"\\newlabel\{%s\}\{\{([^}]*)\}" % re.escape(label), text)
        if not m or m.group(1) != number:
            raise RuntimeError("Figure S1 prints Lemma %s for %s, but the supplement numbers it %s"
                               % (number, label, m.group(1) if m else "(no label)"))


def build(out_main, out_supp):
    check_lemma_numbers()
    out = {"main": Path(out_main), "supp": Path(out_supp)}
    for d in out.values():
        d.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        for name, where in FIGURES:
            try:
                r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-jobname=" + name,
                                    "-output-directory=" + tmp, "diagrams.tex"],
                                   cwd=HERE, capture_output=True, text=True, errors="replace")
            except FileNotFoundError:
                raise RuntimeError("pdflatex is not on the PATH; it typesets Figures 1, 2 and S1") from None
            if r.returncode:
                raise RuntimeError("pdflatex failed for %s:\n%s" % (name, r.stdout[-3000:]))
            log = (Path(tmp) / (name + ".log")).read_text(errors="replace")
            warnings = [ln for ln in log.splitlines() if "Warning" in ln or ln.startswith("Overfull")]
            if warnings:
                raise RuntimeError("%s: TeX warnings:\n%s" % (name, "\n".join(warnings)))
            shutil.copy(Path(tmp) / (name + ".pdf"), out[where] / (name + ".pdf"))


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2])
    print("wrote %s" % ", ".join(n + ".pdf" for n, _ in FIGURES))
