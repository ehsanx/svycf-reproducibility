#!/usr/bin/env python3
"""Build every simulation display of the manuscript from the saved grading outputs.

Nothing here is typed by hand: each number in Figure 1, Table 1, the supplement's tables and the
\\Sim... macros used in the prose is read from the saved summaries in results/simulation/ and
results/application/. Re-running the script on the
same inputs reproduces the same files.

Inputs (read only), under results/simulation/:
    grading/                                                   grades, displays
    execution-check/assemble/assemble-20261016T071738Z/        G5, R-b gate, G7
    preparation/prep-panel/prep-panel-20261012T081105Z/        frozen kappa_AE, G6 limb (a)
    unresolved_t2_by_cell_arm.csv                              unresolved-T2 counts
and the SPA aggregates under results/application/ (app_displays.py).

Outputs (overwritten), under results/displays/ (or $DISPLAY_OUT):
    main/fig1_simulation.pdf
    main/tab1_estimates.tex          (tabular body of Table 1, panel A: estimates, benchmarks)
    main/tab1_coverage.tex           (tabular body of Table 1, panel B)
    main/sim_numbers.tex             (\\newcommand macros quoted in the prose)
    supplement/tabS_*.tex            (supplement tables, complete environments, without the grading labels)
    extended/tabS_*.tex, rec_index.tex   (every table in its original form, with the original classifications and
                                      checks, and the cell index; not printed in the paper)

Run from anywhere:   python src/displays/make_sim_displays.py
Requires: python >= 3.10, pandas, numpy, matplotlib.
"""
import ast
import os
import sys
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
SIM = ROOT / "results" / "simulation"
GRADE = SIM / "grading"
DISP = GRADE / "displays"
R3D = SIM / "execution-check" / "assemble" / "assemble-20261016T071738Z"
KAPPA = SIM / "preparation" / "prep-panel" / "prep-panel-20261012T081105Z" / "kappa_ae.json"
G6 = SIM / "preparation" / "prep-panel" / "prep-panel-20261012T081105Z" / "g6.csv"
UNRESOLVED = SIM / "unresolved_t2_by_cell_arm.csv"

OUT = Path(os.environ.get("DISPLAY_OUT", ROOT / "results" / "displays"))
OUT_MAIN = OUT / "main"
OUT_SUPP = OUT / "supplement"
OUT_REC = OUT / "extended"

# ----------------------------------------------------------------------------------------------
# Formatting helpers
# ----------------------------------------------------------------------------------------------

def num(x, d=3):
    """A number for LaTeX, safe in text and math mode; '---' for a missing value."""
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "---"
    return r"\ensuremath{%.*f}" % (d, x)


def interval(lo, hi, d=3):
    if lo is None or hi is None or any(isinstance(v, float) and math.isnan(v) for v in (lo, hi)):
        return "---"
    return "[%s, %s]" % (num(lo, d), num(hi, d))


def est_int(e, lo, hi, d=3):
    if e is None or (isinstance(e, float) and math.isnan(e)):
        return "---"
    return "%s %s" % (num(e, d), interval(lo, hi, d))


def tex_escape(s):
    rep = {"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}"}
    return "".join(rep.get(c, c) for c in str(s))


def region(lo, hi, d=2):
    lo_inf = isinstance(lo, float) and math.isinf(lo) or lo is None or (isinstance(lo, float) and math.isnan(lo))
    if lo_inf:
        return r"$\le %s$" % ("%.*f" % (d, hi))
    if abs(lo + hi) < 1e-12:
        return r"$|\cdot|\le %s$" % ("%.*f" % (d, hi))
    return r"$[%s, %s]$" % ("%.*f" % (d, lo), "%.*f" % (d, hi))


STAT = {
    "mean_Z": r"mean $Z$",
    "sd_Z_minus_1": r"SD$(Z)-1$",
    "KS_Z": r"Kolmogorov distance",
    "mean_R_var_minus_1": r"mean $R_{\mathrm{var}}-1$",
    "sd_R_var": r"SD$(R_{\mathrm{var}})$",
    "transfer": r"mean transfer",
    "Q90_abs_T1_std": r"$Q_{0.9}(|T_1|/\sigma_{m,\mathrm{wor}})$",
    "Q90_abs_T2_std": r"$Q_{0.9}(|T_2|/\sigma_{m,\mathrm{wor}})$",
    "Q90_abs_paired_psi_std": r"$Q_{0.9}(|\hat\psi_{\mathrm{TMLE}}-\hat\psi_{\mathrm{AIPW}}|/\sigma_{m,\mathrm{wor}})$",
    "Q90_abs_se2_diff_std": r"$Q_{0.9}(|\widehat{\mathrm{se}}^2_{\mathrm{TMLE}}-\widehat{\mathrm{se}}^2_{\mathrm{AIPW}}|/\sigma_m^2)$",
    "Q90_bridge": r"$Q_{0.9}(|\Psi_N-\Psi|/\sigma_{m,\mathrm{wor}})$",
    "mean_Z_shift": r"mean $Z$ (shifted centre)",
    "mean_R_var_nomh": r"mean $R_{\mathrm{var}}$",
    "mean_R_var_unit": r"mean $R_{\mathrm{var}}$",
    "mean_R_var_nostrat": r"mean $R_{\mathrm{var}}$",
    "c3a_signed": r"$\widehat{\mathrm{Var}}(S_{\mathrm{des}})/\sigma_m^2-(1-\Delta_m/\sigma_m^2)$",
}

ARM = {
    "oracle": "oracle",
    "Ra-AIPW": "R-a AIPW",
    "Ra-TMLE": "R-a TMLE",
    "Rb-AIPW": "R-b AIPW",
    "Rb-TMLE": "R-b TMLE",
    "Rd-AIPW-a020": r"R-d AIPW, $\alpha=0.20$",
    "Rd-AIPW-a025": r"R-d AIPW, $\alpha=0.25$",
    "Rd-AIPW-a030": r"R-d AIPW, $\alpha=0.30$",
    "Rd-AIPW-a035": r"R-d AIPW, $\alpha=0.35$",
    "census": "census",
    "Ra-AIPW-nocf": "R-a AIPW, no cross-fitting",
    "Rc-AIPW-psucf": "R-c AIPW, PSU folds",
    "Rc-AIPW-nocf": "R-c AIPW, no cross-fitting",
    "Rc-AIPW-indcf": "R-c AIPW, unit folds",
    "Rc-AIPW-indcf-coarse": r"R-c AIPW, unit folds, coarsened $C_{\mathrm{psu}}$",
}

SETTING = {
    "P-ms-H800": r"many strata, $m=1600$",
    "P-ms-H400": r"many strata, $m=800$",
    "P-fpc-fit-J4": r"$f=0.50$, $m=100$",
    "P-strict-m200": r"$f=0.50$, $m=200$",
    "P-fs-m1600": r"fixed strata, $m=1600$",
    "P-fs-m3200": r"fixed strata, $m=3200$",
    "P-vanish-census": r"$f=0.0025$, $m=400$",
    "P-e3-ctrl-m400": r"fixed strata, $m=400$",
    "P-fpc-or-J4": r"$f=0.50$, $m=100$",
}


def short_disp(d):
    d = str(d)
    if d.startswith("RED-FLAG"):
        return r"\textbf{RED-FLAG}"
    if d.startswith("PASS"):
        return "PASS"
    if d.startswith("BOUNDARY"):
        return "BOUNDARY (descriptive)"
    if d.startswith("FAILS AS PREDICTED"):
        return "fails as predicted"
    if d.startswith("DIAGNOSTIC"):
        return "diagnostic only"
    if d.startswith("UNGRADABLE"):
        if "MC precision" in d and "premise" not in d:
            return r"\textbf{UNGRADABLE} (MC precision)"
        if "STOPPED" in d and "alpha < 1/4" not in d:
            return r"\textbf{UNGRADABLE} (stopped)"
        if "alpha < 1/4" in d:
            return r"\textbf{UNGRADABLE} (no rung)"
        if "matched AIPW" in d:
            return r"\textbf{UNGRADABLE} (premise)"
        if "quarter-rate" in d:
            return r"\textbf{UNGRADABLE} (premise)"
        if "Z_oracle" in d:
            return r"\textbf{UNGRADABLE} (premise)"
        return r"\textbf{UNGRADABLE}"
    return tex_escape(d)


# ----------------------------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------------------------
grades = pd.read_csv(GRADE / "grades.csv")
claims = pd.read_csv(GRADE / "claims.csv")
complete = pd.read_csv(GRADE / "completeness.csv")
reqchk = pd.read_csv(GRADE / "required_checks.csv")
popsum = pd.read_csv(GRADE / "population_indexed_summary.csv")
provchk = pd.read_csv(GRADE / "provenance_checks.csv")
ms = pd.read_csv(DISP / "e1_many_strata_by_arm.csv")
fs = pd.read_csv(DISP / "e1_fixed_strata_by_arm.csv")
e1o = pd.read_csv(DISP / "e1_oracle_by_population.csv")
e2f = pd.read_csv(DISP / "e2_fitted_by_arm.csv")
fpc = pd.read_csv(DISP / "e2_fpc_oracle_by_population.csv")
strict = pd.read_csv(DISP / "e2_strictness_witness.csv")
uneq = pd.read_csv(DISP / "e2_unequal_pair.csv")
vanish = pd.read_csv(DISP / "e2_vanishing_trajectory.csv")
ctrl_rep = pd.read_csv(DISP / "e3_controls_reported.csv")
rc = pd.read_csv(DISP / "e3_rc_foil_by_arm.csv")
vs = pd.read_csv(DISP / "e3_v_sensitivity_by_arm.csv")
e4 = pd.read_csv(DISP / "e4_by_rung.csv")
g5 = pd.read_csv(R3D / "g5.csv")
rbg = pd.read_csv(R3D / "rb_gates.csv")
g7 = pd.read_csv(R3D / "g7.csv")
g6 = pd.read_csv(G6)
kappa = json.loads(KAPPA.read_text())
unresolved = {(r.cell, r.arm): int(r.unresolved) for r in pd.read_csv(UNRESOLVED).itertuples()}


def grow(claim, cell, arm, stat):
    """One graded row of grades.csv."""
    r = grades[(grades.claim == claim) & (grades.cell == cell) & (grades.arm == arm) & (grades.statistic == stat)]
    assert len(r) == 1, (claim, cell, arm, stat, len(r))
    return r.iloc[0]


def disp_row(df, **kw):
    r = df
    for k, v in kw.items():
        r = r[r[k] == v]
    assert len(r) == 1, (kw, len(r))
    return r.iloc[0]


# ----------------------------------------------------------------------------------------------
# Macros quoted in the prose
# ----------------------------------------------------------------------------------------------
M = {}


def mac(name, value):
    assert re.fullmatch(r"(Sim|App)[A-Za-z]+", name), name
    assert name not in M, name
    M[name] = value


def mac_grade(prefix, claim, cell, arm, stat, d=3):
    r = grow(claim, cell, arm, stat)
    mac(prefix, num(r.estimate, d))
    mac(prefix + "Int", interval(r.ci_lo, r.ci_hi, d))
    return r


# C1: many strata, R-a AIPW, m = 1600 (PASS)
mac_grade("SimCoMsMeanZ", "C1", "P-ms-H800", "Ra-AIPW", "mean_Z")
mac_grade("SimCoMsSD", "C1", "P-ms-H800", "Ra-AIPW", "sd_Z_minus_1")
mac_grade("SimCoMsKS", "C1", "P-ms-H800", "Ra-AIPW", "KS_Z")
# C1: f = 0.50, m = 100, R-a AIPW (RED-FLAG)
mac_grade("SimCoFpcMeanZ", "C1", "P-fpc-fit-J4", "Ra-AIPW", "mean_Z")
r_sd = mac_grade("SimCoFpcSD", "C1", "P-fpc-fit-J4", "Ra-AIPW", "sd_Z_minus_1")
mac("SimCoFpcSDZ", num(1 + r_sd.estimate, 3))
mac("SimCoFpcSDZInt", interval(1 + r_sd.ci_lo, 1 + r_sd.ci_hi, 3))
mac("SimCoFpcVarExcessPct", "%.0f" % (100 * ((1 + r_sd.estimate) ** 2 - 1)))
mac("SimCoFpcSDExcessPct", "%.0f" % (100 * r_sd.estimate))
mac_grade("SimCoFpcKS", "C1", "P-fpc-fit-J4", "Ra-AIPW", "KS_Z")
# C1: f = 0.50, m = 200, R-a AIPW (UNGRADABLE, MC precision)
mac_grade("SimCoStrMeanZ", "C1", "P-strict-m200", "Ra-AIPW", "mean_Z")
mac_grade("SimCoStrSD", "C1", "P-strict-m200", "Ra-AIPW", "sd_Z_minus_1")
mac_grade("SimCoStrKS", "C1", "P-strict-m200", "Ra-AIPW", "KS_Z")
# R-b rate gate
rb = rbg.iloc[0]
mac("SimRbSlope", num(rb.slope, 3))
mac("SimRbUpper", num(rb.upper_95, 3))
mac("SimRbMeanZEight", num(disp_row(ms, cell="P-ms-H400", arm="Rb-AIPW").mean_Z, 2))
mac("SimRbMeanZFifty", num(disp_row(ms, cell="P-ms-H25", arm="Rb-AIPW").mean_Z, 2))
# E2 fitted cell, components of coverage
fa = disp_row(e2f, cell="P-fpc-fit-J4", arm="Ra-AIPW")
fo = disp_row(e2f, cell="P-fpc-fit-J4", arm="oracle")
mac("SimEtRvar", num(fa.mean_R_var, 3))
mac("SimEtRvarPct", "%.1f" % (100 * (fa.mean_R_var - 1)))
mac("SimEtRvarSqrt", num(math.sqrt(fa.mean_R_var), 3))
mac("SimEtRvarSqrtPct", "%.0f" % (100 * (math.sqrt(fa.mean_R_var) - 1)))
mac("SimEtCovNinetyFive", num(fa.cov_se_95, 3))
mac("SimEtCovNinety", num(fa.cov_se_90, 3))
mac("SimEtCovEighty", num(fa.cov_se_80, 3))
mac("SimEtOracleCovNinetyFive", num(fo.cov_se_95, 3))
mac("SimEtOracleSDZ", num(fo.sd_Z, 3))
mac("SimEtOracleRvar", num(fo.mean_R_var, 3))
j4 = fpc[fpc.J == 4]
mac("SimEtSigRatioMin", num(j4.sigma_ratio.min(), 2))
mac("SimEtSigRatioMax", num(j4.sigma_ratio.max(), 2))
mac("SimEtSigRatioMed", num(j4.sigma_ratio.median(), 2))
sa = disp_row(e2f, cell="P-strict-m200", arm="Ra-AIPW")
mac("SimStrSDZ", num(sa.sd_Z, 3))
mac("SimStrCovNinetyFive", num(sa.cov_se_95, 3))
mac("SimStrRvar", num(sa.mean_R_var, 3))
# C2
for nm, cell, arm in (("Ra", "P-ms-H800", "Ra-AIPW"), ("Rb", "P-ms-H400", "Rb-AIPW")):
    mac_grade("SimCt%sMean" % nm, "C2", cell, arm, "mean_R_var_minus_1")
    mac_grade("SimCt%sSD" % nm, "C2", cell, arm, "sd_R_var")
    mac_grade("SimCt%sTr" % nm, "C2", cell, arm, "transfer")
rbm = grow("C2", "P-ms-H400", "Rb-AIPW", "mean_R_var_minus_1").estimate
mac("SimCtRbPct", "%.1f" % (100 * rbm))
mac("SimCtRbSEPct", "%.1f" % (100 * (math.sqrt(1 + rbm) - 1)))
ram = grow("C2", "P-ms-H800", "Ra-AIPW", "mean_R_var_minus_1").estimate
mac("SimCtRaPct", "%.1f" % (100 * ram))
# C3a
c3a = grades[grades.claim == "C3a"]
mac("SimCthaPass", str(int((c3a.disposition == "PASS").sum())))
mac("SimCthaN", str(len(c3a)))
mac("SimCthaMin", num(c3a.estimate.min(), 3))
mac("SimCthaMax", num(c3a.estimate.max(), 3))
mac("SimCthaLoMin", num(c3a.ci_lo.min(), 3))
mac("SimCthaHiMax", num(c3a.ci_hi.max(), 3))
# C3b counts by fraction (identical across the three levels; asserted)
JF = {100: "0.02", 20: "0.10", 8: "0.25", 4: "0.50"}
WORD = {100: "Hundred", 20: "Twenty", 8: "Eight", 4: "Four"}
c3b_counts = {}
for J in (100, 20, 8, 4):
    rows = popsum[(popsum.claim == "C3b") & (popsum.cell == "P-fpc-or-J%d" % J)]
    trip = set(zip(rows.PASS, rows.UNGRADABLE, rows.RED_FLAG, rows.populations))
    assert len(rows) == 3 and len(trip) == 1, (J, rows)
    p, u, rf, n = trip.pop()
    c3b_counts[J] = (p, u, rf, n)
    mac("SimCtbPass" + WORD[J], str(p))
    mac("SimCtbUng" + WORD[J], str(u))
mac("SimCtbRedFlag", str(int(popsum[popsum.claim == "C3b"].RED_FLAG.sum())))
mac("SimCtbJFourCovDetMed", num(j4.cov_det_sigma_m_95.median(), 3))
mac("SimCtbJFourCovDetMin", num(j4.cov_det_sigma_m_95.min(), 3))
mac("SimCtbJFourCovDetMax", num(j4.cov_det_sigma_m_95.max(), 3))
mac("SimCtbJFourOverlayMed", num(j4.overlay_95.median(), 3))
mac("SimCtbJFourCovWorMed", num(j4.cov_wor_95.median(), 3))
# C4
fsr = fs[fs.arm == "Ra-AIPW"].sort_values("m")
mac("SimCfRaFirst", num(fsr.q90_abs_T1_std.iloc[0], 3))
mac("SimCfRaEightHundred", num(disp_row(fs, cell="P-fs-m800", arm="Ra-AIPW").q90_abs_T1_std, 3))
mac_grade("SimCfRaQ", "C4a", "P-fs-m1600", "Ra-AIPW", "Q90_abs_T1_std")
ratios = fsr.q90_abs_T1_std.values[:-1] / fsr.q90_abs_T1_std.values[1:]
mac("SimCfRaDropMin", "%.1f" % ratios.min())
mac("SimCfRaDropMax", "%.1f" % ratios.max())
mac("SimCfRaEnvMin", num(fsr.q90_sqrt_m_T1_over_rhat.min(), 2))
mac("SimCfRaEnvMax", num(fsr.q90_sqrt_m_T1_over_rhat.max(), 2))
mac("SimCfRaTtwoTerm", num(disp_row(fs, cell="P-fs-m1600", arm="Ra-AIPW").q90_abs_T2_std, 3))
mac_grade("SimCfDqT", "C4a", "P-fs-m3200", "Rd-AIPW-a025", "Q90_abs_T1_std")
mac_grade("SimCfDtT", "C4a", "P-fs-m3200", "Rd-AIPW-a030", "Q90_abs_T1_std")
mac_grade("SimCfDtTwo", "C4b", "P-fs-m3200", "Rd-AIPW-a030", "Q90_abs_T2_std")
mac_grade("SimCfDqTwo", "C4b", "P-fs-m3200", "Rd-AIPW-a025", "Q90_abs_T2_std")
d30 = fs[fs.arm == "Rd-AIPW-a030"].sort_values("m")
d25 = fs[fs.arm == "Rd-AIPW-a025"].sort_values("m")
mac("SimCfDtTwoFirst", num(d30.q90_abs_T2_std.iloc[0], 3))
mac("SimCfDqTwoEight", num(disp_row(fs, cell="P-fs-m800", arm="Rd-AIPW-a025").q90_abs_T2_std, 3))
mac("SimCfDqTwoSixteen", num(disp_row(fs, cell="P-fs-m1600", arm="Rd-AIPW-a025").q90_abs_T2_std, 3))
mac("SimKappaSigmaMin", num(min(kappa["sigma_wor"]), 4))
mac("SimKappaSigmaMax", num(max(kappa["sigma_wor"]), 4))
mac("SimTolStd", r"\ensuremath{%.0f\times10^{-5}}" % (1e5 * 5e-7 / min(kappa["sigma_wor"])))
mac("SimCfDtTFirst", num(d30.q90_abs_T1_std.iloc[0], 3))
mac("SimCfDqTwoFirst", num(d25.q90_abs_T2_std.iloc[0], 3))
mac("SimCfDtTwoDropPct", "%.0f" % (100 * (1 - d30.q90_abs_T2_std.iloc[-1] / d30.q90_abs_T2_std.iloc[0])))
mac("SimCfDtEnvTwoFirst", num(d30.q90_T2_over_rQ_rg.iloc[0], 2))
mac("SimCfDtEnvTwoLast", num(d30.q90_T2_over_rQ_rg.iloc[-1], 2))
mac("SimCfDtTwoTimes", "%.1f" % (grow("C4b", "P-fs-m3200", "Rd-AIPW-a030", "Q90_abs_T2_std").estimate / 0.10))
mac("SimCfRaTimes", "%.1f" % (grow("C4a", "P-fs-m1600", "Ra-AIPW", "Q90_abs_T1_std").estimate / 0.10))
fa16 = disp_row(fs, cell="P-fs-m1600", arm="Ra-AIPW")
fo16 = disp_row(fs, cell="P-fs-m1600", arm="oracle")
mac("SimFsRaCov", num(fa16.cover95_studentized, 3))
mac("SimFsOrCov", num(fo16.cover95_studentized, 3))
mac("SimFsRaSDZ", num(fa16.sd_Z, 3))
fo32 = disp_row(fs, cell="P-fs-m3200", arm="oracle")
mac("SimFsOrCovThree", num(fo32.cover95_studentized, 3))
mac("SimFsDqCov", num(disp_row(fs, cell="P-fs-m3200", arm="Rd-AIPW-a025").cover95_studentized, 3))
mac("SimFsDtCov", num(disp_row(fs, cell="P-fs-m3200", arm="Rd-AIPW-a030").cover95_studentized, 3))
mac("SimFsDqSDZ", num(disp_row(fs, cell="P-fs-m3200", arm="Rd-AIPW-a025").sd_Z, 3))
mac("SimFsDtSDZ", num(disp_row(fs, cell="P-fs-m3200", arm="Rd-AIPW-a030").sd_Z, 3))
_gaps = [100 * (fo32.cover95_studentized - disp_row(fs, cell="P-fs-m3200", arm=a).cover95_studentized)
         for a in ("Rd-AIPW-a025", "Rd-AIPW-a030")]
mac("SimFsCovGapMin", "%.1f" % min(_gaps))
mac("SimFsCovGapMax", "%.1f" % max(_gaps))
# The same coverages as percentages, for the prose
for nm, val in (("SimEtCovNinetyFivePct", fa.cov_se_95), ("SimEtOracleCovNinetyFivePct", fo.cov_se_95),
                ("SimFsRaCovPct", fa16.cover95_studentized), ("SimFsOrCovPct", fo16.cover95_studentized),
                ("SimFsOrCovThreePct", fo32.cover95_studentized),
                ("SimFsDqCovPct", disp_row(fs, cell="P-fs-m3200", arm="Rd-AIPW-a025").cover95_studentized),
                ("SimFsDtCovPct", disp_row(fs, cell="P-fs-m3200", arm="Rd-AIPW-a030").cover95_studentized),
                ("SimRbCovEightPct", disp_row(ms, cell="P-ms-H400", arm="Rb-AIPW").cover95_studentized)):
    mac(nm, "%.1f" % (100 * float("%.3f" % val)))  # from the 3-decimal value shown in the tables, so they agree
# G5 slopes
for a, nm in ((0.20, "Twenty"), (0.25, "TwentyFive"), (0.30, "Thirty"), (0.35, "ThirtyFive")):
    r = g5[np.isclose(g5.alpha, a)].iloc[0]
    mac("SimGfSlope" + nm, num(r.slope, 3))
    mac("SimGfInt" + nm, interval(r.lo, r.hi, 3))
# C6
mac_grade("SimCsPsi", "C6", "P-ms-H400", "Ra-TMLE", "Q90_abs_paired_psi_std", 4)
mac_grade("SimCsSe", "C6", "P-ms-H400", "Ra-TMLE", "Q90_abs_se2_diff_std", 4)
# C7
for stat, nm in (("mean_Z", "Mean"), ("sd_Z_minus_1", "SD"), ("KS_Z", "KS")):
    r = popsum[(popsum.claim == "C7") & (popsum.statistic == stat)].iloc[0]
    mac("SimCsevPass" + nm, str(r.PASS))
    mac("SimCsevUng" + nm, str(r.UNGRADABLE))
# bridge
mac_grade("SimBridge", "sec. Psi", "P-vanish-census", "census", "Q90_bridge")
vt = vanish.sort_values("f", ascending=False)
mac("SimBridgeFirst", num(vt.q90_bridge.iloc[0], 3))
mac("SimBridgeLastTraj", num(vanish[(vanish.H == 200) & (vanish.J == 80)].q90_bridge.iloc[0], 3))
mac("SimBridgeRatioMin", num(vanish.ratio_mean_Delta_over_var_Smod.min(), 2))
mac("SimBridgeRatioMax", num(vanish.ratio_mean_Delta_over_var_Smod.max(), 2))
# strictness witness
mac("SimStrictMin", num(strict.min_Delta_over_sigma2_m.min(), 2))
mac("SimStrictMax", num(strict.max_Delta_over_sigma2_m.max(), 2))
mac("SimStrictMedMin", num(strict.median_Delta_over_sigma2_m.min(), 2))
mac("SimStrictMedMax", num(strict.median_Delta_over_sigma2_m.max(), 2))
# power controls
for nm, arm_key in (("Shift", "C1 shift"), ("Unit", "C2 wrong clustering: unit-level"),
                    ("Align", "C4 aligned error"), ("Nomh", "C2 wrong clustering: no m_h"),
                    ("Nostrat", "C2 wrong clustering: no stratification")):
    r = grades[(grades.claim == "Power control") & grades.arm.str.startswith(arm_key)]
    assert len(r) == 1, arm_key
    r = r.iloc[0]
    mac("SimPc" + nm, num(r.estimate, 3))
    mac("SimPc" + nm + "Int", interval(r.ci_lo, r.ci_hi, 3))
mac("SimPcGformula", num(ctrl_rep.iloc[0]["mean"], 3))
mac("SimKappaAE", num(kappa["kappa"], 4))
# execution record
mac("SimCells", str(len(complete)))
mac("SimCellsComplete", str(int((complete.status == "COMPLETE").sum())))
mac("SimReqChecks", "{:,}".format(int(reqchk.required_checks.sum())).replace(",", "{,}"))
mac("SimReqFailed", str(int(reqchk.failed.sum())))
mac("SimUnresolvedTotal", str(sum(unresolved.values())))
mac("SimUnresolvedTMLE", str(sum(v for (c, a), v in unresolved.items() if a.endswith("TMLE"))))
mac("SimUnresolvedOther", str(sum(v for (c, a), v in unresolved.items() if not a.endswith("TMLE"))))
mac("SimFailedFits", str(int(grades.failed_replicates.fillna(0).sum())))


# ----------------------------------------------------------------------------------------------
# Table 1, panel A: dispositions in the predeclared order (five columns: cell | statistic | estimate [MC
# interval] | region | disposition), compact enough for a 390pt text block
# ----------------------------------------------------------------------------------------------
rows = []
SHORT = {
    "mean_Z": r"mean $Z$", "sd_Z_minus_1": r"SD$(Z)-1$", "KS_Z": "KS",
    "mean_R_var_minus_1": r"mean $R_{\mathrm{var}}-1$", "sd_R_var": r"SD$(R_{\mathrm{var}})$", "transfer": "transfer",
    "Q90_abs_T1_std": r"$Q_{0.9}|T_1|$", "Q90_abs_T2_std": r"$Q_{0.9}|T_2|$",
    "Q90_abs_paired_psi_std": r"$Q_{0.9}|\Delta\hat\psi|$", "Q90_abs_se2_diff_std": r"$Q_{0.9}|\Delta\widehat{\mathrm{se}}^2|$",
    "Q90_bridge": r"$Q_{0.9}|\Psi_N-\Psi|$", "mean_Z_shift": r"mean $Z$",
    "mean_R_var_nomh": r"mean $R_{\mathrm{var}}$", "mean_R_var_unit": r"mean $R_{\mathrm{var}}$",
    "mean_R_var_nostrat": r"mean $R_{\mathrm{var}}$",
}
CELLTXT = {
    ("Ra-AIPW", "P-ms-H800"): r"R-a AIPW; many strata, $m=1600$",
    ("Ra-AIPW", "P-fpc-fit-J4"): r"R-a AIPW; $f=0.50$, $m=100$",
    ("Ra-AIPW", "P-strict-m200"): r"R-a AIPW; $f=0.50$, $m=200$",
    ("Rb-AIPW", "P-ms-H400"): r"R-b AIPW; many strata, $m=800$",
    ("Rb-AIPW", "P-fpc-fit-J4"): r"R-b AIPW; $f=0.50$, $m=100$",
    ("Rb-AIPW", "P-strict-m200"): r"R-b AIPW; $f=0.50$, $m=200$",
    ("Ra-AIPW", "P-fs-m1600"): r"R-a AIPW; fixed strata, $m=1600$",
    ("Rd-AIPW-a020", "P-fs-m3200"): r"R-d, $\alpha=0.20$; $m=3200$",
    ("Rd-AIPW-a025", "P-fs-m3200"): r"R-d, $\alpha=0.25$; $m=3200$",
    ("Rd-AIPW-a030", "P-fs-m3200"): r"R-d, $\alpha=0.30$; $m=3200$",
    ("Rd-AIPW-a035", "P-fs-m3200"): r"R-d, $\alpha=0.35$; $m=3200$",
    ("Ra-TMLE", "P-ms-H400"): r"R-a TMLE; many strata, $m=800$",
    ("Rb-TMLE", "P-ms-H400"): r"R-b TMLE; many strata, $m=800$",
    ("census", "P-vanish-census"): r"census; $f=0.0025$, $H=200$",
}


def group(title):
    if rows:
        rows.append(r"\midrule")
    rows.append(r"\multicolumn{5}{@{}l}{\textit{%s}}\\" % title)


def line(cell, stat, estimate, reg, disp):
    rows.append("%s & %s & %s & %s & %s\\\\" % (cell, stat, estimate, reg, disp))


def graded_cell(claim, cell, arm, stats):
    first = True
    for s in stats:
        r = grow(claim, cell, arm, s)
        d = 4 if ("paired_psi" in s or "se2_diff" in s) else 3
        line(CELLTXT[(arm, cell)] if first else "", SHORT[s], est_int(r.estimate, r.ci_lo, r.ci_hi, d),
             region(r.region_lo, r.region_hi), short_disp(r.disposition))
        first = False


group(r"C1: point limit about $\Psi_N$, cross-fitted estimator")
graded_cell("C1", "P-ms-H800", "Ra-AIPW", ["mean_Z", "sd_Z_minus_1", "KS_Z"])
graded_cell("C1", "P-fpc-fit-J4", "Ra-AIPW", ["mean_Z", "sd_Z_minus_1", "KS_Z"])
graded_cell("C1", "P-strict-m200", "Ra-AIPW", ["mean_Z", "sd_Z_minus_1", "KS_Z"])
for cell in ("P-ms-H400", "P-fpc-fit-J4", "P-strict-m200"):
    r = grow("C1", cell, "Rb-AIPW", "mean_Z")
    line(CELLTXT[("Rb-AIPW", cell)], "all three", "---", "", short_disp(r.disposition))
group(r"C2: variance estimation")
graded_cell("C2", "P-ms-H800", "Ra-AIPW", ["mean_R_var_minus_1", "sd_R_var", "transfer"])
graded_cell("C2", "P-ms-H400", "Rb-AIPW", ["mean_R_var_minus_1", "sd_R_var", "transfer"])
group(r"C3a: without-replacement variance gap (oracle; counts of populations)")
c3a = grades[grades.claim == "C3a"]
line(r"oracle; $f=0.50$, $m=100$", "gap error",
     "%d/%d (of %d)" % ((c3a.disposition == "PASS").sum(), (c3a.disposition != "PASS").sum(), len(c3a)),
     region(-0.05, 0.05), "PASS (each population)")
group(r"C3b: conservativeness of $\hat\psi_{\mathrm{oracle}}\pm z\sigma_m$ (counts at levels 0.80; 0.90; 0.95)")
for J in (100, 20, 8, 4):
    per = []
    for lev in ("0.80", "0.90", "0.95"):
        rr = popsum[(popsum.claim == "C3b") & (popsum.cell == "P-fpc-or-J%d" % J) &
                    (popsum.statistic == "coverage_minus_%s" % lev)].iloc[0]
        assert rr.RED_FLAG == 0 and rr.BLOCKED == 0 and rr.NOT_GRADED == 0
        per.append("%d/%d" % (rr.PASS, rr.UNGRADABLE))
    n = int(rr.populations)
    line(r"oracle; $f=%s$, $m=100$" % JF[J], r"coverage$-$level", "; ".join(per) + " (of %d)" % n,
         r"$\ge -0.02$", r"mixed (premise); no RED-FLAG" if rr.UNGRADABLE else "PASS")
group(r"C4a: empirical-process remainder $T_1$, terminal rung")
graded_cell("C4a", "P-fs-m1600", "Ra-AIPW", ["Q90_abs_T1_std"])
graded_cell("C4a", "P-fs-m3200", "Rd-AIPW-a025", ["Q90_abs_T1_std"])
graded_cell("C4a", "P-fs-m3200", "Rd-AIPW-a030", ["Q90_abs_T1_std"])
for a in ("Rd-AIPW-a020", "Rd-AIPW-a035"):
    r = grow("C4a", "P-fs-m3200", a, "Q90_abs_T1_std")
    line(CELLTXT[(a, "P-fs-m3200")], SHORT["Q90_abs_T1_std"], "not run", region(r.region_lo, r.region_hi),
         short_disp(r.disposition))
group(r"C4b: product remainder $T_2$, terminal rung")
graded_cell("C4b", "P-fs-m3200", "Rd-AIPW-a025", ["Q90_abs_T2_std"])
graded_cell("C4b", "P-fs-m3200", "Rd-AIPW-a030", ["Q90_abs_T2_std"])
for a in ("Rd-AIPW-a020", "Rd-AIPW-a035"):
    r = grow("C4b", "P-fs-m3200", a, "Q90_abs_T2_std")
    line(CELLTXT[(a, "P-fs-m3200")], SHORT["Q90_abs_T2_std"], "not run", region(r.region_lo, r.region_hi),
         short_disp(r.disposition))
group(r"C5: premise separation")
r = grades[grades.claim == "C5"].iloc[0]
line(r"R-d, $\alpha<1/4$; $m=3200$", "derived", "---", "", short_disp(r.disposition))
group(r"C6: TMLE corollary, paired with the matched AIPW fit")
graded_cell("C6", "P-ms-H400", "Ra-TMLE", ["Q90_abs_paired_psi_std", "Q90_abs_se2_diff_std"])
r = grow("C6", "P-ms-H400", "Rb-TMLE", "Q90_abs_paired_psi_std")
line(CELLTXT[("Rb-TMLE", "P-ms-H400")], "both", "---", "", short_disp(r.disposition))
group(r"C7: design-layer implementation check (oracle; counts of array rungs)")
first = True
for stat in ("mean_Z", "sd_Z_minus_1", "KS_Z"):
    rr = popsum[(popsum.claim == "C7") & (popsum.statistic == stat)].iloc[0]
    reg = grades[(grades.claim == "C7") & (grades.statistic == stat)].iloc[0]
    line(r"deterministic array, $H=25,\dots,400$" if first else "", SHORT[stat],
         "%d/%d (of %d)" % (rr.PASS, rr.UNGRADABLE, rr.populations), region(reg.region_lo, reg.region_hi),
         "PASS" if rr.UNGRADABLE == 0 else r"mixed (MC precision)")
    first = False
group(r"Secondary target: bridge from $\Psi_N$ to $\Psi$")
graded_cell("sec. Psi", "P-vanish-census", "census", ["Q90_bridge"])
group(r"Power controls, each predicted to fail (fixed strata, $m=400$)")
pcs = grades[grades.claim == "Power control"]
PCL = {"C1 shift": r"C1: centre shifted by $0.30\,\sigma_{m,\mathrm{wor}}$",
       "C2 wrong clustering: unit-level": "C2: units as clusters",
       "C2 wrong clustering: no m_h": r"C2: no $m_h/(m_h-1)$ factor",
       "C2 wrong clustering: no stratification": "C2: strata ignored",
       "C4 aligned error": r"C4: aligned nuisance error"}
order = ["C1 shift", "C2 wrong clustering: unit-level", "C4 aligned error", "C2 wrong clustering: no m_h",
         "C2 wrong clustering: no stratification"]
for key in order:
    r = pcs[pcs.arm.str.startswith(key)].iloc[0]
    line(PCL[key], SHORT[r.statistic], est_int(r.estimate, r.ci_lo, r.ci_hi), region(r.region_lo, r.region_hi),
         short_disp(r.disposition))
OUT_MAIN.mkdir(parents=True, exist_ok=True)
OUT_SUPP.mkdir(parents=True, exist_ok=True)
# The classifications, in the order of the margin plan and with its labels, are in the extended table
# tabS_dispositions.tex (results/displays/extended); the main-text Table 1, panel A (below) reports the same estimates and
# intervals beside their benchmarks without a disposition column.
_disp_tab = rows

# ----------------------------------------------------------------------------------------------
# Table 1, panel A (main text): estimates with 95% Monte Carlo intervals beside the study's benchmarks
# ----------------------------------------------------------------------------------------------
est_rows = []


def egroup(title):
    if est_rows:
        est_rows.append(r"\midrule")
    est_rows.append(r"\multicolumn{4}{@{}l}{\textit{%s}}\\" % title)


def eline(setting, stat, estimate, bench):
    est_rows.append("%s & %s & %s & %s\\\\" % (setting, stat, estimate, bench))


def est_cell(claim, cell, arm, stats, setting):
    first = True
    for s in stats:
        r = grow(claim, cell, arm, s)
        d = 4 if ("paired_psi" in s or "se2_diff" in s) else 3
        eline(setting if first else "", SHORT[s], est_int(r.estimate, r.ci_lo, r.ci_hi, d),
              region(r.region_lo, r.region_hi))
        first = False


egroup(r"Normal approximation of $Z$; correctly specified learner (R-a), AIPW")
est_cell("C1", "P-ms-H800", "Ra-AIPW", ["mean_Z", "sd_Z_minus_1", "KS_Z"], r"many strata, $m=1600$")
est_cell("C1", "P-fpc-fit-J4", "Ra-AIPW", ["mean_Z", "sd_Z_minus_1", "KS_Z"], r"$f=0.50$, $m=100$")
est_cell("C1", "P-strict-m200", "Ra-AIPW", ["mean_Z", "sd_Z_minus_1", "KS_Z"], r"$f=0.50$, $m=200$")
egroup(r"Variance estimation, $R_{\mathrm{var}}=\widehat{\mathrm{se}}^2/\sigma_m^2$; AIPW, many strata")
est_cell("C2", "P-ms-H800", "Ra-AIPW", ["mean_R_var_minus_1", "sd_R_var", "transfer"], r"R-a, $m=1600$")
est_cell("C2", "P-ms-H400", "Rb-AIPW", ["mean_R_var_minus_1", "sd_R_var", "transfer"], r"R-b (forest), $m=800$")
egroup(r"Oracle estimator, 12 fixed populations ($H=50$, $m=100$)")
eline(r"$f=0.50$", "gap error, range", "%s to %s [%s, %s]" % (
    num(c3a.estimate.min()), num(c3a.estimate.max()), num(c3a.ci_lo.min()), num(c3a.ci_hi.max())),
    region(-0.05, 0.05))
_low = {}
for J in (100, 20, 8, 4):
    d = fpc[fpc.J == J]
    _low[J] = min((d["cov_det_sigma_m_%s" % z] - lev).min() for lev, z in ((0.80, "80"), (0.90, "90"), (0.95, "95")))
    eline(r"$f=%s$" % JF[J], r"lowest coverage$-$level", num(_low[J]), r"$\ge -0.02$")
mac("SimCtbWorstShortfall", num(-min(_low.values())))
egroup(r"Remainders at the largest sample size, fixed strata ($H=5$); $Q_{0.9}|T_j|$ standardized by $\sigma_{m,\mathrm{wor}}$")
est_cell("C4a", "P-fs-m1600", "Ra-AIPW", ["Q90_abs_T1_std"], r"R-a, $m=1600$")
# Rows without a Monte Carlo interval keep the reference value beside them; the caption gives the reason for each
# kind (R-a's T2 and the larger bridge fractions were summarized without an interval; T2 at alpha = 0.25 is at the
# boundary rate of the product-rate condition).
eline("", SHORT["Q90_abs_T2_std"], "%s (no interval)" % num(disp_row(fs, cell="P-fs-m1600", arm="Ra-AIPW").q90_abs_T2_std),
      r"$\le 0.10$")
est_cell("C4a", "P-fs-m3200", "Rd-AIPW-a025", ["Q90_abs_T1_std"], r"R-d, $\alpha=0.25$, $m=3200$")
_r = grow("C4b", "P-fs-m3200", "Rd-AIPW-a025", "Q90_abs_T2_std")
eline("", SHORT["Q90_abs_T2_std"], est_int(_r.estimate, _r.ci_lo, _r.ci_hi), r"$\le 0.10$; boundary rate")
est_cell("C4a", "P-fs-m3200", "Rd-AIPW-a030", ["Q90_abs_T1_std"], r"R-d, $\alpha=0.30$, $m=3200$")
_r = grow("C4b", "P-fs-m3200", "Rd-AIPW-a030", "Q90_abs_T2_std")
eline("", SHORT["Q90_abs_T2_std"], est_int(_r.estimate, _r.ci_lo, _r.ci_hi), region(_r.region_lo, _r.region_hi))
egroup(r"TMLE versus AIPW on the same samples; R-a, many strata, $m=800$")
est_cell("C6", "P-ms-H400", "Ra-TMLE", ["Q90_abs_paired_psi_std", "Q90_abs_se2_diff_std"], "R-a")
_bridge_start = len(est_rows)
egroup(r"Finite-population versus super-population target (census, $m_h=2$)")
_rb = grow("sec. Psi", "P-vanish-census", "census", "Q90_bridge")


def ftxt(f):
    """A sampling fraction with at least two decimals and no trailing zeros beyond them (0.20, 0.025, 0.0025)."""
    s = ("%.4f" % f).rstrip("0")
    return s + "0" * max(0, 2 - len(s.split(".")[1]))


for _, r in vanish.sort_values("f", ascending=False).iterrows():
    if r.terminal:
        eline(r"$f=%s$" % ftxt(r.f), SHORT["Q90_bridge"], est_int(_rb.estimate, _rb.ci_lo, _rb.ci_hi),
              region(_rb.region_lo, _rb.region_hi))
    else:
        eline(r"$f=%s$" % ftxt(r.f), SHORT["Q90_bridge"], "%s (no interval)" % num(r.q90_bridge), r"$\le 0.10$")
# Layout fix L3: keep the bridge group (heading and all its rows) on one page; "\\*" forbids a longtable page break
# after a row, so the group moves as a whole.
for _i in range(_bridge_start + 1, len(est_rows) - 1):
    if est_rows[_i].endswith("\\\\"):
        est_rows[_i] += "*"
(OUT_MAIN / "tab1_estimates.tex").write_text(
    "% GENERATED by src/displays/make_sim_displays.py from the grading outputs; do not edit by hand.\n"
    + "\n".join(est_rows) + "\n", encoding="utf-8", newline="\n")
# ----------------------------------------------------------------------------------------------
# Table 1, panel B: coverage of the estimator's nominal 95% confidence interval with its components
# ----------------------------------------------------------------------------------------------
cov_rows = []


def cov_line(setting, arm, mz, sdz, rv, cov, n):
    cov_rows.append("%s & %s & %s & %s & %s & %s & %s\\\\" % (setting, arm, num(mz, 3), num(sdz, 3), num(rv, 3),
                                                             num(cov, 3), n))


for cell, m in (("P-ms-H800", 1600), ("P-ms-H400", 800)):
    first = True
    arms = ["oracle", "Ra-AIPW"] if cell == "P-ms-H800" else ["oracle", "Ra-AIPW", "Ra-TMLE", "Rb-AIPW", "Rb-TMLE"]
    for a in arms:
        r = disp_row(ms, cell=cell, arm=a)
        cov_line(SETTING[cell] if first else "", ARM[a], r.mean_Z, r.sd_Z, r.mean_R_var, r.cover95_studentized,
                 "%d" % r.rows)
        first = False
for cell in ("P-fpc-fit-J4", "P-strict-m200"):
    first = True
    for a in ("oracle", "Ra-AIPW", "Rb-AIPW"):
        r = disp_row(e2f, cell=cell, arm=a)
        cov_line(SETTING[cell] if first else "", ARM[a], r.mean_Z, r.sd_Z, r.mean_R_var, r.cov_se_95, "%d" % r.rows)
        first = False
for cell, arms in (("P-fs-m1600", ["oracle", "Ra-AIPW"]), ("P-fs-m3200", ["oracle", "Rd-AIPW-a025", "Rd-AIPW-a030"])):
    first = True
    for a in arms:
        r = disp_row(fs, cell=cell, arm=a)
        cov_line(SETTING[cell] if first else "", ARM[a], r.mean_Z, r.sd_Z, r.mean_R_var, r.cover95_studentized,
                 "%d" % r.rows)
        first = False
(OUT_MAIN / "tab1_coverage.tex").write_text(
    "% GENERATED by src/displays/make_sim_displays.py from the grading outputs; do not edit by hand.\n"
    + "\n".join(cov_rows) + "\n", encoding="utf-8", newline="\n")
mac("SimCovSEBig", num(math.sqrt(0.95 * 0.05 / 2000), 3))
mac("SimCovSESmall", num(math.sqrt(0.95 * 0.05 / 800), 3))

# ----------------------------------------------------------------------------------------------
# Figure 1
# ----------------------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "lines.linewidth": 1.0,
                     "lines.markersize": 3.5, "pdf.fonttype": 42, "axes.spines.top": False,
                     "axes.spines.right": False})
C_OR, C_RA, C_RB, C_D25, C_D30, C_E2 = "#7f7f7f", "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"
BAND = "#dfe8f2"

fig, axes = plt.subplots(3, 3, figsize=(6.6, 6.0))
fig.subplots_adjust(left=0.075, right=0.99, top=0.965, bottom=0.13, wspace=0.36, hspace=0.9)

# (a) point limit, many strata
msx = sorted(ms.m.unique())
def ms_series(arm, col):
    d = ms[ms.arm == arm].sort_values("m")
    return d.m.values, d[col].values

panel_a = [("mean_Z", "(a1) mean $Z$", (-0.10, 0.10), "mean_Z", 0.0),
           ("sd_Z", "(a2) SD$(Z)$", (0.95, 1.05), "sd_Z_minus_1", 1.0),
           ("ks_Z", "(a3) Kolmogorov distance", (0.0, 0.05), "KS_Z", 0.0)]
for j, (col, lab, band, gstat, off) in enumerate(panel_a):
    ax = axes[0, j]
    ax.axhspan(band[0], band[1], color=BAND, lw=0, zorder=0)
    x, y = ms_series("oracle", col); ax.plot(x, y, "-", color=C_OR, marker="s", ms=2.5, label="oracle")
    x, y = ms_series("Ra-AIPW", col); ax.plot(x, y, "-", color=C_RA, marker="o", label="R-a AIPW")
    x, y = ms_series("Rb-AIPW", col); ax.plot(x, y, "--", color=C_RB, marker="o", mfc="white",
                                              label="R-b AIPW (rate condition not established)")
    r = grow("C1", "P-ms-H800", "Ra-AIPW", gstat)
    ax.errorbar([1600 * 1.08], [r.estimate + off], yerr=[[r.estimate - r.ci_lo], [r.ci_hi - r.estimate]],
                fmt="o", color=C_RA, ms=3, capsize=2, lw=0.9)
    for cell, mm, mk in (("P-fpc-fit-J4", 100, "D"), ("P-strict-m200", 200, "D")):
        r = grow("C1", cell, "Ra-AIPW", gstat)
        ax.errorbar([mm * 1.12], [r.estimate + off], yerr=[[r.estimate - r.ci_lo], [r.ci_hi - r.estimate]],
                    fmt=mk, color=C_E2, ms=3, capsize=2, lw=0.9,
                    label="R-a AIPW, $f=0.50$" if mm == 100 else None)
    ax.set_xscale("log", base=2); ax.set_xticks(msx); ax.set_xticklabels([str(v) for v in msx], rotation=0)
    ax.set_xlabel("sampled PSUs $m$"); ax.set_title(lab)
    if j == 1:
        ax.set_ylim(0.9, 1.45)
    if j == 2:
        ax.set_ylim(0, 0.13)
    if j == 0:
        ax.set_ylim(-0.12, 0.30)
h_a, l_a = axes[0, 0].get_legend_handles_labels()
y_leg = axes[0, 0].get_position().y0 - 0.068
fig.legend(h_a, l_a, loc="upper center", bbox_to_anchor=(0.53, y_leg), ncol=4, frameon=False, handlelength=1.8,
           columnspacing=1.2)

# (b) remainders, fixed strata
def fs_series(arm, col):
    d = fs[(fs.arm == arm)].sort_values("m")
    return d.m.values, d[col].values, d.t2_unresolved.values

fsx = [50, 100, 200, 400, 800, 1600, 3200]
specs_b = [
    (axes[1, 0], "q90_abs_T1_std", "(b1) $Q_{0.9}(|T_1|/\\sigma_{m,\\mathrm{wor}})$", "C4a",
     [("Ra-AIPW", C_RA, "o", "-", "R-a AIPW"), ("Rd-AIPW-a025", C_D25, "^", "-", "R-d, $\\alpha=0.25$"),
      ("Rd-AIPW-a030", C_D30, "v", "-", "R-d, $\\alpha=0.30$")]),
    (axes[1, 1], "q90_abs_T2_std", "(b2) $Q_{0.9}(|T_2|/\\sigma_{m,\\mathrm{wor}})$", "C4b",
     [("Ra-AIPW", C_RA, "o", ":", "R-a AIPW (no interval)"),
      ("Rd-AIPW-a025", C_D25, "^", "--", "R-d, $\\alpha=0.25$ (boundary)"),
      ("Rd-AIPW-a030", C_D30, "v", "-", "R-d, $\\alpha=0.30$")]),
]
# Unresolved-remainder labels (layout fix L1): each label is joined to its own point by a thin leader line and
# placed in free space, so its arm is identified by position as well as colour. Counts are unchanged.
UNRES_OFFS = {("q90_abs_T1_std", "Ra-AIPW", 50): (8, 0), ("q90_abs_T1_std", "Rd-AIPW-a025", 200): (12, 6),
              ("q90_abs_T1_std", "Rd-AIPW-a025", 400): (10, 6), ("q90_abs_T1_std", "Rd-AIPW-a030", 200): (-4, -24),
              ("q90_abs_T2_std", "Ra-AIPW", 50): (8, 0), ("q90_abs_T2_std", "Rd-AIPW-a025", 200): (-6, 1),
              ("q90_abs_T2_std", "Rd-AIPW-a025", 400): (8, 6), ("q90_abs_T2_std", "Rd-AIPW-a030", 200): (12, -18)}
for ax, col, title, claim, series in specs_b:
    ax.axhspan(1e-3, 0.10, color=BAND, lw=0, zorder=0)
    ax.axhline(0.10, color="k", lw=0.6, ls=":")
    for arm, colr, mk, ls, lab in series:
        x, y, u = fs_series(arm, col)
        ax.plot(x, y, ls, color=colr, marker=mk, label=lab)
        for xi, yi, ui in zip(x, y, u):
            if ui and ui > 0:
                offs = UNRES_OFFS.get((col, arm, int(xi)), (6, 6))
                ax.annotate("%d unres." % ui, (xi, yi), textcoords="offset points", xytext=offs, fontsize=6,
                            color=colr, ha="center" if offs[0] == 0 else ("left" if offs[0] > 0 else "right"),
                            va="bottom" if offs[1] > 0 else ("top" if offs[1] < 0 else "center"),
                            arrowprops=dict(arrowstyle="-", color=colr, lw=0.5, shrinkA=0, shrinkB=2))
        term = grades[(grades.claim == claim) & (grades.arm == arm) & grades.estimate.notna()]
        for _, r in term.iterrows():
            mm = 1600 if r.cell == "P-fs-m1600" else 3200
            ax.errorbar([mm * 1.1], [r.estimate], yerr=[[r.estimate - r.ci_lo], [r.ci_hi - r.estimate]],
                        fmt=mk, color=colr, ms=3, capsize=2, lw=0.9)
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_ylim(0.05, 2.0)
    ax.set_xticks(fsx); ax.set_xticklabels(["50", "", "200", "", "800", "", "3200"])
    ax.set_xlabel("sampled PSUs $m$"); ax.set_title(title)
ax = axes[1, 2]
for arm, colr, mk, lab in (("Ra-AIPW", C_RA, "o", "R-a, $T_1$"), ("Rd-AIPW-a025", C_D25, "^", "R-d .25, $T_1$"),
                           ("Rd-AIPW-a030", C_D30, "v", "R-d .30, $T_1$")):
    x, y, _ = fs_series(arm, "q90_sqrt_m_T1_over_rhat"); ax.plot(x, y, "-", color=colr, marker=mk, label=lab)
for arm, colr, mk, lab in (("Rd-AIPW-a025", C_D25, "^", "R-d .25, $T_2$"), ("Rd-AIPW-a030", C_D30, "v", "R-d .30, $T_2$")):
    # T2 envelopes are dashed; markers stay filled, because hollow symbols mark the forest throughout the figure
    x, y, _ = fs_series(arm, "q90_T2_over_rQ_rg"); ax.plot(x, y, "--", color=colr, marker=mk, label=lab)
ax.set_xscale("log", base=2); ax.set_xticks(fsx); ax.set_xticklabels(["50", "", "200", "", "800", "", "3200"])
ax.set_ylim(0, 2.7); ax.set_xlabel("sampled PSUs $m$")
ax.set_title("(b3) envelope ratios")
ax.legend(loc="lower left", frameon=False, ncol=2, handlelength=1.6, borderaxespad=0.1, columnspacing=0.8)

# (c) interval conservativeness across the FPC sweep
c3b = grades[grades.claim == "C3b"].copy()
fpos = {100: 0, 20: 1, 8: 2, 4: 3}
for j, (lev, z) in enumerate(((0.80, "80"), (0.90, "90"), (0.95, "95"))):
    ax = axes[2, j]
    ax.axhspan(lev - 0.02 - 1, lev - 0.02, color="#f4dede", lw=0, zorder=0)
    ax.axhline(lev, color="k", lw=0.6, ls="--")
    for J, xp in fpos.items():
        d = fpc[fpc.J == J].sort_values("population")
        jit = np.linspace(-0.22, 0.22, len(d))
        for (_, r), jj in zip(d.iterrows(), jit):
            # every population is drawn alike: coverage is reported for all 12, without a classification
            ax.plot(xp + jj, r["cov_det_sigma_m_%s" % z], "o", color=C_RA, mfc=C_RA, ms=2.6, mew=0.7)
            ax.plot(xp + jj, r["cov_wor_%s" % z], "x", color=C_OR, ms=2.4, mew=0.6)
        ax.plot([xp - 0.3, xp + 0.3], [d["overlay_%s" % z].median()] * 2, "-", color="k", lw=1.2)
    ax.set_xticks(list(fpos.values())); ax.set_xticklabels(["0.02", "0.10", "0.25", "0.50"])
    ax.set_xlabel("first-stage sampling fraction $f$")
    ax.set_title("(c%d) nominal level %.2f" % (j + 1, lev))
    ax.set_xlim(-0.5, 3.5)
    vals = np.concatenate([fpc["cov_det_sigma_m_%s" % z].values, fpc["cov_wor_%s" % z].values,
                           fpc["overlay_%s" % z].values])
    ax.set_ylim(min(lev - 0.035, vals.min() - 0.008), min(1.0, vals.max() + 0.008))
h = [plt.Line2D([], [], color=C_RA, marker="o", ls="", ms=3, label="oracle $\\pm z\\sigma_m$ (one per population)"),
     plt.Line2D([], [], color=C_OR, marker="x", ls="", ms=3, label="oracle $\\pm z\\sigma_{m,\\mathrm{wor}}$"),
     plt.Line2D([], [], color="k", lw=1.2, label="analytic overlay (median)")]
fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.53, 0.0), ncol=4, frameon=False, handlelength=1.4,
           columnspacing=1.2)
h_b = [plt.Line2D([], [], color=C_RA, marker="o", label="R-a AIPW (dotted in b2: no interval)"),
       plt.Line2D([], [], color=C_D25, marker="^", label="R-d, $\\alpha=0.25$ (dashed in b2: boundary rate)"),
       plt.Line2D([], [], color=C_D30, marker="v", label="R-d, $\\alpha=0.30$")]
y_leg_b = axes[1, 0].get_position().y0 - 0.068
fig.legend(handles=h_b, loc="upper center", bbox_to_anchor=(0.53, y_leg_b), ncol=3, frameon=False, handlelength=1.8,
           columnspacing=1.2)
fig.savefig(OUT_MAIN / "fig1_simulation.pdf", metadata={"CreationDate": None})  # no timestamp: reproducible bytes
plt.close(fig)

# ----------------------------------------------------------------------------------------------
# Supplement tables
# ----------------------------------------------------------------------------------------------
OUT_SUPP.mkdir(parents=True, exist_ok=True)
OUT_REC.mkdir(parents=True, exist_ok=True)
HDR = "% GENERATED by src/displays/make_sim_displays.py from the grading outputs; do not edit by hand.\n"
# The tables below are written twice where needed: in their original form, with the grading labels, to the archive
# record (out=OUT_REC), and in the supplement's form, without them (the default, OUT_SUPP). Both read the same data.


def longtable(fname, colspec, header, body, caption, label, size=r"\footnotesize", land=False, out=None):
    txt = [HDR] + (["\\begin{landscape}"] if land else []) + ["{%s\\setlength{\\tabcolsep}{3pt}" % size,
           "\\begin{longtable}{%s}" % colspec,
           "\\caption{%s}\\label{%s}\\\\" % (caption, label),
           "\\toprule", header + "\\\\", "\\midrule", "\\endfirsthead",
           "\\multicolumn{%d}{@{}l}{\\textit{(continued)}}\\\\" % colspec.count("l") if False else "",
           "\\toprule", header + "\\\\", "\\midrule", "\\endhead",
           "\\bottomrule", "\\endlastfoot"]
    txt = [t for t in txt if t != ""]
    txt += body
    txt += ["\\end{longtable}}"] + (["\\end{landscape}"] if land else []) + [""]
    ((out or OUT_SUPP) / fname).write_text("\n".join(txt), encoding="utf-8", newline="\n")


# S: every frozen cell
ROLE = {
    "P-ms-H400": ("E1", "many strata, $H=400$, $m_h=2$, $m=800$", "C1/C2 (R-b), C6; display"),
    "P-ms-H800": ("E1", "many strata, $H=800$, $m_h=2$, $m=1600$", "C1, C2 (R-a)"),
    "P-ms-H25": ("E1", "many strata, $H=25$, $m=50$", "trajectory display"),
    "P-ms-H50": ("E1", "many strata, $H=50$, $m=100$", "trajectory display"),
    "P-ms-H100": ("E1", "many strata, $H=100$, $m=200$", "trajectory display"),
    "P-ms-H200": ("E1", "many strata, $H=200$, $m=400$", "trajectory display"),
    "P-fs-m50": ("E1", "fixed strata, $H=5$, $m=50$", "remainder trajectory"),
    "P-fs-m100": ("E1", "fixed strata, $m=100$", "remainder trajectory"),
    "P-fs-m200": ("E1", "fixed strata, $m=200$", "remainder trajectory"),
    "P-fs-m400": ("E1", "fixed strata, $m=400$", "remainder trajectory"),
    "P-fs-m800": ("E1", "fixed strata, $m=800$", "remainder trajectory"),
    "P-fs-m1600": ("E1", "fixed strata, $m=1600$", "C4a (R-a)"),
    "P-fs-m3200": ("E1", "fixed strata, $m=3200$ (R-d only)", "C4a, C4b, C5 (R-d)"),
    "P-fpc-or-J100": ("E2", "oracle, $H=50$, $J=100$, $f=0.02$", "C3b"),
    "P-fpc-or-J20": ("E2", "oracle, $J=20$, $f=0.10$", "C3b"),
    "P-fpc-or-J8": ("E2", "oracle, $J=8$, $f=0.25$", "C3b"),
    "P-fpc-or-J4": ("E2", "oracle, $J=4$, $f=0.50$", "C3a, C3b"),
    "P-fpc-fit-J4": ("E2", "fitted, $H=50$, $J=4$, $f=0.50$, $m=100$", "C1"),
    "P-strict-m200": ("E2", "fitted, $H=100$, $J=4$, $f=0.50$, $m=200$", "C1"),
    "P-calib-J4": ("E2", "analytic calibration, $f=0.50$, population 1", "implementation check"),
    "P-wr-calib": ("E2", "with-replacement code check", "implementation check"),
    "P-strict-census": ("E2", "fixed-fraction trajectory, census, $H=25,\\dots,400$", "strictness witness"),
    "P-unequal-census": ("E2", "unequal-allocation pair, census", "invariance check"),
    "P-vanish-census": ("E2", "vanishing fraction, census, 5 rungs", "bridge to $\\Psi$"),
    "P-e3-rc-m400": ("E3", "interpolating-forest comparison, $m=400$", "reported only"),
    "P-e3-vsens-V2": ("E3", "$V=2$ sensitivity, $m=200$", "reported only"),
    "P-e3-vsens-V5": ("E3", "$V=5$ sensitivity, $m=200$", "reported only"),
    "P-e3-ctrl-m400": ("E3", "control conditions, $m=400$", "power controls"),
    "P-e4-array": ("E4", "deterministic array, $H=25,\\dots,400$", "C7"),
}
for H, m in ((25, 50), (50, 100), (100, 200), (200, 400), (400, 800), (800, 1600)):
    ROLE["P-ms-or-H%d" % H] = ("E1", "many strata oracle, $H=%d$, $m=%d$" % (H, m),
                                "C1 $Z_{\\mathrm{oracle}}$, C2 oracle limb" if H == 800 else "oracle display")
for m in (50, 100, 200, 400, 800, 1600):
    ROLE["P-fs-or-m%d" % m] = ("E1", "fixed strata oracle, $m=%d$" % m,
                               "C1 $Z_{\\mathrm{oracle}}$" if m == 1600 else "oracle display")
body = []
rq = reqchk.set_index("cell")
for _, r in complete.iterrows():
    exp, desc, role = ROLE[r.cell]
    q = rq.loc[r.cell]
    body.append("%s & %s & %s & %s & %s & %s & %s\\\\" % (
        tex_escape(r.cell), exp, desc, role, tex_escape(r.detail).replace(" x ", r" $\times$ "),
        "{:,}".format(int(q.required_checks)).replace(",", "{,}"), int(q.failed)))
assert set(ROLE) == set(complete.cell), set(ROLE) ^ set(complete.cell)
longtable("tabS_cells.tex", "@{}llp{0.25\\linewidth}p{0.15\\linewidth}p{0.22\\linewidth}rr@{}",
          "cell & exp. & coordinates & role & design (all COMPLETE) & checks & failed", body,
          "Every frozen production cell, its role, its completed design and its required implementation checks. "
          "All %s cells are complete against the frozen design. One score-total equality check "
          "failed in one replicate of the with-replacement code-check cell; its variance comparison "
          "passed. No reported benchmark comparison uses that cell." % len(complete),
          "tab:S-cells", size=r"\scriptsize", land=True, out=OUT_REC)

# S: every graded scientific row
body = []
PC_SHORT = {"C1 shift": r"C1 control: shifted centre", "C2 wrong clustering: unit-level": "C2 control: units as clusters",
            "C2 wrong clustering: no m_h": r"C2 control: no $m_h/(m_h-1)$",
            "C2 wrong clustering: no stratification": "C2 control: strata ignored",
            "C4 aligned error": "C4 control: aligned error"}


def arm_short(a):
    for k, v in PC_SHORT.items():
        if str(a).startswith(k):
            return v
    if str(a).startswith("R-d (alpha"):
        return r"R-d, $\alpha<1/4$"
    return ARM.get(a, tex_escape(a))


sci = grades[~grades.claim.isin(["C3b", "C3b premise: Z_oracle", "C1 Z_oracle", "C2 oracle limb", "C3a", "C7"])]
for _, r in sci.iterrows():
    d = 4 if ("paired" in r.statistic or "se2" in r.statistic) else 3
    body.append("%s & %s & %s & %s & %s & %s & %s\\\\" % (
        tex_escape(r.claim), tex_escape(r.cell), arm_short(r.arm), SHORT.get(r.statistic, STAT.get(r.statistic, tex_escape(r.statistic))),
        num(r.estimate, d), interval(r.ci_lo, r.ci_hi, d), short_disp(r.disposition)))
longtable("tabS_grades.tex", "@{}lllp{0.14\\linewidth}rll@{}",
          "claim & cell & arm & statistic & estimate & 95\\% MC interval & disposition", body,
          "Every graded scientific, boundary and power-control row of the grading pass, with its 95\\% Monte Carlo "
          "interval (percentile interval of 10{,}000 keyed population-block resamples of the saved output; "
          "population resamples for the census bridge). A dash marks a row for which no interval was formed, "
          "because a prerequisite was not established or the trajectory was stopped before production. Reasons "
          "for UNGRADABLE rows are given in Table~\\ref{tab:S-dispositions}, the text and Table~\\ref{tab:S-gates}. "
          "$Q_{0.9}|T_j|$ is the 0.9 quantile of $|T_j|/\\sigma_{m,\\mathrm{wor}}$; $\\Delta\\hat\\psi="
          "(\\hat\\psi_{\\mathrm{TMLE}}-\\hat\\psi_{\\mathrm{AIPW}})/\\sigma_{m,\\mathrm{wor}}$; "
          "$\\Delta\\widehat{\\mathrm{se}}^2=(\\widehat{\\mathrm{se}}^2_{\\mathrm{TMLE}}-\\widehat{\\mathrm{se}}^2_{"
          "\\mathrm{AIPW}})/\\sigma_m^2$; the bridge is standardized by $\\sigma_{m,\\mathrm{wor}}$; KS is the "
          "Kolmogorov distance.",
          "tab:S-grades", size=r"\scriptsize", land=True, out=OUT_REC)

# S: the original dispositions in the analysis plan's order (formerly the main text's Table 1, panel A)
longtable("tabS_dispositions.tex",
          "@{}>{\\raggedright\\arraybackslash}p{0.25\\linewidth}lll>{\\raggedright\\arraybackslash}p{0.19\\linewidth}@{}",
          "cell & statistic & estimate [95\\% MC interval] & region & disposition", _disp_tab,
          "Original dispositions under the analysis plan's classification rule, in the plan's order (C1--C7), "
          "followed by the bridge to the super-population target and the power controls; the labels are those of "
          "the plan and the grading output, unchanged. Each graded statistic is given with its point estimate and "
          "95\\% Monte Carlo interval (percentile interval of 10{,}000 resamples of whole populations; of samples "
          "for the population-indexed oracle cells; of populations for the census bridge), the compliant region "
          "of its margin, and its disposition. PASS: the interval lies inside the region. RED-FLAG: the interval "
          "lies entirely outside it. UNGRADABLE: the interval straddles the boundary (MC precision), or a "
          "prerequisite was not established before grading (premise), or the trajectory was stopped before "
          "production (stopped), or no rung below the quarter rate was retained (no rung). BOUNDARY: for the "
          "product remainder $T_2$ (C4b), $\\alpha=0.25$ lies on the quarter-rate threshold and is reported without "
          "a grade; $T_1$ at that rate is graded (C4a). Population-indexed "
          "rows give counts of populations, or of array rungs for C7, as PASS/UNGRADABLE; for C3b they are given "
          "at the levels 0.80, 0.90 and 0.95 in turn. Arms: R-a, correctly specified random-intercept learner; "
          "R-b, random forest; R-d, R-a refitted on a throttled subset of training PSUs, designed so that its "
          "error decays at rate $m^{-\\alpha}$. KS is the Kolmogorov distance of $Z$ from $\\mathcal N(0,1)$; "
          "$Q_{0.9}|\\cdot|$ is the 0.9 quantile of the absolute value of the quantity standardized by "
          "$\\sigma_{m,\\mathrm{wor}}$ ($\\sigma_m^2$ for the difference in variance estimates); the gap error is "
          "$\\widehat{\\mathrm{Var}}(S_{\\mathrm{des}})/\\sigma_m^2-(1-\\Delta_m/\\sigma_m^2)$. "
          "Table~\\ref{tab:S-crosswalk} gives the plain-language reading used in the main text.",
          "tab:S-dispositions", size=r"\scriptsize", out=OUT_REC)

# S: population-indexed summary
body = []
for _, r in popsum.iterrows():
    body.append("%s & %s & %s & %s & %d & %d & %d & %d\\\\" % (
        tex_escape(r.claim).replace("Z\\_oracle", "$Z_{\\mathrm{oracle}}$"), tex_escape(r.role), tex_escape(r.cell),
        STAT.get(r.statistic, tex_escape(r.statistic)).replace("coverage\\_minus\\_", "coverage $-$ "),
        r.populations, r.PASS, r.UNGRADABLE, r.RED_FLAG))
longtable("tabS_popindexed.tex", "@{}llllrrrr@{}",
          "claim & role & cell & statistic & units & PASS & UNGRADABLE & RED-FLAG", body,
          "Population-indexed dispositions of the design-conditional cells (oracle estimator; each of 12 "
          "populations, or each of the 5 array rungs for C7, graded on its own redraws and never averaged). "
          "A C3b row is UNGRADABLE when that population's own $Z_{\\mathrm{oracle}}$ does not pass all three C1 "
          "margins; every C3b UNGRADABLE here has that premise reason, arising from Monte Carlo straddles of the "
          "$Z_{\\mathrm{oracle}}$ margins. No row is BLOCKED or NOT GRADED.", "tab:S-pop", size=r"\scriptsize", land=True, out=OUT_REC)

# S: many-strata trajectories
body = []
for _, r in ms.sort_values(["m", "arm"]).iterrows():
    body.append("%d & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        r.m, ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.ks_Z), num(r.cover95_studentized),
        num(r.mean_R_var), num(r.sd_R_var), num(r.transfer_mean)))
longtable("tabS_ms.tex", "@{}rlrrrrrrr@{}",
          "$m$ & arm & mean $Z$ & SD$(Z)$ & KS & cov.\\ 95 & mean $R_{\\mathrm{var}}$ & SD$(R_{\\mathrm{var}})$ & transfer",
          body,
          "Many-strata regime ($m_h=2$, $J=20$, $f=0.10$), every arm and rung; 400 populations $\\times$ 5 draws, "
          "populations equally weighted. $Z=(\\hat\\psi-\\Psi_N)/\\sigma_{m,\\mathrm{wor}}$, so mean $Z$ is the "
          "standardized bias; cov.\\ 95 is the coverage of the estimator's nominal 95\\% confidence interval "
          "$\\hat\\psi\\pm1.96\\,\\widehat{\\mathrm{se}}$ about $\\Psi_N$; $R_{\\mathrm{var}}=\\widehat{\\mathrm{se}}^2/"
          "\\sigma_m^2$. Graded rungs are $m=1600$ (R-a) and $m=800$ (R-b, TMLE); other rungs are descriptive. "
          "No fit failed.", "tab:S-ms", size=r"\scriptsize", out=OUT_REC)

# S: fixed-strata trajectories
body = []
for _, r in fs.sort_values(["m", "arm"]).iterrows():
    tmle = isinstance(r.remainder_note, str)
    body.append("%d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        r.m, ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.cover95_studentized), num(r.mean_R_var),
        "excluded" if tmle else ("---" if r.arm == "oracle" else "%d" % r.t2_unresolved),
        num(r.q90_abs_T1_std), num(r.q90_abs_T2_std), num(r.q90_sqrt_m_T1_over_rhat, 2),
        num(r.q90_T2_over_rQ_rg, 2), num(r.mean_r_hat)))
longtable("tabS_fs.tex", "@{}rlrrrrrrrrrr@{}",
          "$m$ & arm & mean $Z$ & SD$(Z)$ & cov.\\ 95 & $R_{\\mathrm{var}}$ & unres. & $Q_{0.9}|T_1|$ & "
          "$Q_{0.9}|T_2|$ & env.\\ $T_1$ & env.\\ $T_2$ & $\\hat r_m$", body,
          "Fixed-strata regime ($H=5$, $f=0.10$), every arm and rung; 400 populations $\\times$ 5 draws. "
          "$Q_{0.9}|T_j|$ is the 0.9 quantile of $|T_j|/\\sigma_{m,\\mathrm{wor}}$ over replicates with a resolved "
          "remainder; ``unres.'' counts replicates whose $T_2$ (and hence $T_1$) could not be resolved "
          "numerically and are kept as missing. The envelope ratios $\\sqrt m|T_1|/\\hat r_m$ and "
          "$|T_2|/(\\hat r_Q\\hat r_g)$ (0.9 quantiles) and the mean out-of-fold error $\\hat r_m$ are displays "
          "without margins. R-a TMLE remainders were not independently checked and are shown as "
          "dashes. R-d runs at $m\\ge200$ only; $\\alpha=0.20$ and $0.35$ were stopped before production.",
          "tab:S-fs", size=r"\scriptsize", land=True, out=OUT_REC)

# S: E2 fitted cells
body = []
for _, r in e2f.iterrows():
    body.append("%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        SETTING[r.cell], ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.ks_Z), num(r.mean_R_var),
        num(r.cov_se_80), num(r.cov_se_90), num(r.cov_se_95)))
longtable("tabS_e2fit.tex", "@{}llrrrrrrr@{}",
          "setting & arm & mean $Z$ & SD$(Z)$ & KS & mean $R_{\\mathrm{var}}$ & cov.\\ 80 & cov.\\ 90 & cov.\\ 95",
          body,
          "Fitted cells at first-stage fraction $f=0.50$ (160 populations $\\times$ 5 draws). Coverage is that of "
          "the estimator's nominal confidence interval $\\hat\\psi\\pm z\\,\\widehat{\\mathrm{se}}$ about $\\Psi_N$; "
          "it is reported beside its components, the standardized bias (mean $Z$), the dispersion relative to the "
          "without-replacement norming (SD$(Z)$) and the variance ratio against the ultimate-cluster target "
          "($R_{\\mathrm{var}}$).", "tab:S-e2fit", size=r"\scriptsize", out=OUT_REC)

# S: E2 FPC oracle, per population
body = []
for _, r in fpc.sort_values(["J", "population"], ascending=[False, True]).iterrows():
    disp = c3b[(c3b.cell == "P-fpc-or-J%d" % r.J) & (c3b.population == r.population)]
    ok = (disp.disposition == "PASS").all()
    body.append("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        JF[int(r.J)], int(r.population), num(r.sigma_ratio), num(r.cov_wor_95), num(r.cov_det_sigma_m_80),
        num(r.cov_det_sigma_m_90), num(r.cov_det_sigma_m_95), num(r.overlay_95), num(r.cov_se_95),
        num(r.c3a_statistic), "PASS" if ok else "UNGR."))
longtable("tabS_fpc.tex", "@{}rrrrrrrrrrl@{}",
          "$f$ & pop. & $\\sigma_m/\\sigma_{m,\\mathrm{wor}}$ & wor 95 & $\\sigma_m$ 80 & $\\sigma_m$ 90 & "
          "$\\sigma_m$ 95 & overlay 95 & $\\widehat{\\mathrm{se}}$ 95 & C3a stat. & C3b", body,
          "Oracle interval hierarchy across the first-stage fraction sweep ($H=50$, $m_h=2$, $m=100$), each of 12 "
          "populations on 2{,}000 redraws (10{,}000 at $f=0.50$ for the C3a statistic). Columns give coverage of "
          "$\\hat\\psi_{\\mathrm{oracle}}\\pm z\\sigma_{m,\\mathrm{wor}}$ (``wor''), of "
          "$\\hat\\psi_{\\mathrm{oracle}}\\pm z\\sigma_m$ (the graded C3b interval) at 0.80, 0.90, 0.95, the analytic "
          "overlay $2\\Phi(z\\sigma_m/\\sigma_{m,\\mathrm{wor}})-1$ at 0.95, and coverage of "
          "$\\hat\\psi_{\\mathrm{oracle}}\\pm z\\widehat{\\mathrm{se}}$. C3b is PASS or UNGRADABLE (premise) per "
          "population, identically at the three levels.", "tab:S-fpc", size=r"\scriptsize", out=OUT_REC)

# S: E1 oracle anchors, summarised by cell
body = []
for cell, d in e1o.groupby("cell", sort=False):
    body.append("%s & %d & %s & %s & %s & %s & %s & %s\\\\" % (
        tex_escape(cell), len(d),
        "%s to %s" % (num(d.mean_Z.min()), num(d.mean_Z.max())), "%s to %s" % (num(d.sd_Z.min()), num(d.sd_Z.max())),
        "%s to %s" % (num(d.ks_Z.min()), num(d.ks_Z.max())), "%s to %s" % (num(d.mean_R_var.min()), num(d.mean_R_var.max())),
        "%s to %s" % (num(d.cover95_studentized.min()), num(d.cover95_studentized.max())),
        "%s to %s" % (num(d.linearization_monitor.min()), num(d.linearization_monitor.max()))))
longtable("tabS_e1oracle.tex", "@{}lrllllll@{}",
          "cell & pop. & mean $Z$ & SD$(Z)$ & KS & mean $R_{\\mathrm{var}}$ & cov.\\ 95 & lin.\\ monitor", body,
          "Design-conditional oracle anchors of E1: range over the 12 populations (4{,}000 redraws each) of the "
          "oracle $Z$ summaries, the variance ratio, the coverage of $\\hat\\psi_{\\mathrm{oracle}}\\pm1.96\\,"
          "\\widehat{\\mathrm{se}}$ and the linearization monitor $\\widehat{\\mathrm{Var}}(S_{\\mathrm{des}})/"
          "\\sigma^2_{m,\\mathrm{wor}}$. The terminal rungs (many strata $H=800$; fixed strata $m=1600$) are graded "
          "per population as implementation checks (Table~\\ref{tab:S-pop}); the others are displays.",
          "tab:S-e1oracle", size=r"\scriptsize", land=True, out=OUT_REC)

# S: E2 census limbs
body = []
for _, r in strict.iterrows():
    body.append("strict, $f=0.50$ & $H=%d$ & %s & %s & %s & %s\\\\" % (
        r.H, num(r.min_Delta_over_sigma2_m), num(r.median_Delta_over_sigma2_m), num(r.max_Delta_over_sigma2_m),
        num(r.median_sigma_ratio)))
longtable("tabS_strict.tex", "@{}llrrrr@{}",
          "trajectory & rung & min $\\Delta_m/\\sigma_m^2$ & median & max & median $\\sigma_m/\\sigma_{m,\\mathrm{wor}}$",
          body,
          "Strictness witness (census algebra, 120 populations per rung, no sampling): the relative "
          "without-replacement gap along the fixed-fraction trajectory $m_h=2$, $J=4$. It stays away from zero "
          "along the frozen rungs, which is compatible with a strict regime; a finite trajectory cannot establish "
          "an asymptotic limit infimum.", "tab:S-strict", size=r"\scriptsize", out=OUT_REC)
body = []
for _, r in uneq.iterrows():
    body.append("%s & %d & %s & %s\\\\" % (tex_escape(r.alloc).replace("\\_", "-"), r.populations,
                                         num(r.median_Delta_over_sigma2_m), "\\ensuremath{%.4e}" % r.median_Delta_m))
longtable("tabS_unequal.tex", "@{}lrrr@{}",
          "allocation & populations & median $\\Delta_m/\\sigma_m^2$ & median $\\Delta_m$", body,
          "Unequal-allocation invariance check (census, same populations under the two allocations): "
          "$\\Delta_m$ is allocation-invariant and agrees between the two allocations, while the relative gap "
          "differs because $\\sigma_m^2$ does. An implementation check.", "tab:S-unequal", size=r"\scriptsize", out=OUT_REC)
body = []
for _, r in vanish.sort_values("f", ascending=False).iterrows():
    body.append("%d & %d & %s & %s & %s & %s & %s%s\\\\" % (
        r.H, r.J, num(r.f, 4), num(r.median_bridge), num(r.q90_bridge), num(r.max_bridge),
        num(r.ratio_mean_Delta_over_var_Smod, 2), " (graded rung)" if r.terminal else ""))
longtable("tabS_vanish.tex", "@{}rrrrrrl@{}",
          "$H$ & $J$ & $f$ & median & $Q_{0.9}$ & max & $\\bar\\Delta_m/\\widehat{\\mathrm{Var}}(\\Psi_N-\\Psi)$",
          body,
          "Vanishing-fraction trajectory of the bridge $|\\Psi_N-\\Psi|/\\sigma_{m,\\mathrm{wor}}$ (census, 120 "
          "populations per rung, $m_h=2$). Only the terminal rung $(H,J)=(200,800)$ is graded; the others are "
          "descriptive. The last column reports the leading-order relation between the mean without-replacement "
          "gap and the across-population variance of $\\Psi_N-\\Psi$, which is outside Theorem~\\ref{thm:clt} and "
          "carries no criterion.", "tab:S-vanish", size=r"\scriptsize", out=OUT_REC)

# S: E3
body = []
for _, r in rc.iterrows():
    body.append("%s & %s & %s & %s & %s & %s & %s\\\\" % (
        "R-c comparison, $m=400$", ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.cover95_studentized), num(r.mean_R_var),
        "reported only" if r.arm.startswith("Rc") else "diagnostic"))
for _, r in vs.iterrows():
    body.append("%s & %s & %s & %s & %s & %s & %s\\\\" % (
        "$V=%s$, $m=200$" % r.cell[-1], ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.cover95_studentized),
        num(r.mean_R_var), "diagnostic"))
longtable("tabS_e3.tex", "@{}llrrrrl@{}",
          "cell & arm & mean $Z$ & SD$(Z)$ & cov.\\ 95 & mean $R_{\\mathrm{var}}$ & status", body,
          "E3 displays (fixed strata, 400 populations $\\times$ 5 draws): the interpolating-forest comparison R-c under "
          "PSU-level folds, no cross-fitting, unit-level folds and unit-level folds with $C_{\\mathrm{psu}}$ "
          "coarsened, beside the R-a controls on the same samples; and the fold-number sensitivity ($V=2$ versus "
          "$V=5$, random forest R-b). R-c's quarter-rate premise is not established, so its rows carry no grade. "
          "E3 carries no C2 grade.", "tab:S-e3", size=r"\scriptsize", land=True, out=OUT_REC)

# S: E4
body = []
for _, r in e4.iterrows():
    body.append("%d & %d & %s & %s & %s & %s & %s & %s\\\\" % (
        r.H, r.draws, num(r.mean_Z), num(r.sd_Z), num(r.ks_Z), num(r.var_ratio_wor), num(r.mean_se2_over_sigma2_m),
        num(r.lindeberg, 4)))
longtable("tabS_e4.tex", "@{}rrrrrrrr@{}",
          "$H$ & draws & mean $Z$ & SD$(Z)$ & KS & $\\widehat{\\mathrm{Var}}/\\sigma^2_{m,\\mathrm{wor}}$ & "
          "mean $\\widehat{\\mathrm{se}}^2/\\sigma_m^2$ & Lindeberg ratio", body,
          "E4 deterministic triangular array (oracle only, one unit per PSU, no outcome or treatment model): the "
          "design algebra on a population for which no super-population law is posited. The Lindeberg ratio is "
          "$\\max_h\\mathrm{Var}_{\\mathrm{des}}(W_h)/\\sum_h\\mathrm{Var}_{\\mathrm{des}}(W_h)$. An implementation "
          "check, never evidence for or against a theorem claim.", "tab:S-e4", size=r"\scriptsize", out=OUT_REC)

# S: gates
body = []
for _, r in g5.iterrows():
    body.append("G5 rate validation & R-d, $\\alpha=%.2f$ & slope %s, 95\\%% interval %s & %s\\\\" % (
        r.alpha, num(r.slope), interval(r.lo, r.hi), "PASS" if r.status == "PASS" else "failed; trajectory stopped"))
r = rbg.iloc[0]
body.append("R-b rate gate & many strata & slope %s, one-sided 95\\%% upper bound %s & "
            "$L^2$ established; quarter rate not established\\\\" % (num(r.slope), num(r.upper_95)))
# Readable names for the control conditions (the raw g6.csv strings are plain-text notation).
G6NAME = {
    "C1 shift": r"C1 shift: mean $Z$ about $\Psi_N+0.30\,\sigma_{m,\mathrm{wor}}$ (R8)",
    "C2 wrong clustering: no m_h": r"C2 wrong clustering: no $m_h/(m_h-1)$ factor, mean $R_{\mathrm{var}}$",
    "C2 wrong clustering: unit-level": r"C2 wrong clustering: units as clusters, mean $R_{\mathrm{var}}$",
    "C2 wrong clustering: no stratification": r"C2 wrong clustering: strata ignored, mean $R_{\mathrm{var}}$",
    "C4 aligned error": r"C4 aligned error: $Q_{0.9}(|T_2|/\sigma_{m,\mathrm{wor}})$ ($\kappa_{\mathrm{AE}}$ fixed; C-AE-3)",
}


def g6_name(c):
    for k, v in G6NAME.items():
        if str(c).startswith(k):
            return v
    raise KeyError(c)


for _, r in g6.iterrows():
    if str(r.control).startswith("U10"):
        continue
    body.append("G6 limb (a) & %s & predicted %s, distance %s, MCSE %s & %s\\\\" % (
        g6_name(r.control),
        num(float(r.predicted)), num(float(r.distance)), num(float(r.mcse), 4),
        "graded power control" if bool(r.graded) else "diagnostic only"))
body.append("$\\kappa_{\\mathrm{AE}}$ calibration & 160 panel populations & $\\kappa_{\\mathrm{AE}}=%s$, "
            "$Q_{0.9}(T_2/\\sigma_{m,\\mathrm{wor}})=%s$ & frozen, reused unchanged\\\\" % (
                num(kappa["kappa"], 4), num(kappa["stat"], 2)))
longtable("tabS_gates.tex", "@{}p{0.15\\linewidth}p{0.27\\linewidth}p{0.32\\linewidth}p{0.2\\linewidth}@{}",
          "gate & object & value & consequence", body,
          "Pre-run and dry-run gates that decide premises before grading. G5 checks that each throttled learner's "
          "out-of-fold error decays at its designed rate $m^{-\\alpha}$ (slope of $\\log\\hat r_m$ on $\\log m$ over "
          "$m=200,\\dots,3200$; population-clustered CR2 standard errors, $t_{24}$ reference); a failure stops the "
          "trajectory. The random-forest gate decides R-b's premises on its own ladder $m=50,\\dots,800$. G6 limb "
          "(a) admits a power control only if its census-predicted effect clears its margin by more than three "
          "Monte Carlo standard errors.", "tab:S-gates", size=r"\scriptsize", out=OUT_REC)

# S: G7 premise audit
body = []
# The raw g7.csv strings use plain-text notation (and one mis-encoded section sign); render them as LaTeX.
G7REP = [(r"cond:mix", r"Condition~\ref{cond:mix}"), (r"cond:nd", r"Condition~\ref{cond:nd}"),
         (r"eta0", r"$\eta_0$"), (r"corr\(M_k, phibar_k \| stratum\)", r"$\mathrm{corr}(M_k,\bar\phi_k\mid\mathrm{stratum})$"),
         (r">=", r"$\ge$"), (r"DEFF_aug", r"DEFF$_{\mathrm{aug}}$"), (r"lambda_min", r"$\lambda_{\min}$"),
         (r"c_H", r"$c_H$"), (r"E1-E3", r"E1--E3"), (r"the .12\.3 generating code", r"the generating code, specification \\S12.3")]


def g7_tex(s):
    s = str(s)
    for a, b in G7REP:
        s = re.sub(a, lambda _m: b.replace("\\\\", "\\"), s)
    return s.replace("%", r"\%")


for _, r in g7.iterrows():
    body.append("%s & %s\\\\" % (g7_tex(r.condition), g7_tex(r.entry)))
longtable("tabS_g7.tex", "@{}p{0.48\\linewidth}p{0.48\\linewidth}@{}", "condition & entry", body,
          "The conditions-to-data-generating-process audit (gate G7), as recorded by the round-3 assembly and "
          "read by the grader. Entries are QUANTITATIVE (a recorded value against its criterion), SATISFIED BY "
          "CONSTRUCTION (a specification clause with the code check that enforces it) or NOT POSITED / NOT "
          "INVOKED.", "tab:S-g7", size=r"\scriptsize", out=OUT_REC)

# S: unresolved T2
body = []
for (cell, arm), v in unresolved.items():
    body.append("%s & %s & %d & %s\\\\" % (tex_escape(cell), ARM[arm], v,
                                          "excluded" if arm.endswith("TMLE") else "display rung; kept as missing"))
longtable("tabS_unresolved.tex", "@{}llrl@{}", "cell & arm & unresolved & use", body,
          "Arm-by-replicate entries whose second-order remainder $T_2$ could not be resolved by the refined "
          "quadrature (production audit; %s entries in total, all at $m\\le400$). A sample can contribute an entry "
          "for more than one arm, so this is not a count of distinct samples; the R-a TMLE rows are already "
          "excluded from remainder reporting because they were not independently checked. $T_2$ and $T_1$ are "
          "recorded as missing and every estimate is "
          "kept; no graded rung is affected." % sum(unresolved.values()),
          "tab:S-unresolved", size=r"\scriptsize", out=OUT_REC)

# ----------------------------------------------------------------------------------------------
# Supplement versions (OUT_SUPP): the same values, without the grading labels, with plain captions
# ----------------------------------------------------------------------------------------------
# Where each cell is reported; "record" = the extended tables only. A* labels are the supplement's (xr prefix in the
# cell index); T1/F1 are the main text's Table 1 and Figure 1.
PURPOSE = {
    "P-ms-H400": "forest and TMLE at their largest size", "P-ms-H800": "normal approximation and variance (R-a)",
    "P-fs-m1600": "remainders, correctly specified learner", "P-fs-m3200": "remainders, throttled learner",
    "P-fpc-or-J4": "oracle variance gap and coverage", "P-fpc-fit-J4": "fitted estimator at $f=0.50$",
    "P-strict-m200": "fitted estimator at $f=0.50$", "P-calib-J4": "check of the analytic normings",
    "P-wr-calib": "with-replacement code check", "P-strict-census": "relative gap at a fixed fraction",
    "P-unequal-census": "allocation-invariance check", "P-vanish-census": r"bridge from $\Psi_N$ to $\Psi$",
    "P-e3-rc-m400": "interpolating forest and fold schemes", "P-e3-vsens-V2": "number of folds",
    "P-e3-vsens-V5": "number of folds", "P-e3-ctrl-m400": "control conditions",
    "P-e4-array": "design computation without a generating model",
}
WHERE = {
    "P-ms-H400": ["S-ms", "T1"], "P-ms-H800": ["S-ms", "T1", "F1"], "P-fs-m1600": ["S-fs", "T1", "F1"],
    "P-fs-m3200": ["S-fs", "T1", "F1"], "P-fpc-fit-J4": ["S-e2fit", "T1", "F1"], "P-strict-m200": ["S-e2fit", "T1", "F1"],
    "P-strict-census": ["S-strict"], "P-vanish-census": ["S-vanish", "T1"], "P-e3-rc-m400": ["S-e3"],
    "P-e3-vsens-V2": ["S-e3"], "P-e3-vsens-V5": ["S-e3"], "P-e3-ctrl-m400": ["S-estimates", "S-rates"],
    "P-calib-J4": ["record"], "P-wr-calib": ["record"], "P-unequal-census": ["record"], "P-e4-array": ["record"],
}
for c in complete.cell:
    if c.startswith("P-ms-H") or c.startswith("P-fs-m"):
        PURPOSE.setdefault(c, "trajectory")
        WHERE.setdefault(c, ["S-ms" if c.startswith("P-ms") else "S-fs", "F1"])
    elif c.startswith("P-fpc-or-J"):
        PURPOSE.setdefault(c, "oracle coverage by sampling fraction")
        WHERE.setdefault(c, ["S-fpc", "F1"])
    elif c.startswith("P-ms-or-H") or c.startswith("P-fs-or-m"):
        PURPOSE.setdefault(c, "oracle at 12 fixed populations")
        WHERE.setdefault(c, ["record"])
WHERE["P-fpc-or-J4"] = ["S-fpc", "T1", "F1"]
assert set(PURPOSE) == set(complete.cell) == set(WHERE), set(PURPOSE) ^ set(complete.cell)


def rec_table(cell):
    """The record-only table that reports a cell (labels as in the record)."""
    if cell.startswith("P-ms-or-H") or cell.startswith("P-fs-or-m"):
        return "S-e1oracle"
    return {"P-unequal-census": "S-unequal", "P-e4-array": "S-e4"}.get(cell, "S-cells")


def where_txt(cell, pre=""):
    out = []
    for k in WHERE[cell]:
        if k == "T1":
            out.append(r"Table~\ref{tab:sim}")
        elif k == "F1":
            out.append(r"Figure~\ref{fig:sim}")
        elif k == "record":
            out.append(r"code repository" if not pre else r"Table~\ref{tab:%s}" % rec_table(cell))
        else:
            out.append(r"Table~\ref{%stab:%s}" % (pre, k))
    return ", ".join(out)


def cell_rows(pre=""):
    rows_ = []
    for _, r in complete.iterrows():
        exp, desc, _role = ROLE[r.cell]
        rows_.append("%s & %s & %s & %s & %s & %s\\\\" % (
            tex_escape(r.cell), exp, desc.replace("rungs", "settings"), PURPOSE[r.cell],
            re.sub(r"^1 populations", "1 population",
                   tex_escape(r.detail).replace(" x ", r" $\times$ ").replace("frozen ", "")
                   .replace("rungs", "settings")), where_txt(r.cell, pre)))
    return rows_


CELLS_HDR = "cell & exp. & coordinates & question & design & reported in"
CELLS_SPEC = ("@{}l@{\\hspace{3pt}}l>{\\raggedright\\arraybackslash}p{0.19\\linewidth}"
              ">{\\raggedright\\arraybackslash}p{0.17\\linewidth}>{\\raggedright\\arraybackslash}p{0.16\\linewidth}"
              ">{\\raggedright\\arraybackslash}p{0.16\\linewidth}@{}")
longtable("tabS_cells.tex", CELLS_SPEC, CELLS_HDR, cell_rows(),
          "The %s simulation cells: experiment, coordinates, the question each answers, its Monte Carlo design and "
          "where its results are reported. Every cell completed, and no learner fit failed. Implementation checks "
          "and the oracle cells at fixed populations are reported in the results files of the code repository."
          % len(complete), "tab:S-cells", size=r"\scriptsize")
# The same index as an extended table, pointing into the supplement (xr prefix A:) and into the extended tables.
(OUT_REC / "rec_index.tex").write_text(HDR + "\n".join([
    "{\\scriptsize\\setlength{\\tabcolsep}{3pt}", "\\begin{longtable}{%s}" % CELLS_SPEC,
    "\\caption{Index: every simulation cell and where its results are reported (supplement tables carry the prefix "
    "S; Table~1 and Figure~1 are in the main text).}\\label{tab:R-index}\\\\",
    "\\toprule", CELLS_HDR + "\\\\", "\\midrule", "\\endfirsthead", "\\toprule", CELLS_HDR + "\\\\", "\\midrule",
    "\\endhead", "\\bottomrule", "\\endlastfoot"] + cell_rows("A:") + ["\\end{longtable}}", ""]), encoding="utf-8", newline="\n")

# Learner rate diagnostics and control calibration (the record keeps the original gate table)
QUESTION = {"C1": "normal approximation", "C2": "variance estimation", "C4a": "remainder $T_1$",
            "C4b": "remainder $T_2$", "C6": "TMLE versus AIPW", "sec. Psi": r"bridge to $\Psi$",
            "Power control": "control condition"}
body = []
for _, r in g5.iterrows():
    passed = r.status == "PASS"
    body.append(r"throttled learner R-d, $\alpha=%.2f$ & slope %s, 95\%% interval %s & %s\\" % (
        r.alpha, num(r.slope), interval(r.lo, r.hi),
        r"compatible with the designed slope $-\alpha$" if passed
        else r"faster than designed ($-\alpha$ outside the interval); not run"))
r = rbg.iloc[0]
body.append(r"random forest R-b & slope %s, one-sided 95\%% upper bound %s & below 0, consistent with $L^2$ "
            r"consistency; not below $-1/4$, so the quarter rate was not established\\" % (num(r.slope), num(r.upper_95)))
PCN = {"C1 shift": r"centre shifted by $0.30\,\sigma_{m,\mathrm{wor}}$ (mean $Z$)",
       "C2 wrong clustering: unit-level": r"units treated as clusters (mean $R_{\mathrm{var}}$)",
       "C2 wrong clustering: no m_h": r"no $m_h/(m_h-1)$ factor (mean $R_{\mathrm{var}}$)",
       "C2 wrong clustering: no stratification": r"strata ignored (mean $R_{\mathrm{var}}$)",
       "C4 aligned error": r"aligned nuisance error ($Q_{0.9}|T_2|$)"}
for _, r in g6.iterrows():
    if str(r.control).startswith("U10"):
        continue
    nm = [v for k, v in PCN.items() if str(r.control).startswith(k)][0]
    body.append(r"%s & census prediction %s; distance from the reference boundary %s (Monte Carlo SE %s) "
                r"& %s\\" % (nm, num(float(r.predicted)), num(float(r.distance)), num(float(r.mcse), 4),
                             "control, designed to leave the reference region" if bool(r.graded)
                             else "control, expected to stay inside the reference region; descriptive"))
body.append(r"aligned-error constant & 160 power-panel populations; $\kappa_{\mathrm{AE}}=%s$, giving "
            r"$Q_{0.9}(|T_2|/\sigma_{m,\mathrm{wor}})=%s$ & computed once on the census and reused\\" % (
                num(kappa["kappa"], 4), num(kappa["stat"], 2)))
longtable("tabS_rates.tex", "@{}>{\\raggedright\\arraybackslash}p{0.27\\linewidth}>{\\raggedright\\arraybackslash}"
          "p{0.4\\linewidth}>{\\raggedright\\arraybackslash}p{0.29\\linewidth}@{}",
          "learner or control & diagnostic & reading", body,
          "Learner rate diagnostics and control calibration, computed on the dry-run replicates before the main "
          "runs. For the throttled learner, the slope of the log out-of-fold error on $\\log m$ over "
          "$m=200,\\dots,3200$ (population-clustered CR2 standard errors, $t_{24}$ reference) is compared with "
          "the designed slope $-\\alpha$; for the forest, the slope over $m=50,\\dots,800$. These finite checks "
          "do not prove asymptotic rates. The control rows give each control's census-predicted value and its "
          "distance from the boundary of the relevant reference region.", "tab:S-rates", size=r"\scriptsize")

# All benchmarked estimates (the record keeps the labelled version, Table S-grades there)
body = []
for _, r in sci[sci.estimate.notna()].iterrows():
    d = 4 if ("paired" in r.statistic or "se2" in r.statistic) else 3
    reg = region(r.region_lo, r.region_hi)
    if str(r.disposition).startswith("BOUNDARY"):
        reg += "; boundary rate"
    if r.claim == "Power control":
        reg = (reg + "; expected inside, descriptive") if str(r.disposition).startswith("DIAGNOSTIC") \
            else ("outside " + reg)
    arm_txt = ([v for k, v in PCN.items() if str(r.arm).startswith(k)] or [arm_short(r.arm)])[0]
    body.append("%s & %s & %s & %s & %s & %s & %s\\\\" % (
        QUESTION.get(r.claim, tex_escape(r.claim)), SETTING.get(r.cell, tex_escape(r.cell)), arm_txt,
        SHORT.get(r.statistic, STAT.get(r.statistic, tex_escape(r.statistic))),
        num(r.estimate, d), interval(r.ci_lo, r.ci_hi, d), reg))
longtable("tabS_estimates.tex", "@{}llllrl>{\\raggedright\\arraybackslash}p{0.17\\linewidth}@{}",
          "question & setting & arm & statistic & estimate & 95\\% MC interval & benchmark", body,
          "Every estimate compared with its benchmark, with its 95\\% Monte Carlo interval (percentile interval "
          "of 10{,}000 resamples of whole populations; of populations for the census bridge). For a control "
          "condition the benchmark column gives the region it was designed to leave. Not listed, because no "
          "estimate exists: the throttled learner at $\\alpha=0.20$ and $0.35$, which was not run. Two variance "
          "controls were expected to stay inside the region and are descriptive. The forest's "
          "normal-approximation and TMLE values, for which the required convergence-rate condition was not "
          "established, are reported descriptively in Table~\\ref{tab:S-ms}. $Q_{0.9}|T_j|$ is the 0.9 quantile "
          "of $|T_j|/\\sigma_{m,\\mathrm{wor}}$; $\\Delta\\hat\\psi=(\\hat\\psi_{\\mathrm{TMLE}}-"
          "\\hat\\psi_{\\mathrm{AIPW}})/\\sigma_{m,\\mathrm{wor}}$; $\\Delta\\widehat{\\mathrm{se}}^2="
          "(\\widehat{\\mathrm{se}}^2_{\\mathrm{TMLE}}-\\widehat{\\mathrm{se}}^2_{\\mathrm{AIPW}})/\\sigma_m^2$; KS is "
          "the Kolmogorov distance.", "tab:S-estimates", size=r"\scriptsize", land=True)

# Trajectories and displays, with plain captions
body = []
for _, r in ms.sort_values(["m", "arm"]).iterrows():
    body.append("%d & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        r.m, ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.ks_Z), num(r.cover95_studentized),
        num(r.mean_R_var), num(r.sd_R_var), num(r.transfer_mean)))
longtable("tabS_ms.tex", "@{}rlrrrrrrr@{}",
          "$m$ & arm & mean $Z$ & SD$(Z)$ & KS & cov.\\ 95 & mean $R_{\\mathrm{var}}$ & SD$(R_{\\mathrm{var}})$ & transfer",
          body,
          "Many-strata regime ($m_h=2$, $J=20$, $f=0.10$), every arm and sample size; 400 populations $\\times$ 5 "
          "draws, populations equally weighted. $Z=(\\hat\\psi-\\Psi_N)/\\sigma_{m,\\mathrm{wor}}$, so mean $Z$ is "
          "the standardized bias; cov.\\ 95 is the coverage of the estimator's nominal 95\\% confidence interval "
          "$\\hat\\psi\\pm1.96\\,\\widehat{\\mathrm{se}}$ about $\\Psi_N$; $R_{\\mathrm{var}}=\\widehat{\\mathrm{se}}^2/"
          "\\sigma_m^2$. Table~\\ref{tab:sim} compares the largest sizes ($m=1600$ for R-a; $m=800$ for the forest "
          "and TMLE) with the benchmarks. For the forest R-b the required convergence-rate condition was not "
          "established, so its values are descriptive.", "tab:S-ms", size=r"\scriptsize")
body = []
for _, r in fs.sort_values(["m", "arm"]).iterrows():
    tmle = isinstance(r.remainder_note, str)
    body.append("%d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        r.m, ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.cover95_studentized), num(r.mean_R_var),
        "n.r." if tmle else ("---" if r.arm == "oracle" else "%d" % r.t2_unresolved),
        num(r.q90_abs_T1_std), num(r.q90_abs_T2_std), num(r.q90_sqrt_m_T1_over_rhat, 2),
        num(r.q90_T2_over_rQ_rg, 2), num(r.mean_r_hat)))
longtable("tabS_fs.tex", "@{}rlrrrrrrrrrr@{}",
          "$m$ & arm & mean $Z$ & SD$(Z)$ & cov.\\ 95 & $R_{\\mathrm{var}}$ & unres. & $Q_{0.9}|T_1|$ & "
          "$Q_{0.9}|T_2|$ & env.\\ $T_1$ & env.\\ $T_2$ & $\\hat r_m$", body,
          "Fixed-strata regime ($H=5$, $f=0.10$), every arm and sample size; 400 populations $\\times$ 5 draws. "
          "$Q_{0.9}|T_j|$ is the 0.9 quantile of $|T_j|/\\sigma_{m,\\mathrm{wor}}$ over replicates with a resolved "
          "remainder; ``unres.'' counts replicates whose $T_2$ (and hence $T_1$) could not be resolved "
          "numerically and are kept as missing. The envelope ratios $\\sqrt m|T_1|/\\hat r_m$ and "
          "$|T_2|/(\\hat r_Q\\hat r_g)$ (0.9 quantiles) and the mean out-of-fold error $\\hat r_m$ are descriptive. "
          "R-a TMLE remainders are not reported (n.r.), because their numerical evaluation was not validated "
          "(Section~\\ref{supp:numres}). The throttled learner R-d runs at $m\\ge200$ only; $\\alpha=0.20$ and $0.35$ "
          "were not run.", "tab:S-fs", size=r"\scriptsize", land=True)
body = []
for _, r in e2f.iterrows():
    body.append("%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        SETTING[r.cell], ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.ks_Z), num(r.mean_R_var),
        num(r.cov_se_80), num(r.cov_se_90), num(r.cov_se_95)))
longtable("tabS_e2fit.tex", "@{}llrrrrrrr@{}",
          "setting & arm & mean $Z$ & SD$(Z)$ & KS & mean $R_{\\mathrm{var}}$ & cov.\\ 80 & cov.\\ 90 & cov.\\ 95",
          body,
          "Fitted cells at first-stage fraction $f=0.50$ (160 populations $\\times$ 5 draws). Coverage is that of "
          "the estimator's nominal confidence interval $\\hat\\psi\\pm z\\,\\widehat{\\mathrm{se}}$ about $\\Psi_N$; "
          "it is reported beside its components, the standardized bias (mean $Z$), the dispersion relative to the "
          "without-replacement norming (SD$(Z)$) and the variance ratio against the ultimate-cluster target "
          "($R_{\\mathrm{var}}$). The forest R-b rows are descriptive.", "tab:S-e2fit", size=r"\scriptsize")
body = []
for _, r in fpc.sort_values(["J", "population"], ascending=[False, True]).iterrows():
    body.append("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s\\\\" % (
        JF[int(r.J)], int(r.population), num(r.sigma_ratio), num(r.cov_wor_95), num(r.cov_det_sigma_m_80),
        num(r.cov_det_sigma_m_90), num(r.cov_det_sigma_m_95), num(r.overlay_95), num(r.cov_se_95),
        num(r.c3a_statistic)))
longtable("tabS_fpc.tex", "@{}rrrrrrrrrr@{}",
          "$f$ & pop. & $\\sigma_m/\\sigma_{m,\\mathrm{wor}}$ & wor 95 & $\\sigma_m$ 80 & $\\sigma_m$ 90 & "
          "$\\sigma_m$ 95 & overlay 95 & $\\widehat{\\mathrm{se}}$ 95 & gap error", body,
          "Oracle interval hierarchy across the first-stage fraction sweep ($H=50$, $m_h=2$, $m=100$), each of 12 "
          "populations on 2{,}000 redraws (10{,}000 at $f=0.50$ for the gap error). Columns give coverage of "
          "$\\hat\\psi_{\\mathrm{oracle}}\\pm z\\sigma_{m,\\mathrm{wor}}$ (``wor''), of "
          "$\\hat\\psi_{\\mathrm{oracle}}\\pm z\\sigma_m$ at 0.80, 0.90, 0.95, the analytic overlay "
          "$2\\Phi(z\\sigma_m/\\sigma_{m,\\mathrm{wor}})-1$ at 0.95, the coverage of "
          "$\\hat\\psi_{\\mathrm{oracle}}\\pm z\\widehat{\\mathrm{se}}$, and the gap error "
          "$\\widehat{\\mathrm{Var}}(S_{\\mathrm{des}})/\\sigma_m^2-(1-\\Delta_m/\\sigma_m^2)$. The binomial standard "
          "error of one population's coverage near 0.95 is about %s." % M["SimCovSEBig"],
          "tab:S-fpc", size=r"\scriptsize")
body = []
for _, r in strict.iterrows():
    body.append("strict, $f=0.50$ & $H=%d$ & %s & %s & %s & %s\\\\" % (
        r.H, num(r.min_Delta_over_sigma2_m), num(r.median_Delta_over_sigma2_m), num(r.max_Delta_over_sigma2_m),
        num(r.median_sigma_ratio)))
longtable("tabS_strict.tex", "@{}llrrrr@{}",
          "trajectory & $H$ & min $\\Delta_m/\\sigma_m^2$ & median & max & median $\\sigma_m/\\sigma_{m,\\mathrm{wor}}$",
          body,
          "Relative without-replacement gap along the fixed-fraction trajectory $m_h=2$, $J=4$ (census algebra, "
          "120 populations per size, no sampling). It stays away from zero, which is compatible with normings that "
          "stay strictly separated; a finite trajectory cannot establish an asymptotic limit.",
          "tab:S-strict", size=r"\scriptsize")
body = []
for _, r in vanish.sort_values("f", ascending=False).iterrows():
    body.append("%d & %d & %s & %s & %s & %s & %s\\\\" % (
        r.H, r.J, num(r.f, 4), num(r.median_bridge), num(r.q90_bridge), num(r.max_bridge),
        num(r.ratio_mean_Delta_over_var_Smod, 2)))
longtable("tabS_vanish.tex", "@{}rrrrrrr@{}",
          "$H$ & $J$ & $f$ & median & $Q_{0.9}$ & max & $\\bar\\Delta_m/\\widehat{\\mathrm{Var}}(\\Psi_N-\\Psi)$",
          body,
          "Vanishing-fraction trajectory of the bridge $|\\Psi_N-\\Psi|/\\sigma_{m,\\mathrm{wor}}$ (census, 120 "
          "populations per fraction, $m_h=2$). The smallest fraction, $f=0.0025$, was included after a "
          "leading-order calculation predicted a 0.9 quantile of about 0.07 there. The quantiles displayed in "
          "this table are observed census summaries. This comparison illustrates the prediction rather than "
          "testing it independently. The last column reports "
          "the leading-order relation between the mean without-replacement gap and the across-population variance "
          "of $\\Psi_N-\\Psi$, which is outside Theorem~\\ref{thm:clt}.", "tab:S-vanish", size=r"\scriptsize")
body = []
for _, r in rc.iterrows():
    body.append("%s & %s & %s & %s & %s & %s\\\\" % (
        "R-c comparison, $m=400$", ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.cover95_studentized), num(r.mean_R_var)))
for _, r in vs.iterrows():
    body.append("%s & %s & %s & %s & %s & %s\\\\" % (
        "$V=%s$, $m=200$" % r.cell[-1], ARM[r.arm], num(r.mean_Z), num(r.sd_Z), num(r.cover95_studentized),
        num(r.mean_R_var)))
longtable("tabS_e3.tex", "@{}llrrrr@{}",
          "cell & arm & mean $Z$ & SD$(Z)$ & cov.\\ 95 & mean $R_{\\mathrm{var}}$", body,
          "Nuisance complexity and folds (fixed strata, 400 populations $\\times$ 5 draws): the interpolating "
          "forest R-c under PSU-level folds, no cross-fitting, unit-level folds and unit-level folds with "
          "$C_{\\mathrm{psu}}$ coarsened, beside the R-a fits on the same samples; and the number of folds ($V=2$ "
          "versus $V=5$, random forest R-b). The required convergence-rate condition was not established for R-c "
          "or R-b, so these rows are descriptive.", "tab:S-e3", size=r"\scriptsize", land=True)
body = []
for (cell, arm), v in unresolved.items():
    body.append("%s & %s & %d & %s\\\\" % (tex_escape(cell), ARM[arm], v,
                                          "not reported (not validated)" if arm.endswith("TMLE") else "kept as missing"))
longtable("tabS_unresolved.tex", "@{}llrl@{}", "cell & arm & unresolved & treatment", body,
          "Arm-by-replicate entries whose second-order remainder $T_2$ could not be resolved by the refined "
          "quadrature (%s entries in total, all at $m\\le400$). A sample can contribute an entry for more than one "
          "arm, so this is not a count of distinct samples. $T_2$ and $T_1$ are recorded as missing and every "
          "estimate is kept; R-a TMLE remainders are not reported in any case. No sample size compared with a "
          "benchmark is affected." % sum(unresolved.values()),
          "tab:S-unresolved", size=r"\scriptsize")
# Remove supplement copies of tables that now live only in the record (written there above).
for _f in ("tabS_grades.tex", "tabS_dispositions.tex", "tabS_popindexed.tex", "tabS_e1oracle.tex",
           "tabS_unequal.tex", "tabS_e4.tex", "tabS_gates.tex", "tabS_g7.tex"):
    (OUT_SUPP / _f).unlink(missing_ok=True)

# ----------------------------------------------------------------------------------------------
# Application (SPA census-resampling study): aggregate outputs only, see app_displays.py
# ----------------------------------------------------------------------------------------------
sys.dont_write_bytecode = True        # no __pycache__ beside the script
import app_displays  # noqa: E402  (same directory; run as a script from anywhere)
app_displays.build(ROOT, OUT_MAIN, OUT_SUPP, mac, num, interval, longtable,
                   "% GENERATED by src/displays/make_sim_displays.py (app_displays.py) from saved SPA aggregate "
                   "outputs; do not edit by hand.\n")

# ----------------------------------------------------------------------------------------------
# Macros file
# ----------------------------------------------------------------------------------------------
lines = ["% GENERATED by src/displays/make_sim_displays.py from the grading outputs; do not edit by hand.",
         "% Grading pass grade-20261017T175903Z (package 3603568) over production package 19b4ae9."]
for k in sorted(M):
    lines.append("\\newcommand{\\%s}{%s}" % (k, M[k]))
(OUT_MAIN / "sim_numbers.tex").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
print("wrote %d macros, Figure 1, Table 1 (2 panels) and %d supplement tables" % (
    len(M), len(list(OUT_SUPP.glob("tabS_*.tex")))))
