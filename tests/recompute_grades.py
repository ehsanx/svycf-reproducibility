"""Independent recomputation of the key graded point estimates from production results.csv, compared with
grades.csv. Plain statistics (every population has 5 replicates, so the grader's population weights are all 1):
type-7 quantile, mean, SD (n - 1), Kolmogorov distance. Plus an approximate population-block bootstrap interval
(numpy RNG, not the grader's keyed stream) for the RED-FLAG remainder rows."""
import csv, glob, math, os, sys
import numpy as np

prod, recov_cell, grades_csv = sys.argv[1], sys.argv[2], sys.argv[3]
csv.field_size_limit(10**9)

def load(cell, arm):
    rows = []
    folders = sorted(glob.glob(os.path.join(prod, cell, "*")))
    if cell == "P-e3-ctrl-m400":   # the task-7 recovery replaces the empty timed-out folder
        folders = [f for f in folders if "t007" not in os.path.basename(f)] + glob.glob(os.path.join(recov_cell, "*t007-20261017*"))
    for f in folders:
        p = os.path.join(f, "results.csv")
        if not os.path.isfile(p): continue
        for r in csv.DictReader(open(p, encoding="utf-8")):
            if r["arm"] == arm and r.get("arm_status", "OK") == "OK": rows.append(r)
    return rows

def col(rows, c): return np.array([float(r[c]) for r in rows if r[c] not in ("", "NA")])
def q7(x, p): return float(np.quantile(x, p))          # numpy's default "linear" is type 7
def ks(z):
    z = np.sort(z); n = len(z); cw = np.arange(1, n + 1) / n
    from math import erf, sqrt
    p = np.array([0.5 * (1 + erf(v / sqrt(2))) for v in z])
    return float(max(np.max(np.abs(cw - p)), np.max(np.abs(np.concatenate([[0], cw[:-1]]) - p))))

G = {}
for r in csv.DictReader(open(grades_csv, encoding="utf-8")):
    G[(r["claim"], r["cell"], r["arm"], r["statistic"])] = r

out = []
def rep(claim, cell, arm, stat, val, boot=None):
    g = G.get((claim, cell, arm, stat), {})
    est = g.get("estimate", "")
    d = abs(val - float(est)) if est not in ("", "NA") else float("nan")
    line = f"{claim:6s} {cell:14s} {arm:44s} {stat:24s} mine {val:.6g}  grader {est[:10]:10s} |diff| {d:.1e}"
    if boot is not None: line += f"  approx-boot95 [{boot[0]:.4f}, {boot[1]:.4f}] grader [{g.get('ci_lo','')[:7]}, {g.get('ci_hi','')[:7]}] -> {g.get('disposition','')[:12]}"
    print(line)

def pboot(rows, c, fun, B=2000, seed=1):
    pops = sorted({r["pop_id"] for r in rows}); by = {p: [] for p in pops}
    for r in rows:
        if r[c] not in ("", "NA"): by[r["pop_id"]].append(float(r[c]))
    arr = [np.array(by[p]) for p in pops]; rng = np.random.default_rng(seed); vals = []
    for _ in range(B):
        idx = rng.integers(0, len(arr), len(arr)); vals.append(fun(np.concatenate([arr[i] for i in idx])))
    return np.quantile(vals, [0.025, 0.975])

# C4a / C4b remainders at the terminal rungs
ra = load("P-fs-m1600", "Ra-AIPW")
x = np.abs(col(ra, "T1_over_sigma_wor"))
rep("C4a", "P-fs-m1600", "Ra-AIPW", "Q90_abs_T1_std", q7(x, .9), pboot(ra, "T1_over_sigma_wor", lambda v: q7(np.abs(v), .9)))
for a in ("Rd-AIPW-a025", "Rd-AIPW-a030"):
    rd = load("P-fs-m3200", a)
    rep("C4a", "P-fs-m3200", a, "Q90_abs_T1_std", q7(np.abs(col(rd, "T1_over_sigma_wor")), .9), pboot(rd, "T1_over_sigma_wor", lambda v: q7(np.abs(v), .9)))
    rep("C4b", "P-fs-m3200", a, "Q90_abs_T2_std", q7(np.abs(col(rd, "T2_over_sigma_wor")), .9), pboot(rd, "T2_over_sigma_wor", lambda v: q7(np.abs(v), .9)))
    print(f"       rows {len(rd)}; unresolved {sum(r.get('T2_status') == 'unresolved' for r in rd)}")
# C1 / C2
for cell, arm in (("P-ms-H800", "Ra-AIPW"), ("P-fpc-fit-J4", "Ra-AIPW"), ("P-strict-m200", "Ra-AIPW")):
    rr = load(cell, arm); z = col(rr, "Z_wor")
    rep("C1", cell, arm, "mean_Z", float(z.mean())); rep("C1", cell, arm, "sd_Z_minus_1", float(z.std(ddof=1) - 1)); rep("C1", cell, arm, "KS_Z", ks(z))
    if cell == "P-fpc-fit-J4": print(f"       {cell} {arm}: n = {len(z)}, populations {len({r['pop_id'] for r in rr})}")
rr = load("P-ms-H800", "Ra-AIPW"); orc = {(r["pop_id"], r["rep_id"]): float(r["R_var"]) for r in load("P-ms-H800", "oracle")}
rv = col(rr, "R_var")
rep("C2", "P-ms-H800", "Ra-AIPW", "mean_R_var_minus_1", float(rv.mean() - 1)); rep("C2", "P-ms-H800", "Ra-AIPW", "sd_R_var", float(rv.std(ddof=1)))
rep("C2", "P-ms-H800", "Ra-AIPW", "transfer", float(np.mean([abs(float(r["R_var"]) - orc[(r["pop_id"], r["rep_id"])]) for r in rr])))
# C6 (Ra-TMLE at m = 800, many-strata)
tm = load("P-ms-H400", "Ra-TMLE")
rep("C6", "P-ms-H400", "Ra-TMLE", "Q90_abs_paired_psi_std", q7(np.abs(np.array([float(r["paired_tmle_minus_aipw"]) / float(r["sigma_m_wor"]) for r in tm])), .9))
rep("C6", "P-ms-H400", "Ra-TMLE", "Q90_abs_se2_diff_std", q7(np.abs(col(tm, "se2_diff_over_sigma2_m")), .9))
# the aligned-error power control
ae = load("P-e3-ctrl-m400", "AE-AIPW")
rep("Power control", "P-e3-ctrl-m400", "C4 aligned error: Q_0.9(|T2|/sigma_m,wor), AE-AIPW", "Q90_abs_T2_std", q7(np.abs(col(ae, "T2_over_sigma_wor")), .9))
print(f"       AE rows {len(ae)} (expected 800)")
