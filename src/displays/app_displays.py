"""Application displays (SPA census-resampling study), built from saved aggregate outputs only.

Called by make_sim_displays.py, which supplies its formatting helpers. Reads only aggregate CSV outputs under
results/application/ (no record-level file is opened, nothing is estimated) and writes:
    results/displays/main/tab2_application.tex       (main-text Table 2 rows)
    results/displays/supplement/tabS_app_*.tex        (supplement S4 tables)
and \\App... macros added to the shared number file.

Method-to-reference mapping (resampling script src/application/R/02-resample.R and the protocol):
    Fully-Aware -> tmle_sf, Partially-Aware -> tmle_sf, Non-Aware -> tmle_sf_unw, Fully-Aware-CV -> tmle_cv,
    Fully-Aware-CF -> tmle_cf, AIPW-CF -> aipw_cf, the two [unw+design C] arms -> tmle_cf_unwd / aipw_cf_unwd.
The mapping is asserted below against the `reference` column that the calibration files themselves carry.
"""
import re
from pathlib import Path

import pandas as pd

SURVEYS = [("MW_SPA_2013-14", "Malawi 2013--14", "Mw"), ("HT_SPA_2017-18", "Haiti 2017--18", "Hts"),
           ("HT_SPA_2013", "Haiti 2013", "Htt")]
FAMILY = {"Fully-Aware": "tmle_sf", "Partially-Aware": "tmle_sf", "Non-Aware": "tmle_sf_unw",
          "Fully-Aware-CV": "tmle_cv", "Fully-Aware-CF": "tmle_cf", "AIPW-CF": "aipw_cf",
          "Fully-Aware-CF [unw+design C]": "tmle_cf_unwd", "AIPW-CF [unw+design C]": "aipw_cf_unwd"}
LABEL = {"Fully-Aware": "TMLE, single fit, design-aware",
         "Partially-Aware": "TMLE, single fit, weights only",
         "Non-Aware": "TMLE, single fit, no design",
         "Fully-Aware-CV": "TMLE, internal CV, design-aware",
         "Fully-Aware-CF": r"\textbf{cross-fitted TMLE}",
         "AIPW-CF": r"\textbf{cross-fitted AIPW}",
         "Fully-Aware-CF [unw+design C]": "cross-fitted TMLE, sensitivity",
         "AIPW-CF [unw+design C]": "cross-fitted AIPW, sensitivity"}
F_NOMINAL = 0.25


def build(root, out_main, out_supp, mac, num, interval, longtable, hdr):
    res = root / "results" / "application"
    frames = pd.read_csv(res / "frames.csv").set_index("survey")      # scripts/application/extract_frames.py
    cal, ref, spl, ovl, diag, frame = {}, {}, {}, {}, {}, {}
    for code, name, key in SURVEYS:
        cal[code] = pd.read_csv(res / ("resample_%s_f250_calibration.csv" % code))
        ref[code] = pd.read_csv(res / ("resample_%s_reference.csv" % code)).set_index("family").pseudo_reference
        spl[code] = pd.read_csv(res / ("resample_%s_reference_splits.csv" % code))
        ovl[code] = pd.read_csv(res / ("resample_%s_census_overlap.csv" % code))
        diag[code] = pd.read_csv(res / ("resample_%s_f250_diagnostics_summary.csv" % code)).set_index("stat")
        audit = [p for p in res.glob("count-audit_%s_*.csv" % code)]
        assert len(audit) == 1, code
        a = pd.read_csv(audit[0]).query("section == 'metric'").set_index("name").value.astype(float)
        fr_row = frames.loc[code]
        frame[code] = dict(records=int(fr_row.records), facilities=int(fr_row.facilities), excluded=int(fr_row.excluded),
                           strata=int(a["n_strata"]), min_stratum=int(a["min_facilities_per_stratum"]))
        assert frame[code]["facilities"] == int(a["n_facilities"]), code
        # every method's reference column must equal its family's pseudo-reference
        for _, r in cal[code].iterrows():
            assert abs(r.reference - ref[code][FAMILY[r.method]]) < 1e-12, (code, r.method)
            assert r.R_ok == r.R_attempted == 500, (code, r.method)

    # ---------------- macros
    for code, name, key in SURVEYS:
        fr = frame[code]
        mac("App%sFacilities" % key, "{:,}".format(fr["facilities"]).replace(",", "{,}"))
        mac("App%sRecords" % key, "{:,}".format(fr["records"]).replace(",", "{,}"))
        mac("App%sExcluded" % key, str(fr["excluded"]))
        mac("App%sPsu" % key, str(int(diag[code].loc["n_psu", "median"])))
        mac("App%sRecMin" % key, str(int(diag[code].loc["n_rec", "min"])))
        mac("App%sRecMax" % key, str(int(diag[code].loc["n_rec", "max"])))
        mac("App%sMinStratum" % key, str(fr["min_stratum"]))
        mac("App%sMinSampled" % key, str(int(round(F_NOMINAL * fr["min_stratum"]))))
        for meth, mk in (("Fully-Aware-CF", "Tmle"), ("AIPW-CF", "Aipw")):
            r = cal[code].query("method == @meth").iloc[0]
            mac("App%s%sVr" % (key, mk), num(r.var_ratio, 3))
            mac("App%s%sVrMcse" % (key, mk), num(r.var_ratio_mcse, 3))
            mac("App%s%sSeSd" % (key, mk), num(r.mean_se_over_sd, 3))
            mac("App%s%sSePct" % (key, mk), "%.0f" % (100 * (r.mean_se_over_sd - 1)))
            mac("App%s%sCov" % (key, mk), "%.1f" % (100 * r.coverage))
            mac("App%s%sRef" % (key, mk), num(100 * r.reference, 2))
        sens = cal[code].query("method == 'Fully-Aware-CF [unw+design C]'").iloc[0]
        mac("App%sSensTmleVr" % key, num(sens.var_ratio, 3))
        mac("App%sSensTmleCov" % key, "%.1f" % (100 * sens.coverage))
        o = ovl[code].set_index(["covariate_set", "stat"]).value
        mac("App%sOverlapSens" % key, "%.1f" % (100 * o[("C + design (sensitivity)", "share_truncated")]))
        mac("App%sOverlapPrim" % key, "%.1f" % (100 * o[("primary C", "share_truncated")]))
    prim = [cal[c].query("method == @m").iloc[0] for c, _, _ in SURVEYS for m in ("Fully-Aware-CF", "AIPW-CF")]
    mac("AppCovMin", "%.1f" % (100 * min(r.coverage for r in prim)))
    mac("AppCovMax", "%.1f" % (100 * max(r.coverage for r in prim)))
    mac("AppSePctMin", "%.0f" % (100 * (min(r.mean_se_over_sd for r in prim) - 1)))
    mac("AppSePctMax", "%.0f" % (100 * (max(r.mean_se_over_sd for r in prim) - 1)))
    bound = 1 / (1 - F_NOMINAL)
    mw = cal["MW_SPA_2013-14"].query("method == 'Fully-Aware-CF'").iloc[0]
    zb = (mw.var_ratio - bound) / mw.var_ratio_mcse
    assert zb > 2
    mac("AppMwTmleBoundZ", "%.1f" % zb)
    others = [r for r in prim if not (r is prim[0])]
    assert all(r.var_ratio < bound for r in others)
    hta = cal["HT_SPA_2013"].query("method == 'AIPW-CF'").iloc[0]
    assert (hta.var_ratio - 1) / hta.var_ratio_mcse < 2
    mac("AppMinSampledAll", str(min(int(round(F_NOMINAL * frame[c]["min_stratum"])) for c, _, _ in SURVEYS)))
    mac("AppFpcBound", num(1 / (1 - F_NOMINAL), 2))
    mac("AppFpcBoundSe", num((1 / (1 - F_NOMINAL)) ** 0.5, 2))
    e1 = all((cal[c].query("method == @m").iloc[0].var_ratio
              >= 1 - 2 * cal[c].query("method == @m").iloc[0].var_ratio_mcse)
             for c, _, _ in SURVEYS for m in ("Fully-Aware-CF", "AIPW-CF"))
    assert e1
    covmc = [cal[c].query("method == @m").iloc[0].coverage_mcse for c, _, _ in SURVEYS for m in ("Fully-Aware-CF", "AIPW-CF")]
    mac("AppCovMcseMin", "%.1f" % (100 * min(covmc)))
    mac("AppCovMcseMax", "%.1f" % (100 * max(covmc)))
    # E2: partially-aware below fully-aware in every survey
    e2 = all(cal[c].query("method == 'Partially-Aware'").iloc[0].var_ratio
             < cal[c].query("method == 'Fully-Aware'").iloc[0].var_ratio for c, _, _ in SURVEYS)
    assert e2

    # ---------------- main Table 2
    rows = []
    for code, name, key in SURVEYS:
        first = True
        for meth, lab in (("Fully-Aware-CF", "cross-fitted TMLE"), ("AIPW-CF", "cross-fitted AIPW")):
            r = cal[code].query("method == @meth").iloc[0]
            rows.append("%s & %s & %s (%s) & %s & %s (%s) & %s\\\\" % (
                name if first else "", lab, num(r.var_ratio, 3), num(r.var_ratio_mcse, 3), num(r.mean_se_over_sd, 3),
                "%.1f" % (100 * r.coverage), "%.1f" % (100 * r.coverage_mcse), num(100 * r.reference, 2)))
            first = False
    (out_main / "tab2_application.tex").write_text(hdr + "\n".join(rows) + "\n", encoding="utf-8", newline="\n")

    # ---------------- S4 tables
    body = []
    for code, name, key in SURVEYS:
        fr, d = frame[code], diag[code]
        body.append("%s & %d & {%s} & %d & %d & %d & %d--%d & %s--%s & %s\\\\" % (
            name, fr["strata"], "{:,}".format(fr["facilities"]).replace(",", "{,}"),
            fr["records"], fr["excluded"], int(d.loc["n_psu", "median"]), int(d.loc["n_rec", "min"]),
            int(d.loc["n_rec", "max"]), num(d.loc["w_min", "min"], 1), num(d.loc["w_max", "max"], 1),
            num(d.loc["kish_deff_w", "median"], 2)))
    longtable("tabS_app_frames.tex", "@{}lrrrrrrrr@{}",
              "survey & strata & facilities & records & excl. & sampled fac. & sampled rec. & weights & Kish deff",
              body,
              "SPA analysis frames and the imposed design. Frames are the eligible facility and client records "
              "(completed facility, observation and exit interview), after excluding records with missing exposure or "
              "outcome (excl.). Each resample draws a quarter of the facilities in every facility-type stratum by "
              "SRS-WOR, then up to three records per sampled facility by SRS-WOR; the last four columns summarize the "
              "500 resamples (range of sampled records, range of design weights, median Kish design effect of the "
              "weights).", "tab:S-app-frames", size=r"\scriptsize")
    body = []
    for code, name, key in SURVEYS:
        first = True
        for _, r in cal[code].iterrows():
            body.append("%s & %s & %s (%s) & %s & %s (%s) & %s (%s) & %s\\\\" % (
                name if first else "", LABEL[r.method], num(r.var_ratio, 3), num(r.var_ratio_mcse, 3),
                num(r.mean_se_over_sd, 3), num(100 * r.bias, 2), num(100 * r.bias_mcse, 2),
                "%.1f" % (100 * r.coverage), "%.1f" % (100 * r.coverage_mcse), num(100 * r.mean_ci_width, 1)))
            first = False
    longtable("tabS_app_arms.tex", "@{}llrrrrr@{}",
              "survey & method & variance ratio (MCSE) & SE/SD & bias, pp (MCSE) & coverage, \\% (MCSE) & CI width, pp",
              body,
              "All eight methods, each on 500 of 500 successful resamples per survey. The "
              "variance ratio is the mean estimated variance divided by the empirical variance of the estimates "
              "across resamples, with a replicate-bootstrap Monte Carlo standard error; SE/SD is the mean estimated "
              "standard error divided by the empirical standard deviation (not the square root of the variance "
              "ratio). Bias and coverage of the nominal 95\\% confidence interval are relative to each method "
              "family's whole-frame reference value (Table~\\ref{tab:S-app-ref}), not to a known truth; their Monte Carlo "
              "standard errors treat that reference as fixed. Single-fit arms fit the nuisances once on the whole "
              "sample; ``weights only'' ignores strata and clusters in the standard error; ``no design'' ignores the "
              "design entirely. The two sensitivity arms use unweighted nuisances and add facility type and log "
              "caseload to the covariates, so they target a differently adjusted contrast.",
              # land=False: supp_app.tex opens the landscape page itself, so the S4.3 heading and introduction share
              # the page with this table instead of standing alone on a portrait page (layout fix L2).
              "tab:S-app-arms", size=r"\scriptsize", land=False)
    body = []
    FAMLAB = {"tmle_sf": "TMLE, single fit (weighted nuisances)", "tmle_sf_unw": "TMLE, single fit (unweighted)",
              "tmle_cv": "TMLE, internal CV", "tmle_cf": "cross-fitted TMLE", "aipw_cf": "cross-fitted AIPW",
              "tmle_cf_unwd": "cross-fitted TMLE, sensitivity", "aipw_cf_unwd": "cross-fitted AIPW, sensitivity"}
    for code, name, key in SURVEYS:
        g = spl[code].groupby("family").value.agg(["count", "min", "max"])
        first = True
        for fam in ["tmle_sf", "tmle_sf_unw", "tmle_cv", "tmle_cf", "aipw_cf", "tmle_cf_unwd", "aipw_cf_unwd"]:
            body.append("%s & %s & %s & %s to %s\\\\" % (name if first else "", FAMLAB[fam], num(100 * ref[code][fam], 2),
                                                         num(100 * g.loc[fam, "min"], 2), num(100 * g.loc[fam, "max"], 2)))
            assert int(g.loc[fam, "count"]) == 20
            first = False
    longtable("tabS_app_ref.tex", "@{}llrl@{}", "survey & estimator family & reference value, pp & range over 20 splits",
              body,
              "Whole-frame reference values: each estimator family run on the whole analysis frame with unit weights, in its "
              "own nuisance mode and covariate set, over 20 fixed-seed fold splits; the median is the reference. It is "
              "an algorithmic reference value for that family, not a known causal effect and not a known finite-population "
              "value of the theorem's estimand. Values are adjusted urban-minus-rural differences in the proportion "
              "``very satisfied'', in percentage points.", "tab:S-app-ref", size=r"\scriptsize")
    body = []
    for code, name, key in SURVEYS:
        o = ovl[code].set_index(["covariate_set", "stat"]).value
        for cs, lab in (("primary C", "primary"), ("C + design (sensitivity)", "sensitivity")):
            body.append("%s & %s & %s & %s & %s & %s\\\\" % (
                name if cs == "primary C" else "", lab, num(o[(cs, "p01")], 3), num(o[(cs, "p50")], 3),
                num(o[(cs, "p99")], 3), "%.1f" % (100 * o[(cs, "share_truncated")])))
    longtable("tabS_app_overlap.tex", "@{}llrrrr@{}",
              "survey & covariates & $P_{0.01}$ & median & $P_{0.99}$ & outside $[0.05,0.95]$, \\%", body,
              "Census overlap diagnostic: out-of-fold SuperLearner estimates of $P(A=1\\mid C)$ on the whole frame, for "
              "the primary and sensitivity covariate sets. A covariate set is flagged as overlap-strained when "
              "more than 10\\% of values fall outside $[0.05,0.95]$; this is a model-based diagnostic, not a proof of "
              "positivity, and it was never used to respecify an analysis.", "tab:S-app-overlap", size=r"\scriptsize")
