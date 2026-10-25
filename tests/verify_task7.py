"""Verify the task-7 recovery: the recovered folder against the design and provenance, its engine identity against task 5
of the same cell, and the recovered replicates against the values the first, timed-out run printed.
    python tests/verify_task7.py <assembled production folder>/P-e3-ctrl-m400 results/simulation/task7_timed_out_printed.csv
(The original check also compared the cell's other 15 task folders byte for byte between the two cluster archives; the
release ships those folders once, in the production archive.)"""
import collections, csv, glob, json, os, re, sys

new_cell, printed_csv = sys.argv[1], sys.argv[2]
old_cell = new_cell
csv.field_size_limit(10**9)
SPEC = "58fe3c9f53dda3c54c0e20c06216ff00e1411541933a23be8a9bbf8f64eead61"

t7 = [d for d in glob.glob(os.path.join(new_cell, "*t007-20261017*"))][0]
m = json.load(open(os.path.join(t7, "manifest.json"), encoding="utf-8"))
print("folder", os.path.basename(t7))
print("mode", m.get("mode"), "| spec ok", (m.get("spec_freeze") or {}).get("sha256") == SPEC,
      "| package", (m.get("git") or {}).get("package_commit", "")[:7], "| exit", (m.get("checks") or {}).get("exit_status"),
      "| n_required", (m.get("checks") or {}).get("n_required"), "| wall s", round(m.get("wall_seconds", 0)))
res = list(csv.DictReader(open(os.path.join(t7, "results.csv"), encoding="utf-8")))
arms = collections.Counter(r["arm"] for r in res)
pops = sorted({int(r["pop_id"]) for r in res}); reps = sorted({int(r["rep_id"]) for r in res})
print("rows", len(res), "| arms", dict(arms), "| pops", pops[0], "-", pops[-1], f"({len(pops)})", "| reps", reps)
print("FAILED", sum(r.get("arm_status") == "FAILED" for r in res), "| unresolved T2", sum(r.get("T2_status") == "unresolved" for r in res))
ck = list(csv.DictReader(open(os.path.join(t7, "checks.csv"), encoding="utf-8")))
req = [r for r in ck if r["required"] == "TRUE"]
print("required checks", len(req), "| failed", sum(r["pass"] != "TRUE" for r in req))

# engine identity with a task of the first run (same package => same engine-file hashes)
ref = json.load(open(glob.glob(os.path.join(old_cell, "*t005-*", "manifest.json"))[0], encoding="utf-8"))
same = ref.get("engine_files_sha256") == m.get("engine_files_sha256") and ref.get("renv_lock") == m.get("renv_lock")
print("engine files + renv.lock identical to production task 5:", same)

# keyed draws: the timed-out run printed oracle / AIPW / AE estimates per replicate
printed = {}
for r in csv.DictReader(open(printed_csv, encoding="utf-8")):
    printed[(int(r["pop_id"]), int(r["rep_id"]))] = dict({k: r[k] for k in ("oracle", "aipw", "aipw_se", "t2")},
                                                         **({"ae": r["ae_aipw"]} if r["ae_aipw"] else {}))
by = {(int(r["pop_id"]), int(r["rep_id"]), r["arm"]): r for r in res}
def f5(x): return f"{float(x):.5f}"
match = mism = 0
for (p, rp), v in sorted(printed.items()):
    o = by.get((p, rp, "oracle")); a = by.get((p, rp, "Ra-AIPW")); e = by.get((p, rp, "AE-AIPW"))
    checks = [(f5(o["psi_hat"]), v["oracle"]), (f5(a["psi_hat"]), v["aipw"]), (f5(a["se_hat"]), v["aipw_se"]),
              (f"{float(a['T2']):+.2e}", v["t2"])]
    if "ae" in v and e is not None: checks.append((f5(e["psi_hat"]), v["ae"]))
    for got, want in checks:
        if got == want: match += 1
        else: mism += 1; print("  MISMATCH", p, rp, got, want)
print(f"timed-out run printed {len(printed)} replicates; recovered values agree to the printed digits: {match} matched, {mism} mismatched")
