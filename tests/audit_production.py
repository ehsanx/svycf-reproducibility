"""Audit of the production tree: completeness against the cell table, required checks, exit status, FAILED fits and
unresolved T2. Reads saved files only.
    python tests/audit_production.py <production folder> config/simulation/production-cells.tsv [package] [--csv FILE]
package: the package identity every run must record (default 19b4ae9, the recorded production; a reproduction
records release:<id>). --csv writes the unresolved-T2 counts by cell and arm."""
import csv, json, os, re, sys, collections

argv = [a for a in sys.argv[1:]]
csv_out = argv[argv.index("--csv") + 1] if "--csv" in argv else None
if csv_out: i = argv.index("--csv"); del argv[i:i + 2]
prod, cells_tsv = argv[0], argv[1]
package = argv[2] if len(argv) > 2 else "19b4ae9"
rows = [l.split("\t") for l in open(cells_tsv, encoding="utf-8").read().splitlines()]
hdr = rows[0]; cells = [dict(zip(hdr, r)) for r in rows[1:] if r and r[0].startswith("P-")]
csv.field_size_limit(10**9)

def ids(spec):
    a, b = spec.split(":") if ":" in spec else (spec, spec)
    return list(range(int(a), int(b) + 1))

tot = collections.Counter(); problems = []; unres_by = collections.Counter(); failed_by = collections.Counter()
reqfail_list = []
print("cell | tasks | folders | exit!=0 (manifest) | result rows / expected | required checks | required failed | FAILED rows | unresolved T2")
for c in cells:
    cid = c["id"]; tasks = int(c["tasks"]); block = int(c["block"]) if c["block"] not in ("", "-") else 0
    args = c["args"].split()
    val = lambda f: args[args.index(f) + 1] if f in args else None
    d = os.path.join(prod, cid)
    folders = sorted(os.path.join(d, x) for x in os.listdir(d)) if os.path.isdir(d) else []
    folders = [f for f in folders if os.path.isdir(f)]
    tags = collections.Counter(re.sub(r"^.*-(t\d{3})-[0-9TZ]+$", r"\1", os.path.basename(f)) for f in folders)
    dup = [t for t, n in tags.items() if n > 1]
    with_man = [f for f in folders if os.path.isfile(os.path.join(f, "manifest.json"))]
    nonzero = 0; res_rows = 0; exp_rows = 0; nreq = 0; nreqf = 0; nfailed = 0; nunres = 0
    for f in folders:
        mp = os.path.join(f, "manifest.json")
        if not os.path.isfile(mp):
            problems.append(f"{cid}/{os.path.basename(f)}: no manifest.json; files: {sorted(os.listdir(f))}")
            continue
        m = json.load(open(mp, encoding="utf-8"))
        ex = (m.get("checks") or {}).get("exit_status", m.get("exit_status"))
        if ex not in (0, None): nonzero += 1
        if m.get("mode") != "prod": problems.append(f"{cid}: mode {m.get('mode')}")
        if (m.get("spec_freeze") or {}).get("sha256", (m.get("spec_freeze") or {}).get("actual")) != "58fe3c9f53dda3c54c0e20c06216ff00e1411541933a23be8a9bbf8f64eead61":
            problems.append(f"{cid}/{os.path.basename(f)}: spec sha differs")
        pc = (m.get("git") or {}).get("package_commit", m.get("package_commit"))
        if pc and not pc.startswith(package): problems.append(f"{cid}: package {pc}")
        rp = os.path.join(f, "results.csv")
        if os.path.isfile(rp):
            rr = list(csv.DictReader(open(rp, encoding="utf-8")))
            res_rows += len(rr)
            nfailed += sum(r.get("arm_status") == "FAILED" for r in rr)
            for r in rr:
                if r.get("arm_status") == "FAILED": failed_by[(cid, r["arm"])] += 1
                if r.get("T2_status") == "unresolved": nunres += 1; unres_by[(cid, r["arm"])] += 1
            arms = (m.get("config") or {}).get("arms")
            pops = sorted({r.get("pop_id") for r in rr}); reps = sorted({r.get("rep_id") for r in rr})
            if arms and pops and reps: exp_rows += len(arms) * len(pops) * len(reps)
        ck = os.path.join(f, "checks.csv")
        if os.path.isfile(ck):
            for r in csv.DictReader(open(ck, encoding="utf-8")):
                if r["required"] == "TRUE":
                    nreq += 1
                    if r["pass"] != "TRUE":
                        nreqf += 1; reqfail_list.append((cid, os.path.basename(f), r["pop_id"], r["rep_id"], r["arm"], r["check"][:110], r["value"][:22], r["reference"][:22], r["abs_diff"][:10], r["tolerance"][:10]))
            nr = (m.get("checks") or {}).get("n_required")
            if nr is not None and int(nr) != sum(1 for r in csv.DictReader(open(ck, encoding="utf-8")) if r["required"] == "TRUE"):
                problems.append(f"{cid}/{os.path.basename(f)}: checks.csv required rows != manifest n_required {nr}")
    # expected design rows from cells.tsv (sampled cells only)
    pops_all = ids(val("--pops")) if val("--pops") else []
    R = int(val("--R")) if val("--R") else 0
    missing_tasks = sorted(set(f"t{i:03d}" for i in range(1, tasks + 1)) - set(tags)) if tasks > 1 else ([] if folders else ["t001"])
    print(f"{cid} | {tasks} | {len(folders)} | {nonzero} | {res_rows} / {exp_rows} | {nreq} | {nreqf} | {nfailed} | {nunres}"
          + (f" | MISSING {missing_tasks}" if missing_tasks else "") + (f" | DUP {dup}" if dup else ""))
    tot["folders"] += len(folders); tot["req"] += nreq; tot["reqfail"] += nreqf; tot["failed"] += nfailed; tot["unres"] += nunres; tot["rows"] += res_rows
print("TOTAL", dict(tot))
print("== problems ==", *problems[:40], sep="\n  ")
print("== failed required checks ==")
for x in reqfail_list[:30]: print("  ", x)
print("== FAILED fits by cell/arm ==", dict(failed_by))
print("== unresolved T2 by cell/arm ==", dict(unres_by))
if csv_out:
    with open(csv_out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("cell,arm,unresolved\n" + "".join("%s,%s,%d\n" % (c, a, v) for (c, a), v in unres_by.items()))
