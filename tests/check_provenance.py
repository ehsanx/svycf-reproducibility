#!/usr/bin/env python3
"""Check this copy against config/provenance.csv, and the simulation code against the recorded runs.

    python tests/check_provenance.py                          # 1 only
    python tests/check_provenance.py <target>/production      # 1 and 2 (after scripts/simulation/assemble_production.py)

1. Every file listed in config/provenance.csv has its recorded SHA-256, and no other file is present. Only what a run
   or setup writes locally is exempt -- the paths .gitignore excludes (outputs/, data/, the R libraries, the frozen
   kappa_AE copy, PACKAGE-COMMIT.txt, config/slurm.env, Python caches) and the .git folder.
   Files whose status is "identical" must have the SHA-256 of their original (recorded_sha256). Files whose status is
   "comments-edited" must have the executable code of their original: the SHA-256 of the file with every R comment
   removed (code_sha256, recomputed here) must equal recorded_code_sha256. Any other change to code or strings fails.
2. Every recorded production run stored the SHA-256 of the engine files it ran (engine_files_sha256 in its
   manifest.json). Each is compared with the file of the same name in src/simulation/: the same bytes, the same
   executable code with comments edited (via the original's recorded_sha256 and the code check of part 1), or a
   difference, named with its status.
Exit status 0 only if every check of part 1 holds.
"""
import collections
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL_DIRS = ("outputs/", "data/", ".git/", "src/simulation/renv/library/", "src/simulation/renv/sandbox/",
              "src/simulation/renv/staging/", "src/simulation/runs/frozen/", "src/application/library/")
LOCAL_FILES = {"PACKAGE-COMMIT.txt", "config/slurm.env", ".Rhistory", "config/provenance.csv"}


def local_only(rel):
    parts = rel.split("/")
    return (rel.startswith(LOCAL_DIRS) or rel in LOCAL_FILES or "__pycache__" in parts or rel.endswith(".pyc")
            or parts[-1] == ".Rhistory")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def r_code_sha256(path):
    """SHA-256 of an R file's executable code: comments (a '#' outside a string, to the end of the line) removed,
    trailing blanks removed, blank lines dropped -- the same canonical form as code_sha256() in R/grade_run.R."""
    text = Path(path).read_bytes().decode("utf-8").replace("\r\n", "\n")
    lines, quote, i, line, col0 = text.split("\n"), None, 0, 0, 0
    cuts = {}
    while i < len(text):
        c = text[i]
        if c == "\n":
            line += 1; col0 = i + 1; i += 1; continue
        if quote:
            if c == "\\":
                i += 2; continue
            if c == quote:
                quote = None
        elif c in "'\"`":
            quote = c
        elif c == "#":
            cuts[line] = i - col0
            j = text.find("\n", i)
            i = len(text) if j < 0 else j
            continue
        i += 1
    for ln, col in cuts.items():
        lines[ln] = lines[ln][:col]
    keep = [l.rstrip() for l in lines if l.strip()]
    return hashlib.sha256(("\n".join(keep) + "\n").encode("utf-8")).hexdigest()


prov = {r["path"]: r for r in csv.DictReader(open(ROOT / "config" / "provenance.csv", encoding="utf-8"))}
bad = [p for p, r in prov.items() if not (ROOT / p).is_file() or sha(ROOT / p) != r["sha256"]]
present = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file()}
extra = sorted(p for p in present - set(prov) if not local_only(p))
origin = []
status_n = collections.Counter(r["status"] for r in prov.values())
for p, r in prov.items():
    if p in bad:
        continue
    if r["status"] == "identical" and r["recorded_sha256"] and r["sha256"] != r["recorded_sha256"]:
        origin.append("%s: not the bytes of its original" % p)
    if r["status"] == "comments-edited":
        got = r_code_sha256(ROOT / p)
        if not (got == r["code_sha256"] == r["recorded_code_sha256"]):
            origin.append("%s: executable code differs from its original" % p)
print("1. %d files listed (%s); %d differ or are missing; %d unexpected; %d fail the check against their original"
      % (len(prov), ", ".join("%d %s" % (n, s) for s, n in sorted(status_n.items())), len(bad), len(extra), len(origin)))
for p in bad:
    print("   differs or missing:", p)
for p in extra:
    print("   unexpected:", p)
for p in origin:
    print("   " + p)
status = 1 if bad or extra or origin else 0

if len(sys.argv) > 1:
    prod = Path(sys.argv[1])
    recorded = collections.defaultdict(collections.Counter)
    runs = 0
    for m in prod.glob("*/*/manifest.json"):
        ef = json.load(open(m, encoding="utf-8")).get("engine_files_sha256") or {}
        runs += 1
        for k, v in ef.items():
            if isinstance(v, str):
                recorded[k][v] += 1
    print("2. engine files hashed by %d recorded run folders, compared with src/simulation/:" % runs)
    for k in sorted(recorded):
        row = prov.get("src/simulation/" + k, {})
        mine = ROOT / "src" / "simulation" / k
        h = sha(mine) if mine.is_file() else "missing"
        hashes = set(recorded[k])
        if hashes == {h}:
            verdict = "identical"
        elif row.get("status") == "comments-edited" and hashes == {row.get("recorded_sha256")}:
            verdict = "identical code (comments edited)"
        else:
            verdict = "differs (%s)" % row.get("status", "not in release")
        print("   %-28s %-34s in %d runs%s" % (k, verdict, sum(recorded[k].values()),
                                               "" if len(hashes) == 1 else "; %d distinct recorded hashes" % len(hashes)))
sys.exit(status)
