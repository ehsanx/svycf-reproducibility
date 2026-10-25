#!/usr/bin/env python3
"""Assemble the recorded per-replicate simulation output from the two release archives, and verify it.

    python scripts/simulation/assemble_production.py <folder holding the two archives> <target folder>

The archives (names, sizes and checksums in config/assets.csv) are checked first. Both are then extracted into
<target>/production/: simulation-production-outputs.tar.gz holds every run folder of the 41 production cells except
task 7 of P-e3-ctrl-m400, whose first run timed out before writing anything; simulation-task7-recovery.tar.gz holds
that task's rerun (same keyed random-number streams, same code). Every extracted file is checked against the
archive's MANIFEST-*.sha256, and the assembled tree is counted: run folders, result rows and cells.
The result is the tree the grading pass read: tests/audit_production.py and src/simulation/runs/run_grade.R
(--in <target>/production) run on it.
"""
import csv
import hashlib
import os
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARCHIVES = ["simulation-production-outputs.tar.gz", "simulation-task7-recovery.tar.gz"]


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(src, target):
    expected = {r["file"]: r for r in csv.DictReader(open(ROOT / "config" / "assets.csv", encoding="utf-8"))}
    target = Path(target).resolve()
    for name in ARCHIVES:
        p = Path(src) / name
        if not p.is_file():
            sys.exit("missing %s" % p)
        got = sha256_file(p)
        if got != expected[name]["sha256"]:
            sys.exit("%s: SHA-256 %s, config/assets.csv records %s" % (name, got, expected[name]["sha256"]))
        print("%s: SHA-256 verified" % name)
    for name in ARCHIVES:
        with tarfile.open(Path(src) / name, "r:gz") as t:
            for m in t:
                parts = Path(m.name).parts
                if not m.isfile() or m.name.startswith("/") or ".." in parts or parts[0] != "production":
                    sys.exit("%s: unexpected member %s" % (name, m.name))
                dest = target / m.name
                if dest.exists():
                    sys.exit("%s: %s already exists (extract into an empty folder)" % (name, dest))
                dest.parent.mkdir(parents=True, exist_ok=True)
                with t.extractfile(m) as fin, open(dest, "wb") as fout:
                    fout.write(fin.read())
        man = target / "production" / ("MANIFEST-%s.sha256" % name.replace(".tar.gz", ""))
        bad = n = 0
        for line in open(man, encoding="utf-8"):
            h, rel = line.rstrip("\n").split("  ", 1)
            n += 1
            if sha256_file(target / rel) != h:
                bad += 1
                print("  differs:", rel)
        if bad:
            sys.exit("%s: %d of %d files differ from %s" % (name, bad, n, man.name))
        print("%s: %d files extracted and verified" % (name, n))
    prod = target / "production"
    cells = sorted(d.name for d in prod.iterdir() if d.is_dir())
    folders = [f for c in cells for f in (prod / c).iterdir() if (f / "manifest.json").is_file()]
    rows = 0
    csv.field_size_limit(10 ** 9)
    for f in folders:
        r = f / "results.csv"
        if r.is_file():
            with open(r, encoding="utf-8", newline="") as fh:
                rows += sum(1 for _ in csv.reader(fh)) - 1
    print("assembled %s: %d cells, %d run folders, %d result rows" % (prod, len(cells), len(folders), rows))
    with open(ROOT / "config" / "simulation" / "production-cells.tsv", encoding="utf-8") as fh:
        want = {line.split("\t")[0] for line in fh.read().splitlines()[1:] if line}
    if set(cells) != want:
        sys.exit("cells differ from config/simulation/production-cells.tsv: %s" % sorted(set(cells) ^ want))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
