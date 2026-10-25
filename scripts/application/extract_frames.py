#!/usr/bin/env python3
"""Write the analysis-frame sizes of the SPA surveys as CSV, from the run log that 02-resample.R appends.

    python scripts/application/extract_frames.py outputs/application/RESAMPLE-RUNS.md > outputs/application/frames.csv

The census pseudo-reference run (--reference) logs its frame as "<n> records; <n> facilities; <n> excluded".
One row per survey: survey, records, facilities, excluded (records with a missing exposure or outcome).
"""
import re
import sys

SURVEYS = ["MW_SPA_2013-14", "HT_SPA_2017-18", "HT_SPA_2013"]


def frames(text):
    out = []
    for code in SURVEYS:
        m = re.search(r"\| %s \| reference \|[^\n]*?(\d+) records; (\d+) facilities; (\d+) excluded" % re.escape(code), text)
        if not m:
            raise SystemExit("%s: no reference run in the log" % code)
        out.append((code,) + tuple(int(x) for x in m.groups()))
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    rows = frames(open(sys.argv[1], encoding="utf-8").read())
    print("survey,records,facilities,excluded")
    for r in rows:
        print("%s,%d,%d,%d" % r)
