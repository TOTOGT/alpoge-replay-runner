#!/usr/bin/env python3
"""Fail unless each expected line appears in a log and no violation is reported.
Usage: gate_py.py LOG KEY [KEY ...]   (each KEY must appear in LOG; every
'violations=N' in LOG must have N == 0; 'MISMATCH' must not appear)."""
import re, sys
log = open(sys.argv[1]).read(); bad = []
for k in sys.argv[2:]:
    if k not in log: bad.append(f"missing expected output: {k!r}")
for n in re.findall(r"violations=(\d+)", log):
    if int(n) != 0: bad.append(f"violations={n}")
if "MISMATCH" in log: bad.append("MISMATCH reported")
for b in bad: print("::error::", b)
sys.exit(1 if bad else 0)
