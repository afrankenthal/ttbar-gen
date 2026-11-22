#!/usr/bin/env python3
"""Try every matching rule on the stored quarks and print the label rates next to HyPER's.

usage: rematch_scan.py 'raw/*.root' [--entries N] [--match-dr R]

Per-task files keep the six decay quarks, so labels can be recomputed without new
events. HyPER's published file shows top 2 matched more often than top 1; this scan
shows which rule, if any, reproduces that. Then: merge.py --rematch --conflict <rule>.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402
import matching  # noqa: E402

HYPER = {"top1": 57.75, "top2": 61.82, "all": 26.80, "labels": [88.5, 84.2, 78.9, 89.0, 81.3, 86.7]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", help="per-task files (glob)")
    ap.add_argument("--entries", type=int, default=1000000)
    ap.add_argument("--match-dr", type=float, default=0.4)
    a = ap.parse_args()

    chunks, n = [], 0
    for f in sorted(glob.glob(a.files)):
        chunks.append(hyperio.read(f))
        n += len(chunks[-1]["jet_pt"])
        if n >= a.entries:
            break
    if not chunks:
        raise SystemExit("rematch_scan: no files match %s" % a.files)
    arrays = hyperio.concatenate(chunks) if len(chunks) > 1 else chunks[0]
    if "quark_eta" not in arrays:
        raise SystemExit("rematch_scan: these files have no stored quarks (merged files drop them)")
    print("%d events passing the jet and b-tag cuts, dR < %.2f" % (n, a.match_dr))
    print("%-24s %9s %7s %7s %7s   %s" % ("rule", "matched W", "top 1", "top 2", "both", "label present, 1..6"))
    print("%-24s %9s %6.1f%% %6.1f%% %6.1f%%   %s" % ("HyPER (published file)", "", HYPER["top1"], HYPER["top2"],
                                                    HYPER["all"], "  ".join("%.1f" % x for x in HYPER["labels"])))
    rng = np.random.default_rng(1)
    for top_order in ["t-first", "random"]:
        for conflict in matching.CONFLICTS:
            s = matching.summary(matching.relabel(arrays, a.match_dr, conflict, top_order, rng))
            sel = float(max(s["n_selected"], 1))
            print("%-24s %8.1f%% %6.1f%% %6.1f%% %6.1f%%   %s"
                  % ("%s, %s" % (conflict, top_order), 100 * s["n_selected"] / float(s["n_events"]),
                     100 * s["n_top1_matched"] / sel, 100 * s["n_top2_matched"] / sel, 100 * s["n_all_matched"] / sel,
                     "  ".join("%.1f" % (100 * x / sel) for x in s["n_label_present"])))
    print("\n'matched W' is the share of stored events kept by the matched-W requirement; the other")
    print("columns are relative to the kept events, as in HyPER's file.")


if __name__ == "__main__":
    main()
