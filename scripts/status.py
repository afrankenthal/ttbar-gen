#!/usr/bin/env python3
"""Summarise finished tasks from their JSON files and list the ones still missing.

usage: status.py <raw dir>   (reads <raw dir>/*.json and <raw dir>/*/*.json)
       status.py <raw dir> [--tasks N] [--first A --last B] [--target-selected N] [--missing FILE]

The pilot numbers to read off: selection efficiency (HyPER: 14.67%), fully matched
fraction of selected (HyPER: 26.7%), seconds per event, and peak scratch use.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402

HYPER_EFF = 8801060.0 / 60000000.0
HYPER_ALL_MATCHED = (2116231 + 117487 + 116960) / 8801060.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raw")
    ap.add_argument("--tasks", type=int, help="total number of tasks expected (ids 1..N)")
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--last", type=int)
    ap.add_argument("--target-selected", type=int, default=79212990)
    ap.add_argument("--missing", help="write missing task ids here, for submit.sh --list")
    a = ap.parse_args()

    metas = []
    for path in hyperio.task_files(a.raw, ".json"):
        with open(path) as f:
            metas.append(json.load(f))
    if not metas:
        raise SystemExit("status: no finished tasks in %s" % a.raw)

    tot = lambda k: sum(m.get(k, 0) for m in metas)  # noqa: E731
    n_read, n_sel, n_all = tot("n_read"), tot("n_selected"), tot("n_all_matched")
    short = [m["task"] for m in metas if m.get("n_requested") and m["n_read"] < 0.98 * m["n_requested"]]
    eff = n_sel / float(n_read)
    print("finished tasks        %d" % len(metas))
    print("events generated      %d" % n_read)
    print("events selected       %d   efficiency %.2f%%  (HyPER %.2f%%)" % (n_sel, 100 * eff, 100 * HYPER_EFF))
    print("fully matched         %d   %.1f%% of selected  (HyPER %.1f%%)"
          % (n_all, 100.0 * n_all / max(n_sel, 1), 100 * HYPER_ALL_MATCHED))
    print("top 1 / top 2 matched %.1f%% / %.1f%% of selected  (HyPER 57.8%% / 61.8%%)"
          % (100.0 * tot("n_top1_matched") / max(n_sel, 1), 100.0 * tot("n_top2_matched") / max(n_sel, 1)))
    present = [sum(m.get("n_label_present", [0] * 6)[k] for m in metas) for k in range(6)]
    print("label present, 1..6   " + "  ".join("%.1f%%" % (100.0 * x / max(n_sel, 1)) for x in present)
          + "  (HyPER 88.5 84.2 78.9 89.0 81.3 86.7)")
    print("stored for re-matching %d events (passed the jet and b-tag cuts)" % tot("n_prematch"))
    print("negative weights      %.1f%% of generated" % (100.0 * tot("n_negative_weight") / n_read))
    print("truth-record failures %d" % tot("n_truth_fail"))
    timed = [m for m in metas if "seconds_generation" in m]
    if timed:
        sec = sum(m["seconds_generation"] + m["seconds_delphes"] + m.get("seconds_convert", 0) for m in timed)
        ev = sum(m["n_read"] for m in timed)
        cores = max(m.get("cores", 1) for m in timed)
        print("wall time per event   %.3f s on %d cores  (generation %.0f%%, Delphes %.0f%%)"
              % (sec / ev, cores, 100.0 * sum(m["seconds_generation"] for m in timed) / sec,
                 100.0 * sum(m["seconds_delphes"] for m in timed) / sec))
        print("longest task          %.1f h" % (max(m["seconds_generation"] + m["seconds_delphes"] for m in timed) / 3600.0))
        print("peak scratch per task %.1f GB (after generation; the uncompressed HepMC file is larger while it is written)"
              % (max(m.get("scratch_after_generation_kb", 0) for m in timed) / 1e6))
        need_gen = a.target_selected / eff
        print("for %d selected       %.0fM generated, about %.0f core-hours"
              % (a.target_selected, need_gen / 1e6, need_gen * sec / ev * cores / 3600.0))
    if short:
        print("WARNING: %d tasks read fewer events than requested: %s" % (len(short), short[:20]))

    last = a.last or a.tasks
    if last:
        done = set(m["task"] for m in metas if "task" in m)
        missing = [t for t in range(a.first, last + 1) if t not in done]
        print("missing tasks         %d of %d" % (len(missing), last - a.first + 1))
        if a.missing:
            with open(a.missing, "w") as f:
                f.write("".join("%d\n" % t for t in missing))
            print("wrote %s" % a.missing)


if __name__ == "__main__":
    main()
