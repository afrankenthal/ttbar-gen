#!/usr/bin/env python3
"""Merge per-task files into shards in HyPER's layout, in task order.

usage: merge.py <raw dir> <out dir> [--shard-size N] [--holdout N] [--prefix NAME]
                [--rematch --match-dr R --conflict RULE --top-order ORDER] [--keep-quarks]

Per-task files hold every event that passed the jet and b-tag cuts. This step keeps
the ones with at least one matched W and drops the stored quarks, which is HyPER's
layout. --rematch first recomputes the labels from the stored quarks under the given
rule (see pylib/matching.py), so a different matching convention needs no new events.

--holdout N writes the first N events to <prefix>_holdout.root, for the "new test
set" of the closure test; the rest go to <prefix>_train_000.root, _001, ...
A manifest JSON records the settings and the event counts.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402
import matching  # noqa: E402


def take(arrays, start, stop):
    return {k: v[start:stop] for k, v in arrays.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("raw")
    ap.add_argument("out")
    ap.add_argument("--shard-size", type=int, default=7921299, help="events per training shard (default: HyPER's train size)")
    ap.add_argument("--holdout", type=int, default=0)
    ap.add_argument("--prefix", default="ttbar_new")
    ap.add_argument("--buffer", type=int, default=200000, help="events held in memory before each write")
    ap.add_argument("--rematch", action="store_true", help="recompute labels from the stored quarks")
    ap.add_argument("--match-dr", type=float, default=0.4)
    ap.add_argument("--conflict", choices=matching.CONFLICTS, default="closest")
    ap.add_argument("--top-order", choices=matching.TOP_ORDERS, default="t-first")
    ap.add_argument("--keep-quarks", action="store_true", help="keep the quark_* branches in the output")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(a.raw, "*.root")))
    if not files:
        raise SystemExit("merge: no .root files in %s" % a.raw)

    plan = []                                       # (name, capacity)
    if a.holdout:
        plan.append(("%s_holdout.root" % a.prefix, a.holdout))
    state = {"i": -1, "writer": None, "room": 0, "name": None, "pending": [], "n_pending": 0}
    manifest = {"shards": [], "source_files": len(files), "rematch": a.rematch}
    if a.rematch:
        manifest.update({"match_dr": a.match_dr, "conflict": a.conflict, "top_order": a.top_order})
    rng = np.random.default_rng(12345)
    seen = 0

    def flush():
        if state["pending"]:
            state["writer"].write(hyperio.concatenate(state["pending"]))
            state["pending"], state["n_pending"] = [], 0

    def close():
        if state["writer"] is not None:
            flush()
            state["writer"].close()
            manifest["shards"].append({"file": state["name"], "events": state["writer"].entries})
            print("merge: %s  %d events" % (state["name"], state["writer"].entries))

    def next_shard():
        close()
        state["i"] += 1
        if state["i"] < len(plan):
            name, cap = plan[state["i"]]
        else:
            name, cap = "%s_train_%03d.root" % (a.prefix, state["i"] - len(plan)), a.shard_size
        state.update(writer=hyperio.Writer(os.path.join(a.out, name)), room=cap, name=name)

    next_shard()
    for path in files:
        arrays = hyperio.read(path)
        seen += len(arrays["jet_pt"])
        if a.rematch:
            if "quark_eta" not in arrays:
                raise SystemExit("merge: %s has no stored quarks; --rematch needs per-task files" % path)
            arrays = matching.relabel(arrays, a.match_dr, a.conflict, a.top_order, rng)
        arrays = matching.select(arrays, np.asarray(arrays["leastOneMatchedW"]) > 0)
        if not a.keep_quarks:
            arrays = hyperio.drop_quarks(arrays)
        n, pos = len(arrays["jet_pt"]), 0
        while pos < n:
            if state["room"] == 0:
                next_shard()
            k = min(n - pos, state["room"])
            state["pending"].append(take(arrays, pos, pos + k))
            state["n_pending"] += k
            state["room"] -= k
            pos += k
            if state["n_pending"] >= a.buffer:
                flush()
    close()
    manifest["events"] = sum(s["events"] for s in manifest["shards"])
    manifest["events_before_matched_w_requirement"] = seen
    with open(os.path.join(a.out, "%s_manifest.json" % a.prefix), "w") as f:
        json.dump(manifest, f, indent=1)
    print("merge: %d events in %d shards" % (manifest["events"], len(manifest["shards"])))


if __name__ == "__main__":
    main()
