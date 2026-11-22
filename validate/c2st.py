#!/usr/bin/env python3
"""Tier 2 closure: can a classifier tell the new sample from the HyPER files?

usage: c2st.py --ref ttbar_train.root --new 'merged/ttbar_new_train_*.root'
               [--n N] [--with-labels] [--out FILE]

Trains a gradient-boosted classifier on N events from each side and reports the AUC
on held-out events. The same procedure on two halves of the HyPER file gives the AUC
expected when the distributions are identical; the new-vs-HyPER AUC should sit inside
that spread. If it does not, the permutation-importance list says which inputs give
the samples away. This is a diagnostic: the pass/fail closure test is tier 3 (train
on one sample, test on the other).

Inputs per event: jet and b-tag multiplicity, HT, and pT, eta, phi, mass, b-tag of
the 10 leading jets. --with-labels adds the truth-match label of each jet, which
tests the matching convention as well.
"""
import argparse
import glob
import json
import os
import sys

import awkward as ak
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402

NJ = 10


def load(pattern, entries):
    chunks, n = [], 0
    for f in sorted(glob.glob(pattern)):
        chunks.append(hyperio.read(f, entry_stop=entries - n))
        n += len(chunks[-1]["jet_pt"])
        if n >= entries:
            break
    if not chunks:
        raise SystemExit("c2st: no files match %s" % pattern)
    return hyperio.concatenate(chunks) if len(chunks) > 1 else chunks[0]


def features(d, with_labels):
    def pad(x):
        return np.asarray(ak.fill_none(ak.pad_none(ak.values_astype(x, np.float64), NJ, axis=1, clip=True), np.nan))

    cols = [np.asarray(ak.num(d["jet_pt"]), dtype=float), np.asarray(ak.sum(d["jet_bTag"] > 0, axis=1), dtype=float),
            np.asarray(ak.sum(d["jet_pt"], axis=1), dtype=float)]
    names = ["njet", "nbtag", "HT"]
    fields = ["pt", "eta", "phi", "m", "bTag"] + (["truthmatch"] if with_labels else [])
    for fld in fields:
        p = pad(d["jet_" + fld])
        for j in range(NJ):
            cols.append(p[:, j])
            names.append("jet%d_%s" % (j + 1, fld))
    return np.stack(cols, axis=1), names


def auc_of(x0, x1, seed, importance=False, names=None):
    rng = np.random.default_rng(seed)
    x = np.concatenate([x0, x1])
    y = np.concatenate([np.zeros(len(x0)), np.ones(len(x1))])
    idx = rng.permutation(len(x))
    cut = int(0.7 * len(x))
    tr, te = idx[:cut], idx[cut:]
    clf = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, random_state=seed)
    clf.fit(x[tr], y[tr])
    auc = roc_auc_score(y[te], clf.predict_proba(x[te])[:, 1])
    top = []
    if importance:
        sub = te[:50000]
        imp = permutation_importance(clf, x[sub], y[sub], scoring="roc_auc", n_repeats=3, random_state=seed)
        order = np.argsort(-imp.importances_mean)[:10]
        top = [(names[i], float(imp.importances_mean[i])) for i in order]
    return float(auc), top


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--n", type=int, default=500000, help="events per side")
    ap.add_argument("--with-labels", action="store_true")
    ap.add_argument("--null-repeats", type=int, default=3)
    ap.add_argument("--out", default="closure_tier2.json")
    a = ap.parse_args()

    ref = load(a.ref, 2 * a.n)
    new = load(a.new, a.n)
    xr, names = features(ref, a.with_labels)
    xn, _ = features(new, a.with_labels)
    n = min(a.n, len(xn), len(xr) // 2)
    if n < a.n:
        print("c2st: using %d events per side (the files are smaller than --n)" % n)
    rng = np.random.default_rng(0)
    xr = xr[rng.permutation(len(xr))]
    xn = xn[rng.permutation(len(xn))[:n]]

    null = [auc_of(*np.array_split(xr[rng.permutation(len(xr))[:2 * n]], 2), seed=s)[0] for s in range(a.null_repeats)]
    auc, top = auc_of(xr[:n], xn, seed=100, importance=True, names=names)
    spread = max(float(np.std(null)), 1e-4)
    excess = (auc - float(np.mean(null))) / spread

    print("AUC, HyPER vs HyPER (identical distributions): %s  mean %.4f" % (["%.4f" % v for v in null], np.mean(null)))
    print("AUC, HyPER vs new:                             %.4f  (%.1f null standard deviations above the null mean)" % (auc, excess))
    print("inputs that separate the samples most (drop in AUC when shuffled):")
    for name, val in top:
        print("  %-18s %+.4f" % (name, val))
    with open(a.out, "w") as f:
        json.dump({"n_per_side": n, "with_labels": a.with_labels, "auc_null": null, "auc": auc,
                   "null_sigmas_above": excess, "top_features": top}, f, indent=1)
    print("wrote %s" % a.out)


if __name__ == "__main__":
    main()
