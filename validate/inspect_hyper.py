#!/usr/bin/env python3
"""Print the layout and conventions of a HyPER-layout file.

Run it on the published ttbar_train.root (or ttbar_reduced.root) before production:
it answers the questions the converter had to assume, and prints a suggested setting
where it can. Run it on a converted file too and compare the two printouts.

usage: inspect_hyper.py <file> [--entries N]
"""
import argparse
import collections
import os
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402


def frac(x):
    return 100.0 * float(np.mean(x)) if len(x) else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--entries", type=int, default=500000)
    a = ap.parse_args()

    with uproot.open(a.file) as f:
        print("objects:", {k: f.classname_of(k) for k in f.keys()})
        tree = f[hyperio.TREE]
        print("entries: %d" % tree.num_entries)
        print("branches:")
        for k in tree.keys():
            print("  %-20s %s" % (k, tree[k].typename))
    d = hyperio.read(a.file, entry_stop=a.entries)
    n = len(d["jet_pt"])
    lab = d["jet_truthmatch"]
    njet = np.asarray(ak.num(d["jet_pt"]))
    print("\n--- selection (first %d events) ---" % n)
    print("njet: min %d max %d mean %.2f" % (njet.min(), njet.max(), njet.mean()))
    print("jet pT min %.2f   |eta| max %.3f" % (ak.min(d["jet_pt"]), ak.max(abs(d["jet_eta"]))))
    ncentral = np.asarray(ak.sum(abs(d["jet_eta"]) < 2.5, axis=1))
    print("jets with |eta| < 2.5 per event: min %d;  events with fewer than 6 of them: %.2f%%"
          % (ncentral.min(), frac(ncentral < 6)))
    print("=> above 0%: the >= 6 jet requirement counted all jets;  0%: it may have counted central jets only")
    print("jets sorted by pT in %.2f%% of events"
          % frac(np.asarray(ak.all(d["jet_pt"][:, :-1] >= d["jet_pt"][:, 1:], axis=1))))
    nb = np.asarray(ak.sum(d["jet_bTag"] > 0, axis=1))
    print("b-tagged jets: |eta| max %.3f" % ak.max(abs(d["jet_eta"][d["jet_bTag"] > 0])))
    print("b-tagged jets: min %d mean %.2f;  jet_bTag values %s" % (nb.min(), nb.mean(), np.unique(np.asarray(ak.flatten(d["jet_bTag"])))[:8]))
    if "nbTagged" in d:
        print("nbTagged equals count of jet_bTag>0 in %.2f%% of events" % frac(np.asarray(d["nbTagged"]) == nb))

    print("\n--- labels ---")
    has = np.stack([np.asarray(ak.any(lab == k, axis=1)) for k in range(1, 7)], axis=1)
    dup = np.stack([np.asarray(ak.sum(lab == k, axis=1)) > 1 for k in range(1, 7)], axis=1)
    print("label present, 1..6: " + "  ".join("%.1f%%" % frac(has[:, k]) for k in range(6)))
    print("a label used on more than one jet: %.3f%% of events" % frac(dup.any(axis=1)))
    top1, top2 = has[:, :3].all(axis=1), has[:, 3:].all(axis=1)
    w1, w2 = has[:, 1:3].all(axis=1), has[:, 4:6].all(axis=1)
    print("top 1 fully matched %.2f%%   top 2 fully matched %.2f%%   both %.2f%%" % (frac(top1), frac(top2), frac(top1 & top2)))
    print("W 1 matched %.2f%%   W 2 matched %.2f%%   at least one %.2f%%" % (frac(w1), frac(w2), frac(w1 | w2)))
    if "allMatchedEvent" in d:
        print("allMatchedEvent agrees with labels in %.2f%% of events" % frac((np.asarray(d["allMatchedEvent"]) > 0) == (top1 & top2)))
    asym = frac(top1) - frac(top2)
    print("=> top-1/top-2 asymmetry %.2f points: %s" % (
        asym, "the two label blocks are NOT interchangeable; TOP_ORDER=random would not reproduce this" if abs(asym) > 0.5
        else "label blocks look symmetric"))
    for k, name in [(1, "b1"), (4, "b2"), (2, "W1q1"), (3, "W1q2"), (5, "W2q1"), (6, "W2q2")]:
        tagged = ak.flatten(d["jet_bTag"][lab == k]) > 0
        print("  jets labelled %d (%s): b-tagged %.1f%%" % (k, name, frac(np.asarray(tagged))))
    both = w1
    q1 = np.asarray(ak.firsts(d["jet_pt"][lab == 2])[both]); q2 = np.asarray(ak.firsts(d["jet_pt"][lab == 3])[both])
    print("W1q1 jet has higher pT than W1q2 jet in %.1f%% of events (50%% = unordered)" % frac(q1 > q2))
    lead = np.asarray(ak.firsts(lab[lab > 0]))
    print("highest-pT matched jet belongs to top 1 in %.1f%% of events" % frac(lead[lead > 0] <= 3))

    if "particle_pid" in d:
        print("\n--- stored particles ---")
        npart = np.asarray(ak.num(d["particle_pid"]))
        print("particles per event: min %d max %d" % (npart.min(), npart.max()))
        pats = collections.Counter(tuple(x) for x in ak.to_list(d["particle_pid"][:20000]))
        print("most common PID sequences (first 20000 events):")
        for pat, c in pats.most_common(8):
            print("  %6d  %s" % (c, list(pat)))
        first_top = ak.firsts(d["particle_pid"][abs(d["particle_pid"]) == 6])
        ft = np.asarray(ak.fill_none(first_top, 0))
        if np.any(ft != 0):
            print("first top in the list is t (not tbar) in %.1f%% of events" % frac(ft[ft != 0] == 6))
            tops_pt = d["particle_pt"][abs(d["particle_pid"]) == 6]
            two = np.asarray(ak.num(tops_pt)) >= 2
            a, b = np.asarray(tops_pt[two][:, 0]), np.asarray(tops_pt[two][:, 1])
            same = np.isclose(a, b, rtol=1e-4)
            print("pT of the first top vs the second: higher %.1f%%, lower %.1f%%, equal %.1f%%"
                  % (frac((a > b) & ~same), frac((a < b) & ~same), frac(same)))
            print("=> a large 'equal' share means the stored tops are matrix-element tops (Born-like events")
            print("   have back-to-back tops), not tops after the shower has recoiled them")
            print("=> t first ~100%: TOP_ORDER=t-first;  ~0%: tbar-first;  first top higher ~100%: pt")
        else:
            print("no top quarks stored; the top ordering has to be inferred from the label asymmetry above")
    if "eventWeight" in d:
        print("\nnegative event weights: %.1f%%" % frac(np.asarray(d["eventWeight"]) < 0))


if __name__ == "__main__":
    main()
