#!/usr/bin/env python3
"""End-to-end test of convert/delphes_to_hyper.py on a synthetic Delphes-like file.

Each synthetic event is built so the right answer is known by construction: jets are
placed on top of chosen quarks, and everything else is kept well away from them.
Run:  python3 tests/test_convert.py
"""
import json
import os
import subprocess
import sys
import tempfile

import awkward as ak
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pylib"))
import hyperio  # noqa: E402
import matching  # noqa: E402

rng = np.random.default_rng(7)


def dr(e1, p1, e2, p2):
    dp = abs(p1 - p2)
    dp = 2 * np.pi - dp if dp > np.pi else dp
    return float(np.hypot(e1 - e2, dp))


def spread_points(n, min_sep, eta_max=2.0, avoid=()):
    pts = list(avoid)
    while len(pts) < len(avoid) + n:
        c = (rng.uniform(-eta_max, eta_max), rng.uniform(-np.pi, np.pi))
        if all(dr(c[0], c[1], p[0], p[1]) > min_sep for p in pts):
            pts.append(c)
    return pts[len(avoid):]


def make_record(quark_dirs, hadronic=True):
    """Generator record with shower copies of the tops and one W. Returns rows and quark rows."""
    rows = []   # (pid, m1, m2, pt, eta, phi, m)

    def add(pid, m1, m2=-1, pt=50.0, eta=0.0, phi=0.0, m=0.0):
        rows.append((pid, m1, m2, pt, eta, phi, m))
        return len(rows) - 1

    add(2212, -1); add(2212, -1)
    g1 = add(21, 0); g2 = add(21, 1)
    t = add(6, g1, g2, pt=120.0, m=172.5)
    tb = add(-6, g1, g2, pt=110.0, m=172.5)
    add(211, 3, pt=1.0)                                   # padding, shifts indices
    t2 = add(6, t, pt=118.0, m=172.5)                     # recoil copy
    glu = add(21, t, pt=10.0)
    add(5, glu, pt=8.0, eta=3.5); add(-5, glu, pt=7.0, eta=-3.5)   # g -> b bbar, not from a top
    tb2 = add(-6, tb, pt=109.0, m=172.5)
    tb3 = add(-6, tb2, pt=108.0, m=172.5)
    wp = add(24, t2, m=80.4)
    b = add(5, t2, eta=quark_dirs[0][0], phi=quark_dirs[0][1], m=4.7)
    bbar = add(-5, tb3, eta=quark_dirs[3][0], phi=quark_dirs[3][1], m=4.7)
    wm = add(-24, tb3, m=80.4)
    wp2 = add(24, wp, m=80.4)                             # W copy
    q1 = add(2, wp2, eta=quark_dirs[1][0], phi=quark_dirs[1][1])
    q2 = add(-1, wp2, eta=quark_dirs[2][0], phi=quark_dirs[2][1])
    if hadronic:
        # antiquark first in the record: the converter must still give the quark label 5
        q4 = add(-4, wm, eta=quark_dirs[5][0], phi=quark_dirs[5][1])
        q3 = add(3, wm, eta=quark_dirs[4][0], phi=quark_dirs[4][1])
    else:
        q3 = add(11, wm); q4 = add(-12, wm)
    add(2, q1, pt=30.0, eta=quark_dirs[1][0] + 1.0)       # later copy of a decay quark
    for _ in range(20):
        add(211, q2, pt=2.0)
    return rows, [b, q1, q2, bbar, q3, q4]


def make_event(kind):
    """Returns (particle rows, jets [(pt, eta, phi, mass, btag)], expected labels by jet or None if dropped)."""
    qd = spread_points(6, 1.2)
    if kind == "wrap":
        qd[1] = (0.3, 3.10)
        qd = [qd[1]] + spread_points(5, 1.2, avoid=[qd[1]])
        qd = [qd[1], qd[0]] + qd[2:]                       # the wrapped quark is W1q1
    rows, _ = make_record(qd, hadronic=(kind != "semilep"))
    jets, labels = [], []

    def jet_on(iq, pt, btag, off=0.05, eta=None):
        e = qd[iq][0] + off if eta is None else eta
        p = qd[iq][1] - off
        if kind == "wrap" and iq == 1:
            p = -3.10                                     # across the phi boundary, dR ~ 0.1
        jets.append((pt, e, p, 8.0, btag)); labels.append(iq + 1)

    matched = {"full": [0, 1, 2, 3, 4, 5], "wrap": [0, 1, 2, 3, 4, 5], "semilep": [0, 1, 2, 3, 4, 5],
               "oneW": [0, 1, 2, 3, 4], "noW": [0, 1, 3, 4], "fewjets": [0, 1, 2, 3],
               "onebtag": [0, 1, 2, 3, 4, 5], "lowpt": [0, 1, 2, 3, 4, 5], "conflict": [0, 1, 2, 3, 4, 5]}[kind]
    for k, iq in enumerate(matched):
        is_b = iq in (0, 3)
        btag = 1 if is_b else (2 if k % 2 else 0)          # bit 1 set on some light jets: must not count
        if kind == "onebtag" and iq == 3:
            btag = 0
        pt = 200.0 - 20.0 * k
        if kind == "lowpt" and iq == 2:
            pt = 22.0                                     # below the 25 GeV cut -> W1q2 unmatched
        jet_on(iq, pt, btag)
    if kind == "lowpt":
        labels[matched.index(2)] = None                   # that jet is removed by the pT cut
    n_extra = 0 if kind == "fewjets" else 4
    for c in spread_points(n_extra, 0.9, avoid=qd):
        jets.append((rng.uniform(26, 60), c[0], c[1], 5.0, 0)); labels.append(0)
    if kind != "fewjets":
        c = spread_points(1, 0.9, avoid=qd)[0]
        jets.append((90.0, 5.2, c[1], 5.0, 0)); labels.append(None)      # beyond |eta| 5, cut
        jets.append((45.0, 3.2, c[1], 5.0, 0)); labels.append(0)         # forward jet, kept
    if kind == "conflict":
        # move W2q2 to 0.1 from W2q1; its own jet is gone, and W2q1's jet is closer to W2q1
        jets.pop(matched.index(5)); labels.pop(matched.index(5))
        rows, _ = make_record(qd[:5] + [(qd[4][0] + 0.02, qd[4][1] + 0.10)])
    perm = rng.permutation(len(jets))                      # the converter must sort by pT itself
    jets = [jets[i] for i in perm]; labels = [labels[i] for i in perm]

    kept = [(j, l) for j, l in zip(jets, labels) if l is not None]
    kept.sort(key=lambda x: -x[0][0])
    exp_labels = [l for _, l in kept]
    nb = sum(j[4] & 1 for j, _ in kept)
    has_w = (2 in exp_labels and 3 in exp_labels) or (5 in exp_labels and 6 in exp_labels)
    if kind == "semilep" or len(kept) < 6 or nb < 2:
        return rows, jets, None
    return rows, jets, {"labels": exp_labels, "pt": [j[0] for j, _ in kept], "nb": nb, "has_w": has_w, "kind": kind,
                        "all": all(k in exp_labels for k in range(1, 7))}


def main():
    kinds = ["full", "oneW", "noW", "fewjets", "onebtag", "lowpt", "conflict", "wrap", "semilep"]
    events = [make_event(kinds[i % len(kinds)]) for i in range(900)]
    tmp = tempfile.mkdtemp(prefix="ttgen_test_")
    src, dst, meta = (os.path.join(tmp, n) for n in ["delphes.root", "out.root", "out.json"])

    def jag(col, src_idx, dtype):
        return ak.values_astype(ak.Array([[r[col] for r in ev[src_idx]] for ev in events]), dtype)

    with uproot.recreate(src) as f:
        data = {
            "Jet.PT": jag(0, 1, np.float32), "Jet.Eta": jag(1, 1, np.float32), "Jet.Phi": jag(2, 1, np.float32),
            "Jet.Mass": jag(3, 1, np.float32), "Jet.BTag": jag(4, 1, np.uint32),
            "Particle.PID": jag(0, 0, np.int32), "Particle.M1": jag(1, 0, np.int32), "Particle.M2": jag(2, 0, np.int32),
            "Particle.PT": jag(3, 0, np.float32), "Particle.Eta": jag(4, 0, np.float32),
            "Particle.Phi": jag(5, 0, np.float32), "Particle.Mass": jag(6, 0, np.float32),
            "Event.Weight": ak.Array([[1.0 if i % 5 else -1.0] for i in range(len(events))]),
        }
        f.mktree("Delphes", {k: v.type.content for k, v in data.items()})
        f["Delphes"].extend(data)
    cmd = [sys.executable, os.path.join(HERE, "..", "convert", "delphes_to_hyper.py"),
           "--input", src, "--output", dst, "--meta", meta, "--step", "250",
           "--extra", "task=3", "--extra", "grid_mode=reuse"]
    # one event in nine is semileptonic: the default guard must stop the conversion
    guard = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert guard.returncode != 0 and b"generator record" in guard.stderr, "the all-hadronic guard did not fire"
    subprocess.run(cmd + ["--max-truth-fail-frac", "0.5"], check=True)

    expected = [e[2] for e in events if e[2] is not None]
    got = hyperio.read(dst)
    m = json.load(open(meta))
    with uproot.open(dst) as f:
        assert f.classname_of("Delphes") == "TTree", f.classname_of("Delphes")
        keys = set(f["Delphes"].keys())
        njet_branch = f["Delphes"]["njet"].array(library="np")
    need = {"njet", "nbTagged", "allMatchedEvent", "leastOneMatchedW", "jet_pt", "jet_eta", "jet_phi", "jet_e",
            "jet_m", "jet_bTag", "jet_truthmatch", "particle_pid", "particle_pt", "particle_eta", "particle_phi", "particle_m",
            "quark_pid", "quark_pt", "quark_eta", "quark_phi", "quark_m"}
    assert need <= keys, "missing branches: %s" % (need - keys)
    # every event passing the jet and b-tag cuts is written; the matched-W requirement is a flag
    assert len(got["jet_pt"]) == len(expected) == m["n_prematch"], (len(got["jet_pt"]), len(expected), m["n_prematch"])
    assert m["n_selected"] == sum(e["has_w"] for e in expected)
    for i, exp in enumerate(expected):
        assert got["jet_truthmatch"][i].tolist() == exp["labels"], (i, got["jet_truthmatch"][i].tolist(), exp["labels"])
        assert np.allclose(got["jet_pt"][i], exp["pt"]), i
        assert int(got["nbTagged"][i]) == exp["nb"], i
        assert int(got["allMatchedEvent"][i]) == int(exp["all"]), i
        assert int(got["leastOneMatchedW"][i]) == int(exp["has_w"]), i
        assert got["quark_pid"][i].tolist() == [5, 2, -1, -5, 3, -4], got["quark_pid"][i].tolist()
        assert njet_branch[i] == len(exp["labels"]), i
        assert got["particle_pid"][i].tolist() == [6, -6, 24, -24], got["particle_pid"][i].tolist()
        assert np.allclose(got["particle_pt"][i].tolist()[:2], [120.0, 110.0])    # first copies of the tops
    e = np.asarray(ak.flatten(got["jet_e"])); pt = np.asarray(ak.flatten(got["jet_pt"]))
    eta = np.asarray(ak.flatten(got["jet_eta"])); mass = np.asarray(ak.flatten(got["jet_m"]))
    assert np.allclose(e, np.sqrt((pt * np.cosh(eta)) ** 2 + mass ** 2), rtol=1e-5)
    assert pt.min() >= 25 and 3.0 < np.abs(eta).max() < 5.0
    per_kind = len(events) // len(kinds)
    assert m["n_read"] == len(events) and m["n_truth_fail"] == per_kind, m
    assert m["n_negative_weight"] == len(events) // 5, m
    assert m["n_all_matched"] == sum(e["all"] for e in expected)
    assert m["task"] == 3 and m["grid_mode"] == "reuse"
    by_kind = {k: sum(1 for e in expected if e["kind"] == k and e["has_w"]) for k in kinds}
    print("with a matched W, per kind:", by_kind)
    assert by_kind == {"full": per_kind, "oneW": per_kind, "noW": 0, "fewjets": 0, "onebtag": 0,
                       "lowpt": per_kind, "conflict": per_kind, "wrap": per_kind, "semilep": 0}, by_kind
    assert sum(1 for e in expected if e["kind"] == "noW") == per_kind      # written, flagged 0

    # ---- re-matching from the stored quarks ----
    again = matching.relabel(got, 0.4, "closest", "t-first")
    assert ak.to_list(again["jet_truthmatch"]) == ak.to_list(got["jet_truthmatch"])
    swapped = matching.relabel(got, 0.4, "closest", "tbar-first")
    remap = {0: 0, 1: 4, 2: 5, 3: 6, 4: 1, 5: 2, 6: 3}
    assert ak.to_list(swapped["jet_truthmatch"]) == [[remap[x] for x in ev] for ev in ak.to_list(got["jet_truthmatch"])]
    is_conf = np.array([e["kind"] == "conflict" for e in expected])
    last = matching.relabel(got, 0.4, "last-wins", "t-first")["jet_truthmatch"]
    first = matching.relabel(got, 0.4, "first-wins", "t-first")["jet_truthmatch"]
    # in the "conflict" events labels 5 and 6 compete for one jet, and 5 is the closer quark
    assert all(5 in ev and 6 not in ev for ev in ak.to_list(got["jet_truthmatch"][is_conf]))
    assert all(5 in ev and 6 not in ev for ev in ak.to_list(first[is_conf]))
    assert all(6 in ev and 5 not in ev for ev in ak.to_list(last[is_conf]))
    assert ak.to_list(last[~is_conf]) == ak.to_list(got["jet_truthmatch"][~is_conf])
    # greedy: the loser of a conflict takes its next-closest free jet
    q_eta = np.zeros((1, 6)); q_phi = np.array([[0.0, 0.05, 1.0, 2.0, 3.0, -2.0]])
    j_eta = np.array([[0.0, 0.0, 0.0]]); j_phi = np.array([[0.02, 0.30, 1.05]])
    assert matching.labels_for(q_eta, q_phi, j_eta, j_phi, 0.4, "closest").tolist() == [[1, 0, 3]]
    assert matching.labels_for(q_eta, q_phi, j_eta, j_phi, 0.4, "greedy").tolist() == [[1, 2, 3]]
    assert matching.labels_for(q_eta, q_phi, j_eta, j_phi, 0.4, "last-wins").tolist() == [[2, 0, 3]]
    print("test_convert: OK (%d events in, %d written, %d with a matched W, files in %s)"
          % (len(events), len(expected), m["n_selected"], tmp))
    return dst


if __name__ == "__main__":
    main()
