#!/usr/bin/env python3
"""Check the two validation scripts on synthetic HyPER-layout files with a known answer:
two samples drawn from one distribution must agree, and a sample with a 4% jet-energy
shift and a lower b-tag rate must be flagged.   Run:  python3 tests/test_validate.py
"""
import json
import os
import subprocess
import sys
import tempfile

import awkward as ak
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pylib"))
import hyperio  # noqa: E402


def make(path, n, seed, pt_scale=1.0, btag_eff=0.75):
    rng = np.random.default_rng(seed)
    njet = np.clip(6 + rng.poisson(1.3, n), 6, 20)
    tot = int(njet.sum())
    pt = (25 + rng.exponential(55, tot)) * pt_scale
    eta = rng.uniform(-2.5, 2.5, tot)
    phi = rng.uniform(-np.pi, np.pi, tot)
    m = rng.uniform(3, 20, tot)
    un = lambda x: ak.unflatten(x, njet)  # noqa: E731
    order = ak.argsort(un(pt), axis=1, ascending=False)
    pt, eta, phi, m = (un(x)[order] for x in (pt, eta, phi, m))
    # labels: six distinct jets get 1..6, each kept with probability 0.8
    lab = np.zeros(tot, dtype=np.int32)
    starts = np.concatenate([[0], np.cumsum(njet)[:-1]])
    for s, k in zip(starts, njet):
        pos = rng.permutation(k)[:6]
        keep = rng.random(6) < 0.8
        lab[s + pos[keep]] = np.arange(1, 7)[keep]
    is_b = (lab == 1) | (lab == 4)
    btag = np.where(is_b, rng.random(tot) < btag_eff, rng.random(tot) < 0.03).astype(np.int32)
    lab, btag = un(lab), un(btag)
    arrays = {"jet_pt": pt, "jet_eta": eta, "jet_phi": phi, "jet_m": m,
              "jet_e": np.sqrt((pt * np.cosh(eta)) ** 2 + m ** 2), "jet_bTag": btag, "jet_truthmatch": lab,
              "nbTagged": np.asarray(ak.sum(btag, axis=1)),
              "allMatchedEvent": np.asarray(ak.sum(lab > 0, axis=1) == 6).astype(int),
              "leastOneMatchedW": np.ones(n, dtype=int), "eventWeight": np.ones(n)}
    w = hyperio.Writer(path)
    w.write(arrays)
    w.close()


def run(script, *args):
    out = subprocess.run([sys.executable, os.path.join(HERE, "..", "validate", script)] + list(args),
                         check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout
    return out


def main():
    tmp = tempfile.mkdtemp(prefix="ttgen_val_")
    ref, same, shifted = (os.path.join(tmp, n) for n in ["ref.root", "same.root", "shifted.root"])
    make(ref, 120000, 1)
    make(same, 60000, 2)
    make(shifted, 60000, 3, pt_scale=1.04, btag_eff=0.70)

    for tag, new in [("same", same), ("shifted", shifted)]:
        print(run("compare_spectra.py", "--ref", ref, "--new", new, "--out", os.path.join(tmp, "t1_" + tag)).split("\n\n")[1])
        run("c2st.py", "--ref", ref, "--new", new, "--n", "60000", "--out", os.path.join(tmp, "t2_%s.json" % tag))
    t1 = {t: json.load(open(os.path.join(tmp, "t1_%s" % t, "spectra.json"))) for t in ("same", "shifted")}
    t2 = {t: json.load(open(os.path.join(tmp, "t2_%s.json" % t))) for t in ("same", "shifted")}

    worst_same = max(r["chi2"] / r["ndf"] for r in t1["same"]["shapes"])
    pt_shift = [r for r in t1["shifted"]["shapes"] if r["observable"] == "jet pT [GeV]"][0]
    btag_same = [r for r in t1["same"]["rates"] if r["rate"].startswith("b-tag rate, jets matched to b")][0]
    btag_shift = [r for r in t1["shifted"]["rates"] if r["rate"].startswith("b-tag rate, jets matched to b")][0]
    print("tier 1: worst chi2/ndf, same distribution %.2f; jet pT chi2/ndf after shift %.1f" % (worst_same, pt_shift["chi2"] / pt_shift["ndf"]))
    print("tier 1: b-tag rate difference %+.3f (same) %+.3f (shifted, true -0.05)" % (btag_same["difference"], btag_shift["difference"]))
    print("tier 2: AUC same %.4f (null %s), shifted %.4f; top inputs %s"
          % (t2["same"]["auc"], ["%.4f" % v for v in t2["same"]["auc_null"]], t2["shifted"]["auc"],
             [n for n, _ in t2["shifted"]["top_features"][:3]]))
    assert worst_same < 2.5, worst_same
    assert pt_shift["chi2"] / pt_shift["ndf"] > 10
    assert abs(btag_same["sigma"]) < 4 and abs(btag_shift["difference"] + 0.05) < 0.01
    assert abs(t2["same"]["auc"] - 0.5) < 0.02 and t2["shifted"]["auc"] > 0.55
    assert os.path.getsize(os.path.join(tmp, "t1_same", "spectra.pdf")) > 10000
    print("test_validate: OK (files in %s)" % tmp)
    return tmp


if __name__ == "__main__":
    main()
