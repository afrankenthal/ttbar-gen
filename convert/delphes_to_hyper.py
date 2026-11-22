#!/usr/bin/env python3
"""Convert a Delphes file of all-hadronic ttbar events to the HyPER dataset layout.

Steps, per event:
  1. keep jets with pT >= --jet-pt-min and |eta| < --jet-abs-eta-max
     (HyPER's files hold jets out to |eta| = 4.98, so the default applies no eta cut)
  2. require >= 6 such jets, >= 2 of them b-tagged
  3. find the six decay quarks (b, q, qbar of each top) in the generator record
  4. match each quark to its closest jet if dR < --match-dr (pylib/matching.py)
Every event passing steps 1-3 is written, with its labels, the flag leastOneMatchedW
and the six quarks. scripts/merge.py keeps the events with a matched W and drops the
quarks, which gives HyPER's layout; because the quarks are kept here, labels can be
recomputed under a different matching rule without generating events again.
Labels follow HyPER: 0 unmatched, 1 = b of top 1, 2 and 3 = quark and antiquark of
its W, 4, 5, 6 the same for top 2.

Assumptions that the HyPER record does not spell out are options or marked ASSUMPTION
below; validate/inspect_hyper.py shows what the published files actually do.
"""
import argparse
import json
import os
import sys
import time

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402
import matching  # noqa: E402

JET_BRANCHES = ["Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass", "Jet.BTag"]
PART_BRANCHES = ["Particle.PID", "Particle.M1", "Particle.M2",
                 "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass"]
WEIGHT_BRANCH = "Event.Weight"


def find_decay(pid, m1, m2, orig):
    """Locate t -> b W(-> q qbar) for both tops in one event's generator record.

    The arguments are parallel arrays for the quarks, tops and W bosons of the
    event; `orig` holds each particle's index in the full record, which is what
    the mother links m1/m2 refer to.

    Returns {+1: (t_first, b, W, q, qbar), -1: (...)} as positions into these
    arrays, or None if the record does not have the expected structure.

    The decay is read off the last copy of the top (the shower makes several) and
    the last copy of the W, so b, q and qbar are the quarks produced at the decay
    vertices, before they radiate. q is the quark and qbar the antiquark of the W
    decay: HyPER's files put the quark first (labels 2 and 5). t_first is the first
    copy of the top, i.e. the matrix-element top, which is what HyPER stores.
    """
    def children(k):
        return np.flatnonzero((m1 == orig[k]) | (m2 == orig[k]))

    def first_copy(k):
        for _ in range(100):
            parent = np.flatnonzero(orig == m1[k])
            if parent.size == 0 or pid[parent[0]] != pid[k]:
                break
            k = parent[0]
        return k

    out = {}
    for sign in (+1, -1):
        found = None
        for k in np.flatnonzero(pid == 6 * sign):          # ascending record order
            kids = children(k)
            if kids.size == 0 or np.any(pid[kids] == pid[k]):
                continue                                  # not the last copy
            w = kids[np.abs(pid[kids]) == 24]
            b = kids[np.abs(pid[kids]) == 5]
            if w.size != 1 or b.size != 1:
                continue
            w_first = w = w[0]
            for _ in range(100):                          # follow the W to its last copy
                wk = children(w)
                same = wk[pid[wk] == pid[w]]
                if same.size == 0:
                    break
                w = same[0]
            wk = children(w)
            q = wk[(np.abs(pid[wk]) >= 1) & (np.abs(pid[wk]) <= 5)]
            if q.size != 2:
                continue
            q, qbar = (q[0], q[1]) if pid[q[0]] > 0 or pid[q[1]] < 0 else (q[1], q[0])
            # if the record holds the decay twice, the later one is the showered one
            found = (first_copy(k), b[0], w_first, q, qbar)
        if found is None:
            return None
        out[sign] = found
    return out


def convert_chunk(arr, opts, rng, counts):
    n = len(arr["Jet.PT"])
    counts["n_read"] += n
    weight = np.asarray(ak.fill_none(ak.firsts(arr[WEIGHT_BRANCH]), 1.0)) if WEIGHT_BRANCH in arr else np.ones(n)
    counts["n_negative_weight"] += int(np.sum(weight < 0))

    # ---- jets ----
    keep = (arr["Jet.PT"] >= opts.jet_pt_min) & (abs(arr["Jet.Eta"]) < opts.jet_abs_eta_max)
    order = ak.argsort(arr["Jet.PT"][keep], axis=1, ascending=False)
    jet = {k: arr[k][keep][order] for k in JET_BRANCHES}
    btag = (ak.values_astype(jet["Jet.BTag"], np.int64) & 1)
    njet = np.asarray(ak.num(jet["Jet.PT"]))
    nb = np.asarray(ak.sum(btag, axis=1))
    pass_njet = njet >= 6
    pass_b = pass_njet & (nb >= 2)
    counts["n_pass_njet"] += int(pass_njet.sum())
    counts["n_pass_btag"] += int(pass_b.sum())

    # ---- generator record, reduced to quarks, tops and W bosons ----
    apid = abs(arr["Particle.PID"])
    sel = (apid <= 6) | (apid == 24)
    part = {k: arr[k][sel] for k in PART_BRANCHES}
    part_orig = ak.local_index(arr["Particle.PID"], axis=1)[sel]
    offs = np.concatenate([[0], np.cumsum(np.asarray(ak.num(part["Particle.PID"])))])
    flat = {k: np.asarray(ak.flatten(v)) for k, v in part.items()}
    flat_orig = np.asarray(ak.flatten(part_orig))
    joffs = np.concatenate([[0], np.cumsum(njet)])
    jflat = {k: np.asarray(ak.flatten(v)) for k, v in jet.items()}
    jb_flat = np.asarray(ak.flatten(btag))

    out = {k: [] for k in ["jet_pt", "jet_eta", "jet_phi", "jet_e", "jet_m", "jet_bTag",
                           "particle_pid", "particle_pt", "particle_eta", "particle_phi", "particle_m",
                           "quark_pid", "quark_pt", "quark_eta", "quark_phi", "quark_m"]}
    ev = {k: [] for k in ["nbTagged", "eventWeight"]}

    for i in np.flatnonzero(pass_b):
        a, b = offs[i], offs[i + 1]
        pid = flat["Particle.PID"][a:b]
        dec = find_decay(pid, flat["Particle.M1"][a:b], flat["Particle.M2"][a:b], flat_orig[a:b])
        if dec is None:
            counts["n_truth_fail"] += 1
            continue
        p_pt, p_eta, p_phi, p_m = (flat[k][a:b] for k in ["Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass"])
        ja, jb = joffs[i], joffs[i + 1]
        j_pt, j_eta, j_phi, j_m = (jflat[k][ja:jb] for k in ["Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"])

        out["jet_pt"].append(j_pt)
        out["jet_eta"].append(j_eta)
        out["jet_phi"].append(j_phi)
        out["jet_m"].append(j_m)
        out["jet_e"].append(np.sqrt((j_pt * np.cosh(j_eta)) ** 2 + j_m ** 2))
        out["jet_bTag"].append(jb_flat[ja:jb])
        # stored particles, as in HyPER's files: t, tbar, W+, W-
        keep_p = [dec[+1][0], dec[-1][0], dec[+1][2], dec[-1][2]]
        # the six decay quarks: b, q, qbar of the top, then of the antitop
        keep_q = [dec[+1][1], dec[+1][3], dec[+1][4], dec[-1][1], dec[-1][3], dec[-1][4]]
        for name, keep in (("particle", keep_p), ("quark", keep_q)):
            out[name + "_pid"].append(pid[keep])
            out[name + "_pt"].append(p_pt[keep])
            out[name + "_eta"].append(p_eta[keep])
            out[name + "_phi"].append(p_phi[keep])
            out[name + "_m"].append(p_m[keep])
        ev["nbTagged"].append(nb[i])
        ev["eventWeight"].append(weight[i])

    if not ev["nbTagged"]:
        return None
    arrays = {k: ak.Array(v) for k, v in out.items()}
    arrays.update({k: np.asarray(v) for k, v in ev.items()})
    # Labels and the matched-W flag. Every event passing the jet and b-tag cuts is
    # written, flagged; scripts/merge.py applies the matched-W requirement.
    arrays = matching.relabel(arrays, opts.match_dr, opts.conflict, opts.top_order, rng)
    for key, val in matching.summary(arrays).items():
        if key == "n_events":
            counts["n_prematch"] += val
        elif key == "n_label_present":
            counts[key] = [x + y for x, y in zip(counts[key], val)]
        else:
            counts[key] += val
    return arrays


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Delphes ROOT file")
    ap.add_argument("--output", required=True, help="output ROOT file: HyPER's layout plus quark_* branches")
    ap.add_argument("--meta", required=True, help="output JSON with counts for bookkeeping")
    ap.add_argument("--jet-pt-min", type=float, default=25.0)
    ap.add_argument("--jet-abs-eta-max", type=float, default=5.0)
    ap.add_argument("--match-dr", type=float, default=0.4)
    ap.add_argument("--top-order", choices=matching.TOP_ORDERS, default="t-first", help="which top gets labels 1-3")
    ap.add_argument("--conflict", choices=matching.CONFLICTS, default="closest",
                    help="what happens when two quarks are closest to the same jet (see pylib/matching.py)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--step", type=int, default=2000, help="events read per chunk")
    ap.add_argument("--max-truth-fail-frac", type=float, default=0.01,
                    help="exit with an error if more than this fraction of events lacks the t -> b W(-> q q) record")
    ap.add_argument("--extra", action="append", default=[], help="key=value copied into the JSON")
    opts = ap.parse_args()

    t0 = time.time()
    rng = np.random.default_rng(opts.seed)
    counts = {k: 0 for k in ["n_read", "n_negative_weight", "n_pass_njet", "n_pass_btag", "n_truth_fail", "n_prematch",
                             "n_selected", "n_all_matched", "n_top1_matched", "n_top2_matched"]}
    counts["n_label_present"] = [0] * 6

    writer = hyperio.Writer(opts.output)
    with uproot.open(opts.input) as f:
        tree = f["Delphes"]
        # Look each branch up by name (uproot searches sub-branches such as Jet/Jet.PT)
        # instead of passing the names as expressions, where "Jet.PT" could be parsed.
        branches = {name: tree[name] for name in JET_BRANCHES + PART_BRANCHES}
        try:
            branches[WEIGHT_BRANCH] = tree[WEIGHT_BRANCH]
        except uproot.KeyInFileError:
            print("convert: no %s branch; event weights set to 1" % WEIGHT_BRANCH)
        for start in range(0, tree.num_entries, opts.step):
            stop = min(start + opts.step, tree.num_entries)
            arr = {name: br.array(entry_start=start, entry_stop=stop) for name, br in branches.items()}
            chunk = convert_chunk(arr, opts, rng, counts)
            if chunk is not None:
                writer.write(chunk)
    writer.close()
    if counts["n_prematch"] == 0:
        os.remove(opts.output)        # an empty tree has no types; the JSON records the zero

    meta = dict(counts)
    meta.update({"seed": opts.seed, "jet_pt_min": opts.jet_pt_min, "jet_abs_eta_max": opts.jet_abs_eta_max,
                 "match_dr": opts.match_dr, "top_order": opts.top_order, "conflict": opts.conflict,
                 "seconds_convert": round(time.time() - t0, 1),
                 "uproot": uproot.__version__, "awkward": ak.__version__})
    for item in opts.extra:
        key, _, value = item.partition("=")
        try:
            meta[key] = int(value)
        except ValueError:
            meta[key] = value
    with open(opts.meta, "w") as f:
        json.dump(meta, f, indent=1, sort_keys=True)

    n = max(counts["n_read"], 1)
    print("convert: read %d, >=6 jets %d, >=2 b-tags %d, written %d, with a matched W %d (%.2f%%), all matched %d, "
          "truth failures %d" % (counts["n_read"], counts["n_pass_njet"], counts["n_pass_btag"], counts["n_prematch"],
                                 counts["n_selected"], 100.0 * counts["n_selected"] / n, counts["n_all_matched"],
                                 counts["n_truth_fail"]))
    if counts["n_pass_btag"] and counts["n_truth_fail"] > opts.max_truth_fail_frac * counts["n_pass_btag"]:
        sys.exit("convert: %d of %d events lacked the expected t -> b W(-> q q) generator record; "
                 "the sample may not be all-hadronic, or the record layout differs from what find_decay expects"
                 % (counts["n_truth_fail"], counts["n_pass_btag"]))


if __name__ == "__main__":
    main()
