"""Truth matching of jets to the six decay quarks, applied to whole arrays of events.

The same code labels events at conversion time and re-labels them later from the
stored quarks (scripts/merge.py --rematch, validate/rematch_scan.py), so a different
matching rule never needs new events.

Quarks are stored per event as [b, q, qbar] of the top, then [b, q, qbar] of the
antitop. Labels: 1 = b of top 1, 2 = quark and 3 = antiquark of its W, 4-6 the same
for top 2. Which physical top is "top 1" is the top_order setting.

Rules for two quarks whose closest jet is the same jet ("conflict"):
  closest      the closer quark keeps the jet, the other is unmatched
  first-wins   the quark with the lower label keeps it
  last-wins    the quark with the higher label keeps it (a loop that overwrites)
  greedy       pairs are assigned in order of increasing dR, so the loser can still
               take its next-closest free jet
"""
import awkward as ak
import numpy as np

CONFLICTS = ["closest", "first-wins", "last-wins", "greedy"]
TOP_ORDERS = ["t-first", "tbar-first", "pt", "random"]
QUARK_FIELDS = ["pid", "pt", "eta", "phi", "m"]


def labels_for(q_eta, q_phi, j_eta, j_phi, dr_max, conflict):
    """q_*: (N, 6) in label order; j_*: (N, J), NaN where an event has fewer jets.
    Returns integer labels (N, J)."""
    n, nj = j_eta.shape
    rows = np.arange(n)
    with np.errstate(invalid="ignore"):
        dphi = np.abs(q_phi[:, :, None] - j_phi[:, None, :])
        dphi = np.where(dphi > np.pi, 2 * np.pi - dphi, dphi)
        dr = np.hypot(q_eta[:, :, None] - j_eta[:, None, :], dphi)
    dr = np.where(np.isnan(dr), np.inf, dr)
    labels = np.zeros((n, nj), dtype=np.int32)

    if conflict == "greedy":
        d = dr.copy()
        for _ in range(6):
            flat = d.reshape(n, -1)
            k = np.argmin(flat, axis=1)
            ok = flat[rows, k] < dr_max
            qi, ji = k // nj, k % nj
            labels[rows[ok], ji[ok]] = qi[ok] + 1
            d[rows[ok], qi[ok], :] = np.inf
            d[rows[ok], :, ji[ok]] = np.inf
        return labels

    best = np.argmin(dr, axis=2)
    best_dr = np.take_along_axis(dr, best[:, :, None], axis=2)[:, :, 0]
    valid = best_dr < dr_max
    same = (best[:, :, None] == best[:, None, :]) & valid[:, :, None] & valid[:, None, :]
    i = np.arange(6)[None, :, None]
    k = np.arange(6)[None, None, :]
    if conflict == "closest":
        di, dk = best_dr[:, :, None], best_dr[:, None, :]
        beaten_by = (dk < di) | ((dk == di) & (k < i))
    elif conflict == "first-wins":
        beaten_by = np.broadcast_to(k < i, same.shape)
    elif conflict == "last-wins":
        beaten_by = np.broadcast_to(k > i, same.shape)
    else:
        raise ValueError("unknown conflict rule %r" % conflict)
    win = valid & ~np.any(same & beaten_by & (i != k), axis=2)
    for q in range(6):
        w = win[:, q]
        labels[rows[w], best[w, q]] = q + 1
    return labels


def relabel(arrays, dr_max=0.4, conflict="closest", top_order="t-first", rng=None):
    """Return a copy of `arrays` with jet_truthmatch, allMatchedEvent and leastOneMatchedW
    recomputed from the stored quarks."""
    n = len(arrays["jet_pt"])
    out = dict(arrays)
    if n == 0:
        return out
    njet = np.asarray(ak.num(arrays["jet_pt"]))
    nj = int(njet.max())

    def pad(x):
        return np.asarray(ak.fill_none(ak.pad_none(ak.values_astype(x, np.float64), nj, axis=1, clip=True), np.nan))

    q_eta = np.asarray(ak.to_numpy(ak.values_astype(arrays["quark_eta"], np.float64)))
    q_phi = np.asarray(ak.to_numpy(ak.values_astype(arrays["quark_phi"], np.float64)))
    if q_eta.shape[1] != 6:
        raise ValueError("expected six stored quarks per event")

    if top_order == "t-first":
        swap = np.zeros(n, dtype=bool)
    elif top_order == "tbar-first":
        swap = np.ones(n, dtype=bool)
    elif top_order == "pt":
        tops = np.asarray(ak.to_numpy(arrays["particle_pt"][:, :2]))      # t, tbar
        swap = tops[:, 1] > tops[:, 0]
    elif top_order == "random":
        swap = (rng or np.random.default_rng(0)).random(n) < 0.5
    else:
        raise ValueError("unknown top order %r" % top_order)
    idx = np.where(swap[:, None], np.array([3, 4, 5, 0, 1, 2])[None, :], np.arange(6)[None, :])
    q_eta = np.take_along_axis(q_eta, idx, axis=1)
    q_phi = np.take_along_axis(q_phi, idx, axis=1)

    labels = labels_for(q_eta, q_phi, pad(arrays["jet_eta"]), pad(arrays["jet_phi"]), dr_max, conflict)
    real = np.arange(nj)[None, :] < njet[:, None]
    out["jet_truthmatch"] = ak.unflatten(labels[real], njet)
    has = np.stack([np.any(labels == k, axis=1) for k in range(1, 7)], axis=1)
    out["allMatchedEvent"] = has.all(axis=1).astype(np.int32)
    out["leastOneMatchedW"] = ((has[:, 1] & has[:, 2]) | (has[:, 4] & has[:, 5])).astype(np.int32)
    return out


def summary(arrays):
    """Counts that describe a labelling; rates are relative to events with a matched W."""
    lab = arrays["jet_truthmatch"]
    has = np.stack([np.asarray(ak.any(lab == k, axis=1)) for k in range(1, 7)], axis=1)
    sel = np.asarray(arrays["leastOneMatchedW"]) > 0
    return {
        "n_events": int(len(sel)),
        "n_selected": int(sel.sum()),
        "n_all_matched": int((has.all(axis=1) & sel).sum()),
        "n_top1_matched": int((has[:, :3].all(axis=1) & sel).sum()),
        "n_top2_matched": int((has[:, 3:].all(axis=1) & sel).sum()),
        "n_label_present": [int((has[:, k] & sel).sum()) for k in range(6)],
    }


def select(arrays, mask):
    return {k: v[mask] for k, v in arrays.items()}
