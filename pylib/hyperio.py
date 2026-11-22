"""Read and write files in the layout of the HyPER ttbar dataset (Zenodo 10653837).

Tree "Delphes" with
  event level:  njet, nbTagged, allMatchedEvent, leastOneMatchedW
  per jet:      jet_pt, jet_eta, jet_phi, jet_e, jet_m, jet_bTag, jet_truthmatch
  per particle: particle_pid, particle_pt, particle_eta, particle_phi, particle_m
The published files store the jagged branches as std::vector and, per event, four
particles: t, tbar, W+, W-. Files written here hold the same values as arrays with a
counter branch (uproot cannot write std::vector); uproot and awkward read both the
same way. They also carry nparticle (a counter uproot adds) and eventWeight.
"""
import awkward as ak
import numpy as np
import uproot

TREE = "Delphes"
JET_FIELDS = ["pt", "eta", "phi", "e", "m", "bTag", "truthmatch"]
PARTICLE_FIELDS = ["pid", "pt", "eta", "phi", "m"]
EVENT_FIELDS = ["nbTagged", "allMatchedEvent", "leastOneMatchedW"]
EXTRA_EVENT_FIELDS = ["eventWeight"]
# Per-task files also hold the six decay quarks (see pylib/matching.py), so labels can
# be recomputed later. They are dropped when merging into HyPER's layout.
QUARK_FIELDS = ["pid", "pt", "eta", "phi", "m"]

# Output types. HyPER's own types are printed by validate/inspect_hyper.py; change
# these if a loader turns out to depend on them.
DTYPES = {
    "jet_pt": np.float32, "jet_eta": np.float32, "jet_phi": np.float32,
    "jet_e": np.float32, "jet_m": np.float32,
    "jet_bTag": np.int32, "jet_truthmatch": np.int32,
    "particle_pid": np.int32, "particle_pt": np.float32, "particle_eta": np.float32,
    "particle_phi": np.float32, "particle_m": np.float32,
    "quark_pid": np.int32, "quark_pt": np.float32, "quark_eta": np.float32,
    "quark_phi": np.float32, "quark_m": np.float32,
    "nbTagged": np.int32, "allMatchedEvent": np.int32, "leastOneMatchedW": np.int32,
    "eventWeight": np.float32,
}


def read(path, entry_stop=None, branches=None):
    """Return {branch name: array} for one file. Jagged branches come back as awkward arrays."""
    with uproot.open(path) as f:
        tree = f[TREE]
        want = ["jet_" + k for k in JET_FIELDS] + ["particle_" + k for k in PARTICLE_FIELDS]
        want += ["quark_" + k for k in QUARK_FIELDS]
        want += EVENT_FIELDS + EXTRA_EVENT_FIELDS
        if branches is not None:
            want = [b for b in want if b in branches]
        have = [b for b in want if b in tree.keys()]
        return dict(zip(have, [tree[b].array(entry_stop=entry_stop) for b in have]))


def n_entries(path):
    with uproot.open(path) as f:
        return f[TREE].num_entries if TREE in f else 0


def to_tree_dict(arrays):
    """Group flat {branch: array} into what uproot needs to write jet_* with counter njet."""
    out = {}
    jets = {k: ak.values_astype(arrays["jet_" + k], DTYPES["jet_" + k]) for k in JET_FIELDS}
    out["jet"] = ak.zip(jets)
    if all(("particle_" + k) in arrays for k in PARTICLE_FIELDS):
        parts = {k: ak.values_astype(arrays["particle_" + k], DTYPES["particle_" + k]) for k in PARTICLE_FIELDS}
        out["particle"] = ak.zip(parts)
    if all(("quark_" + k) in arrays for k in QUARK_FIELDS):
        quarks = {k: ak.values_astype(arrays["quark_" + k], DTYPES["quark_" + k]) for k in QUARK_FIELDS}
        out["quark"] = ak.zip(quarks)
    for k in EVENT_FIELDS + EXTRA_EVENT_FIELDS:
        if k in arrays:
            out[k] = np.asarray(arrays[k]).astype(DTYPES[k])
    return out


class Writer(object):
    """Append chunks of {branch: array} to one output file."""

    def __init__(self, path):
        self.file = uproot.recreate(path)
        self.started = False
        self.entries = 0

    def write(self, arrays):
        n = len(arrays["jet_pt"])
        if n == 0:
            return
        data = to_tree_dict(arrays)
        if not self.started:
            # mktree, not file[TREE] = data: recent uproot versions turn the latter
            # into an RNTuple, which TTree-based loaders cannot read.
            types = {k: (v.type.content if isinstance(v, ak.Array) else v.dtype) for k, v in data.items()}
            self.file.mktree(TREE, types)
            self.started = True
        self.file[TREE].extend(data)
        self.entries += n

    def close(self):
        self.file.close()


def drop_quarks(arrays):
    return {k: v for k, v in arrays.items() if not k.startswith("quark_")}


def concatenate(chunks):
    keys = chunks[0].keys()
    return {k: ak.concatenate([c[k] for c in chunks]) for k in keys}
