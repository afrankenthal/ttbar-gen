#!/usr/bin/env python3
"""Per-task files -> status, merge, re-match. Uses the synthetic Delphes file that
tests/test_convert.py leaves behind.   Run:  python3 tests/test_merge.py <delphes.root>
"""
import glob
import json
import os
import subprocess
import sys
import tempfile

import awkward as ak
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "pylib"))
import hyperio  # noqa: E402
import matching  # noqa: E402


def run(*cmd):
    return subprocess.run([sys.executable] + list(cmd), check=True, stdout=subprocess.PIPE, universal_newlines=True).stdout


def main(delphes):
    tmp = tempfile.mkdtemp(prefix="ttgen_merge_")
    raw = os.path.join(tmp, "raw")
    os.makedirs(raw)
    for task in (1, 2, 3, 5):
        run(os.path.join(ROOT, "convert", "delphes_to_hyper.py"), "--input", delphes,
            "--output", os.path.join(raw, "ttbar_allhad_%06d.root" % task),
            "--meta", os.path.join(raw, "ttbar_allhad_%06d.json" % task), "--max-truth-fail-frac", "0.5",
            "--extra", "task=%d" % task, "--extra", "n_requested=900", "--extra", "seconds_generation=1800",
            "--extra", "seconds_delphes=900", "--extra", "cores=2", "--extra", "scratch_after_generation_kb=9000000")
    status = run(os.path.join(ROOT, "scripts", "status.py"), raw, "--tasks", "6", "--missing", os.path.join(tmp, "missing.txt"))
    print(status)
    assert open(os.path.join(tmp, "missing.txt")).read().split() == ["4", "6"]
    assert "events selected       2000" in status and "stored for re-matching 2400" in status

    rawdata = hyperio.concatenate([hyperio.read(f) for f in sorted(glob.glob(os.path.join(raw, "*.root")))])
    keep = np.asarray(rawdata["leastOneMatchedW"]) > 0

    out = os.path.join(tmp, "merged")
    run(os.path.join(ROOT, "scripts", "merge.py"), raw, out, "--shard-size", "700", "--holdout", "450", "--buffer", "300")
    files = [os.path.join(out, "ttbar_new_holdout.root")] + sorted(glob.glob(os.path.join(out, "ttbar_new_train_*.root")))
    merged = hyperio.concatenate([hyperio.read(f) for f in files])
    assert [hyperio.n_entries(f) for f in files] == [450, 700, 700, 150]
    with uproot.open(files[1]) as f:
        assert f.classname_of("Delphes") == "TTree"
        assert not [k for k in f["Delphes"].keys() if "quark" in k], "merged files must not carry the quark branches"
    expect = hyperio.drop_quarks(matching.select(rawdata, keep))
    assert set(merged) == set(expect)
    for k in expect:
        assert ak.to_list(merged[k]) == ak.to_list(expect[k]), k
    assert np.all(np.asarray(merged["leastOneMatchedW"]) == 1)

    out2 = os.path.join(tmp, "merged_lastwins")
    run(os.path.join(ROOT, "scripts", "merge.py"), raw, out2, "--rematch", "--conflict", "last-wins", "--keep-quarks")
    re = hyperio.read(os.path.join(out2, "ttbar_new_train_000.root"))
    ref = matching.select(matching.relabel(rawdata, 0.4, "last-wins", "t-first"), keep)
    assert ak.to_list(re["jet_truthmatch"]) == ak.to_list(ref["jet_truthmatch"]) and "quark_eta" in re
    assert ak.to_list(re["jet_truthmatch"]) != ak.to_list(expect["jet_truthmatch"])
    man = json.load(open(os.path.join(out2, "ttbar_new_manifest.json")))
    assert man["conflict"] == "last-wins" and man["events"] == 2000 and man["events_before_matched_w_requirement"] == 2400

    print(run(os.path.join(ROOT, "validate", "rematch_scan.py"), os.path.join(raw, "*.root")))
    print("test_merge: OK (files in %s)" % tmp)


if __name__ == "__main__":
    main(sys.argv[1])
