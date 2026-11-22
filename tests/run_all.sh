#!/bin/bash
# Runs every test that does not need MadGraph, a container or a cluster.
# Needs: python3 with uproot, awkward, numpy, scikit-learn, matplotlib.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_convert.py | tee /tmp/ttgen_test_convert.log
DELPHES="$(sed -n 's/.*files in \(.*\))$/\1/p' /tmp/ttgen_test_convert.log)/delphes.root"
python3 tests/test_merge.py "${DELPHES}"
python3 tests/test_validate.py | tail -4
