#!/bin/bash
# Put back the cards that MG5 removed from the process directory during the grid-setup
# run, without redoing the integration. The originals are in gridpack_fresh.tar.gz.
# Works on node-local scratch; run it in an interactive session.
#
# usage: scripts/repair_gridpack.sh
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

FIX=$(mktemp -d "${SCRATCH_BASE}/ttgen_repair_XXXXXX")
trap 'rm -rf "${FIX}"' EXIT
cd "${FIX}"
mkdir fresh
tar -xzf "${WORK}/gridpack_fresh.tar.gz" -C fresh "${PROC_NAME}/Cards"
tar -xzf "${WORK}/gridpack_reuse.tar.gz"

echo "repair: cards missing from the reuse gridpack:"
comm -23 <(cd "fresh/${PROC_NAME}/Cards" && find . -type f | sort) \
         <(cd "${PROC_NAME}/Cards" && find . -type f | sort) | sed 's/^/  /'
echo "repair: cards that differ:"
for f in $(cd "fresh/${PROC_NAME}/Cards" && find . -type f | sort); do
    if [[ -f "${PROC_NAME}/Cards/${f}" ]] && ! cmp -s "fresh/${PROC_NAME}/Cards/${f}" "${PROC_NAME}/Cards/${f}"; then
        echo "  ${f}"
    fi
done

cp -a "fresh/${PROC_NAME}/Cards/." "${PROC_NAME}/Cards/"
for f in run_card.dat shower_card.dat madspin_card.dat param_card.dat; do
    [[ -s "${PROC_NAME}/Cards/${f}" ]] || { echo "repair: ${f} is still missing" >&2; exit 1; }
done
tar -czf gridpack_reuse.tar.gz "${PROC_NAME}"
cp gridpack_reuse.tar.gz "${WORK}/gridpack_reuse.tar.gz.new"
mv "${WORK}/gridpack_reuse.tar.gz.new" "${WORK}/gridpack_reuse.tar.gz"
echo "repair: wrote ${WORK}/gridpack_reuse.tar.gz"
