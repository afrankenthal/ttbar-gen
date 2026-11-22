#!/bin/bash
# Submit production tasks as slurm arrays, ARRAY_CHUNK tasks per sbatch call.
#
# usage: scripts/submit.sh                 all tasks 1..N (N from N_GEN_TOTAL / EVENTS_PER_JOB)
#        scripts/submit.sh <first> <last>  a range, e.g. 1 200 for a first tranche
#        scripts/submit.sh --list <file>   task ids from a file (written by status.py --missing)
# Finished tasks exit immediately, so resubmitting a range is safe.
# Add DRYRUN=1 to print the sbatch commands without submitting.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
[[ -n "${SLURM_ACCOUNT}" ]] || { echo "submit: set SLURM_ACCOUNT in config.sh" >&2; exit 1; }
mkdir -p "${WORK}/logs" "${OUTDIR}/raw"

SBATCH=(sbatch --parsable --job-name=ttgen --account="${SLURM_ACCOUNT}" --partition="${SLURM_PARTITION}"
        --nodes=1 --ntasks=1 --cpus-per-task="${CORES_PER_JOB}" --mem-per-cpu="${MEM_PER_CPU}"
        --time="${WALLTIME}" --output="${WORK}/logs/%A_%a.out")
run() { if [[ -n "${DRYRUN:-}" ]]; then echo "$*"; else "$@"; fi; }

if [[ "${1:-}" == "--list" ]]; then
    # arbitrary ids: one array per chunk-sized block of the 1..N range they fall in
    mapfile -t IDS < <(sort -n "$2" | uniq)
    declare -A BLOCKS=()
    for id in "${IDS[@]}"; do
        block=$(( (id - 1) / ARRAY_CHUNK ))
        local_id=$(( id - block * ARRAY_CHUNK ))
        BLOCKS[${block}]="${BLOCKS[${block}]:+${BLOCKS[${block}]},}${local_id}"
    done
    for block in "${!BLOCKS[@]}"; do
        run "${SBATCH[@]}" --array="${BLOCKS[${block}]}%${ARRAY_THROTTLE}" \
            --export="ALL,TTGEN_ROOT=${TTGEN_ROOT},OFFSET=$(( block * ARRAY_CHUNK ))" "${TTGEN_ROOT}/scripts/job.sh"
    done
    exit 0
fi

FIRST="${1:-1}"
LAST="${2:-$(ttgen_njobs)}"
echo "submit: tasks ${FIRST}..${LAST} of $(ttgen_njobs), ${EVENTS_PER_JOB} events each"
block_first=$(( (FIRST - 1) / ARRAY_CHUNK ))
block_last=$(( (LAST - 1) / ARRAY_CHUNK ))
for (( block = block_first; block <= block_last; block++ )); do
    offset=$(( block * ARRAY_CHUNK ))
    lo=$(( FIRST > offset + 1 ? FIRST - offset : 1 ))
    hi=$(( LAST < offset + ARRAY_CHUNK ? LAST - offset : ARRAY_CHUNK ))
    run "${SBATCH[@]}" --array="${lo}-${hi}%${ARRAY_THROTTLE}" \
        --export="ALL,TTGEN_ROOT=${TTGEN_ROOT},OFFSET=${offset}" "${TTGEN_ROOT}/scripts/job.sh"
done
