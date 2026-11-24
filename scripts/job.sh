#!/bin/bash
# One production task: matrix element + MadSpin + Pythia (MG5_aMC), Delphes, conversion.
# Runs in node-local scratch and copies two small files back: <tag>.root and <tag>.json.
#
# As an array task:  sbatch --array=... scripts/job.sh        (task = OFFSET + SLURM_ARRAY_TASK_ID)
# By hand (pilot):   scripts/job.sh <task id> [events]
set -euo pipefail
source "${TTGEN_ROOT:-$(dirname "${BASH_SOURCE[0]}")/..}/config.sh"
ttgen_modules

if [[ $# -ge 1 ]]; then
    TASK="$1"
else
    TASK=$(( ${OFFSET:-0} + ${SLURM_ARRAY_TASK_ID:?give a task id or run as an array job} ))
fi
NEVENTS="${2:-${EVENTS_PER_JOB}}"
NCORES="${SLURM_CPUS_PER_TASK:-${CORES_PER_JOB}}"
SEED=$(( SEED_BASE + TASK ))
TAG=$(printf "%s_%06d" "${PROC_NAME}" "${TASK}")
# Output goes into subfolders of RAW_BUCKET tasks (raw/000, raw/001, ...), so no single
# directory on the shared filesystem collects tens of thousands of files.
RAWDIR="${OUTDIR}/raw/$(printf '%03d' $(( (TASK - 1) / ${RAW_BUCKET:-1000} )))"
FINAL_ROOT="${RAWDIR}/${TAG}.root"
FINAL_META="${RAWDIR}/${TAG}.json"

if [[ -s "${FINAL_META}" ]]; then
    echo "job ${TAG}: already done"
    exit 0
fi
mkdir -p "${RAWDIR}"
echo "job ${TAG}: seed ${SEED}, ${NEVENTS} events, ${NCORES} cores, host $(hostname)"
T0=${SECONDS}

JOBDIR=$(mktemp -d "${SCRATCH_BASE}/ttgen_${TAG}_XXXXXX")
trap 'rm -rf "${JOBDIR}"' EXIT
cd "${JOBDIR}"

# ---- stage in ----------------------------------------------------------------
# Everything the task needs is copied to node-local disk first, so the shared
# filesystem is touched only here and in the final copy-back. The container image and
# the gridpack go to a per-node cache shared by this user's tasks on the node.
NODE_CACHE="${NODE_CACHE_DIR:-/tmp/${USER}/ttgen_cache}"
mkdir -p "${NODE_CACHE}" 2>/dev/null || NODE_CACHE="${JOBDIR}/cache"
mkdir -p "${NODE_CACHE}"

# usage: stage <file on shared filesystem>  -> prints the local path
stage() {
    local src="$1" dst
    dst="${NODE_CACHE}/$(basename "$1")"
    (
        flock 9
        # same size and modification time as the source means the cached copy is current
        if [[ ! -s "${dst}" || "$(stat -c '%s %Y' "${src}")" != "$(stat -c '%s %Y' "${dst}")" ]]; then
            cp -p "${src}" "${dst}.tmp.$$" && mv "${dst}.tmp.$$" "${dst}"
        fi
    ) 9> "${dst}.lock"
    [[ -s "${dst}" ]] || { echo "job ${TAG}: could not stage ${src}" >&2; return 1; }
    echo "${dst}"
}
[[ -s "${WORK}/gridpack_${GRID_MODE}.tar.gz" && -s "${WORK}/ATLAS_sim_FASTJET_slim.tcl" && -s "${SIF}" ]] \
    || { echo "job ${TAG}: run scripts/00_setup.sh first" >&2; exit 1; }
GRIDPACK=$(stage "${WORK}/gridpack_${GRID_MODE}.tar.gz")
CARD_SRC=$(stage "${WORK}/ATLAS_sim_FASTJET_slim.tcl")
LOCAL_SIF=$(stage "${SIF}")
mkdir -p code
cp -r "${TTGEN_ROOT}/scripts" "${TTGEN_ROOT}/convert" "${TTGEN_ROOT}/pylib" code/
T_STAGE=$(( SECONDS - T0 ))
echo "job ${TAG}: staged in ${T_STAGE} s, dir ${JOBDIR}, cache ${NODE_CACHE}"

# Run a command in the local copy of the container, seeing only the job directory.
jexec() {
    "${CONTAINER_CMD}" exec --cleanenv --home "${JOBDIR}" --bind "${JOBDIR}" --pwd "${JOBDIR}" "${LOCAL_SIF}" "$@"
}

# ---- cards -------------------------------------------------------------------
tar -xzf "${GRIDPACK}"
python3 code/scripts/set_card.py run "${PROC_NAME}/Cards/run_card.dat" \
    "nevents=${NEVENTS}" "iseed=${SEED}"
# MadLAD leaves the shower seed at its default, which gives every job the same
# Pythia random-number stream. Give each job its own.
python3 code/scripts/set_card.py shower "${PROC_NAME}/Cards/shower_card.dat" "rnd_seed=${SEED}"

if [[ "${GRID_MODE}" == "reuse" ]]; then
    # -o: skip the grid setup and generate from the stored grids; -n: do not recompile.
    printf 'set run_mode 2\nset nb_core %s\nlaunch %s -i\ngenerate_events -o -n -f\n' "${NCORES}" "${PROC_NAME}" > mg5_exec_card
else
    # exactly what MadLAD writes for an NLO run with shower
    printf 'set run_mode 2\nset nb_core %s\nlaunch %s\n' "${NCORES}" "${PROC_NAME}" > mg5_exec_card
fi

# ---- generation ----------------------------------------------------------------
TG=${SECONDS}
jexec /app/MG5_aMC/bin/mg5_aMC mg5_exec_card
T_GEN=$(( SECONDS - TG ))

# The showered file must come from the MadSpin-decayed run. A HepMC file anywhere
# else means the tops were not decayed by MadSpin and the sample is not all-hadronic.
mapfile -t HEPMC < <(find "${PROC_NAME}/Events" -path '*_decayed_*' -name '*.hepmc.gz')
if [[ ${#HEPMC[@]} -ne 1 ]]; then
    echo "job ${TAG}: expected one showered file under Events/*_decayed_*, found ${#HEPMC[@]}" >&2
    find "${PROC_NAME}/Events" -maxdepth 2 >&2 || true
    exit 1
fi
HEPMC_BYTES=$(stat -c %s "${HEPMC[0]}")
SCRATCH_PEAK_KB=$(du -sk "${JOBDIR}" | cut -f1)

# ---- detector simulation ---------------------------------------------------------
{ echo "set RandomSeed ${SEED}"; cat "${CARD_SRC}"; } > delphes_card.tcl
T1=${SECONDS}
# shellcheck disable=SC2016
jexec /bin/bash -c 'set -o pipefail; gunzip -c "$1" | DelphesHepMC2 delphes_card.tcl delphes.root' _ "${HEPMC[0]}"
T_DELPHES=$(( SECONDS - T1 ))
rm -rf "${PROC_NAME}"
[[ -s delphes.root ]] || { echo "job ${TAG}: Delphes wrote no output" >&2; exit 1; }

# ---- conversion: HyPER's layout plus the decay quarks ---------------------------------------------------
T2=${SECONDS}
jexec python3 code/convert/delphes_to_hyper.py \
    --input delphes.root --output "${TAG}.root" --meta "${TAG}.json" \
    --jet-pt-min "${JET_PT_MIN}" --jet-abs-eta-max "${JET_ABS_ETA_MAX}" \
    --match-dr "${MATCH_DR}" --top-order "${TOP_ORDER}" --conflict "${MATCH_CONFLICT}" --seed "${SEED}" \
    --extra "task=${TASK}" --extra "n_requested=${NEVENTS}" --extra "grid_mode=${GRID_MODE}" \
    --extra "seconds_stage_in=${T_STAGE}" --extra "seconds_generation=${T_GEN}" --extra "seconds_delphes=${T_DELPHES}" \
    --extra "hepmc_gz_bytes=${HEPMC_BYTES}" --extra "scratch_after_generation_kb=${SCRATCH_PEAK_KB}" \
    --extra "delphes_root_bytes=$(stat -c %s delphes.root)" --extra "cores=${NCORES}"
T_CONVERT=$(( SECONDS - T2 ))

# ---- copy back; the json goes last and marks the task as done -----------------------
if [[ -s "${TAG}.root" ]]; then
    cp "${TAG}.root" "${RAWDIR}/.${TAG}.root.tmp"
    mv "${RAWDIR}/.${TAG}.root.tmp" "${FINAL_ROOT}"
fi
cp "${TAG}.json" "${RAWDIR}/.${TAG}.json.tmp"
mv "${RAWDIR}/.${TAG}.json.tmp" "${FINAL_META}"
echo "job ${TAG}: done in $(( SECONDS - T0 )) s (stage-in ${T_STAGE}, generation ${T_GEN}, delphes ${T_DELPHES}, convert ${T_CONVERT})"
