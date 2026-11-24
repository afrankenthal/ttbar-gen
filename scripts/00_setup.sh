#!/bin/bash
# One-time setup: fetch MadLAD, create the MG5 process directory with MadLAD's own
# card-writing code, integrate the NLO grids once, and pack everything for the jobs.
#
# Run this on a compute node (interactive session with ~8 cores), not a login node:
# it compiles the process and, in GRID_MODE=reuse, runs the NLO integration.
#
# usage: scripts/00_setup.sh [ncores]      build (compute node)
#        scripts/00_setup.sh --fetch       download step only (login node)
#
# On clusters whose compute nodes have no internet access, run --fetch on a login node
# first: it downloads MadLAD and the small Python environment into $WORK. The build then
# uses those copies and needs no network. Where compute nodes do have internet, the
# build downloads what it needs itself and --fetch is unnecessary.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
ttgen_modules
mkdir -p "${WORK}" "${OUTDIR}/raw" "${WORK}/logs"

if [[ "${1:-}" == "--fetch" ]]; then
    if [[ ! -d "${WORK}/src/MadLAD/.git" ]]; then
        git clone "${MADLAD_REPO}" "${WORK}/src/MadLAD"
    fi
    git -C "${WORK}/src/MadLAD" checkout --quiet "${MADLAD_COMMIT}"
    if [[ ! -x "${WORK}/venv/bin/python" ]]; then
        python3 -m venv "${WORK}/venv"
        "${WORK}/venv/bin/pip" install --quiet hydra-core omegaconf pyyaml
    fi
    [[ -s "${SIF}" ]] || echo "fetch: reminder: pull the container image to ${SIF} (also from a login node)"
    echo "fetch: done. Next, on a compute node: scripts/00_setup.sh 8"
    exit 0
fi
NCORES="${1:-8}"

command -v "${CONTAINER_CMD}" >/dev/null || { echo "setup: '${CONTAINER_CMD}' not found; load its module or set MODULE_LOADS" >&2; exit 1; }
[[ -s "${SIF}" ]] || { echo "setup: container image ${SIF} not found; see the README, step 3" >&2; exit 1; }

# Build on node-local scratch, not on the shared filesystem: compiling and integrating
# the process touches tens of thousands of small files, which is slow there and loads
# it for everyone. Only the finished tarballs and the Delphes card are copied to WORK.
# The build directory is deleted with the job's scratch; set SETUP_BUILD_DIR to keep it.
BUILD="$(readlink -m "${SETUP_BUILD_DIR:-${SCRATCH_BASE}/ttgen_setup_${USER}}")"
mkdir -p "${BUILD}"
echo "setup: building in ${BUILD}"

# MadLAD calls a binary named "singularity"; where only apptainer exists, give it that name.
if ! command -v singularity >/dev/null; then
    mkdir -p "${BUILD}/bin"
    ln -sf "$(command -v "${CONTAINER_CMD}")" "${BUILD}/bin/singularity"
    export PATH="${BUILD}/bin:${PATH}"
    echo "setup: using ${CONTAINER_CMD} under the name singularity"
fi

# ---- 1. MadLAD at the pinned commit -------------------------------------------
if [[ ! -d "${BUILD}/MadLAD/.git" ]]; then
    if [[ -d "${WORK}/src/MadLAD/.git" ]]; then
        cp -a "${WORK}/src/MadLAD" "${BUILD}/MadLAD"          # fetched earlier, no network needed
    else
        git clone "${MADLAD_REPO}" "${BUILD}/MadLAD"
    fi
fi
git -C "${BUILD}/MadLAD" checkout --quiet "${MADLAD_COMMIT}"
cp "${TTGEN_ROOT}/process/ttbar-allhad-hyper.yaml" "${BUILD}/MadLAD/processes/"

# ---- 2. Delphes card: MadLAD's ATLAS card with a slimmed output tree -----------
python3 "${TTGEN_ROOT}/scripts/make_slim_card.py" \
    "${BUILD}/MadLAD/delphes/ATLAS_sim_FASTJET.tcl" "${WORK}/ATLAS_sim_FASTJET_slim.tcl"

# ---- 3. host-side Python for MadLAD's driver (it calls singularity itself) -----
if [[ -x "${WORK}/venv/bin/python" ]]; then
    VENV_PY="${WORK}/venv/bin/python"                         # fetched earlier
else
    python3 -m venv "${BUILD}/venv"
    "${BUILD}/venv/bin/pip" install --quiet hydra-core omegaconf pyyaml
    VENV_PY="${BUILD}/venv/bin/python"
fi

# ---- 4. process directory, cards written by MadLAD ----------------------------
cd "${BUILD}/MadLAD"
if [[ -d "${PROC_NAME}" ]]; then
    echo "setup: ${BUILD}/MadLAD/${PROC_NAME} already exists; remove it to start over" >&2
    exit 1
fi
"${VENV_PY}" -m madlad.generate --config-name=ttbar-allhad-hyper \
    "run.image=${SIF}" "gen.block_model.save_dir=${PROC_NAME}" \
    "gen.block_settings.nb_core=${NCORES}"

# MadLAD only warns when a step fails, so check what it produced.
for f in Cards/run_card.dat Cards/madspin_card.dat Cards/param_card.dat Cards/shower_card.dat SubProcesses/setscales.f; do
    [[ -s "${PROC_NAME}/${f}" ]] || { echo "setup: ${PROC_NAME}/${f} is missing; process creation failed" >&2; exit 1; }
done
grep -q "mt\*\*2 + 0.5\*(pt\*\*2+ptbar\*\*2)" "${PROC_NAME}/SubProcesses/setscales.f" \
    || { echo "setup: the custom scale was not written into setscales.f" >&2; exit 1; }
grep -q "decay t > w+ b, w+ > j j" "${PROC_NAME}/Cards/madspin_card.dat" \
    || { echo "setup: the MadSpin card does not hold the all-hadronic decays" >&2; exit 1; }
# These must be editable per job; fail now rather than in 20,000 jobs.
cp "${PROC_NAME}/Cards/run_card.dat" /tmp/ttgen_run_card_check.$$
python3 "${TTGEN_ROOT}/scripts/set_card.py" run /tmp/ttgen_run_card_check.$$ nevents=1 iseed=1
cp "${PROC_NAME}/Cards/shower_card.dat" /tmp/ttgen_shower_card_check.$$
python3 "${TTGEN_ROOT}/scripts/set_card.py" shower /tmp/ttgen_shower_card_check.$$ rnd_seed=1
rm -f /tmp/ttgen_run_card_check.$$ /tmp/ttgen_shower_card_check.$$

tar -czf "${BUILD}/gridpack_fresh.tar.gz" "${PROC_NAME}"
cp "${BUILD}/gridpack_fresh.tar.gz" "${WORK}/gridpack_fresh.tar.gz"
echo "setup: wrote ${WORK}/gridpack_fresh.tar.gz"

# ---- 5. integrate the grids once (GRID_MODE=reuse) ----------------------------
if [[ "${GRID_MODE}" == "reuse" ]]; then
    python3 "${TTGEN_ROOT}/scripts/set_card.py" run "${PROC_NAME}/Cards/run_card.dat" \
        "nevents=${GRID_SETUP_EVENTS}" "iseed=${SEED_BASE}"
    # MG5 removes some cards from Cards/ during a parton-level run (seen on HPC3: the
    # shower card was gone afterwards). Keep the originals and put them back.
    rm -rf "${BUILD}/Cards.orig"
    cp -a "${PROC_NAME}/Cards" "${BUILD}/Cards.orig"
    cat > "${BUILD}/MadLAD/mg5_setup_card" <<CARD
set run_mode 2
set nb_core ${NCORES}
launch ${PROC_NAME} -i
generate_events -p -f
CARD
    cexec "${BUILD}/MadLAD" /app/MG5_aMC/bin/mg5_aMC mg5_setup_card
    ls "${PROC_NAME}"/Events/*/events.lhe.gz >/dev/null 2>&1 \
        || { echo "setup: the grid-setup run produced no events; see the MG5 log above" >&2; exit 1; }
    grep -i "cross.section" -A 3 "${PROC_NAME}"/Events/run_01/summary.txt 2>/dev/null || true
    rm -rf "${PROC_NAME}/Events"/*
    cp -a "${BUILD}/Cards.orig/." "${PROC_NAME}/Cards/"
    for f in run_card.dat shower_card.dat madspin_card.dat param_card.dat; do
        [[ -s "${PROC_NAME}/Cards/${f}" ]] || { echo "setup: ${f} is missing after the grid-setup run" >&2; exit 1; }
    done
    tar -czf "${BUILD}/gridpack_reuse.tar.gz" "${PROC_NAME}"
    cp "${BUILD}/gridpack_reuse.tar.gz" "${WORK}/gridpack_reuse.tar.gz"
    echo "setup: wrote ${WORK}/gridpack_reuse.tar.gz"
fi

echo "setup: done. Next: scripts/job.sh 1 (pilot), then scripts/submit.sh"
