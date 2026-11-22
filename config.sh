#!/bin/bash
# Site and production settings. Every script sources this file.
# Lines marked CHECK are guesses about HPC3 that you must confirm before running.

# ---- paths -----------------------------------------------------------------
TTGEN_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TTGEN_ROOT
export WORK="/pub/${USER}/ttbar-gen/work"        # CHECK: persistent area for the container, gridpack, logs
export OUTDIR="/pub/${USER}/ttbar-gen/out"       # CHECK: where converted files land (about 20 GB for 10x)
export SCRATCH_BASE="${TMPDIR:-/tmp}"            # CHECK: node-local scratch; needs ~15 GB free per running job
export SIF="${WORK}/madlad-hyper.sif"

# ---- cluster ---------------------------------------------------------------
export SLURM_ACCOUNT="asterenb_lab"                             # CHECK: required, your HPC3 account
export SLURM_PARTITION="standard"                # CHECK
export CONTAINER_CMD="singularity"               # MadLAD itself looks for a binary called "singularity"
export MODULE_LOADS="singularity/3.11.3"                           # CHECK: e.g. "singularity/3.11.3"; space-separated, may be empty
export CORES_PER_JOB=2
export MEM_PER_CPU="3G"
export WALLTIME="06:00:00"                       # revise after the pilot
export ARRAY_CHUNK=1000                          # tasks per sbatch call (stay under the site's MaxArraySize)
export ARRAY_THROTTLE=200                        # simultaneous tasks per chunk

# ---- production ------------------------------------------------------------
# Target: 10x HyPER's 7,921,299 selected training events.
# HyPER selected 8,801,060 of 60,000,000 generated (14.67%), so 10x needs ~540M generated.
export N_GEN_TOTAL=540000000
export EVENTS_PER_JOB=25000
export SEED_BASE=10000                           # job seed = SEED_BASE + task id; use a new base for a new campaign
export PROC_NAME="ttbar_allhad"
export MADLAD_REPO="https://github.com/tzuhanchang/MadLAD.git"
export MADLAD_COMMIT="4b0983c3fda9d95c6aaed17d58224d9c8073456b"   # the commit this pipeline was written against

# reuse: integrate the NLO grids once in setup, jobs only generate events (fast; untested here).
# fresh: every job compiles and integrates from scratch, exactly as MadLAD does (slow; safest).
export GRID_MODE="reuse"
export GRID_SETUP_EVENTS=10000

# ---- conversion (assumptions about HyPER's conventions; see README) ----------
export JET_PT_MIN=25.0                           # no-op: the Delphes card already cuts at 25 GeV before the energy scale
export JET_ABS_ETA_MAX=5.0                       # HyPER's files keep jets to |eta| = 4.98: no cut beyond Delphes acceptance
export MATCH_DR=0.4
export TOP_ORDER="t-first"                       # t-first | tbar-first | pt | random
export MATCH_CONFLICT="closest"                   # closest | first-wins | last-wins | greedy; can be changed later with merge.py --rematch

# ---- helpers ---------------------------------------------------------------
ttgen_modules() {
    local m
    for m in ${MODULE_LOADS}; do module load "$m"; done
}

# Run a command inside the container with the given directory as cwd and home.
# usage: cexec <dir> <cmd> [args...]
cexec() {
    local dir="$1"; shift
    "${CONTAINER_CMD}" exec --cleanenv --home "${dir}" \
        --bind "${dir}" --bind "${TTGEN_ROOT}" --bind "${WORK}" --pwd "${dir}" \
        "${SIF}" "$@"
}

ttgen_njobs() {
    echo $(( (N_GEN_TOTAL + EVENTS_PER_JOB - 1) / EVENTS_PER_JOB ))
}
