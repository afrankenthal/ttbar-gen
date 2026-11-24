# ttbar-gen: a 10x all-hadronic ttbar sample matched to the HyPER benchmark

Generates all-hadronic ttbar events with MadGraph5_aMC@NLO + MadSpin + Pythia 8 + Delphes
on a slurm cluster and writes them in the layout of the HyPER dataset
(Zenodo 10.5281/zenodo.10653837), so a loader written for `ttbar_train.root` reads them
unchanged.

**Status (23 Nov 2026): the conversion, matching, merging, bookkeeping and validation code
is tested on synthetic files. The container is built (GitHub Actions) and pulled to HPC3,
with Delphes, MadGraph, Pythia and LHAPDF confirmed inside. No MadGraph step and no slurm
job has run yet.** The setup and the first pilot job will need fixes; the section
"Not yet run" lists where to expect them.

## Target

| | |
|---|---|
| HyPER, selected training events | 7,921,299 |
| HyPER, selection efficiency | 8,801,060 of 60,000,000 generated = 14.67% |
| Target, selected events | 79.2M (10x) |
| Events to generate | 540M |
| Tasks | 21,600 x 25,000 events |
| Final size | about 20 GB merged (HyPER's 7.9M-event file is 1.9 GB), plus an estimated 30 GB of per-task files kept for re-matching |

The 540M figure assumes this setup reproduces HyPER's efficiency. The pilot measures it.

## What is taken from where

| Piece | Source |
|---|---|
| Process, PDF, scales, MadSpin decays, top mass | MadLAD `processes/ttbar-allhad.yaml`, commit `4b0983c` |
| Delphes card | MadLAD `delphes/ATLAS_sim_FASTJET.tcl` (anti-kt R=0.4, jet pT > 25 GeV, pT-dependent b-tag efficiency) |
| Software versions | HyPER's Zenodo record: MG5_aMC 2.9.16, Pythia 8.306, Delphes 3.5.0, FastJet 3.4.0, LHAPDF 6.3.0, NNPDF3.0nlo |
| Selection and labels | HyPER's Zenodo record and the paper draft's Methods |

MadLAD is written by HyPER's contact author and its Delphes card matches what HyPER
describes, but **nobody has confirmed that this config produced the HyPER files**. MadLAD's
example container pins newer versions (MG5 2.9.26, Pythia 8.313), so the versions here
come from the record, not from MadLAD.

### Deliberate differences from running MadLAD as-is

1. **Scale and PDF reweighting off.** It only attaches extra weights to each event.
2. **NLO grids integrated once** (`GRID_MODE=reuse`) and reused by every task with
   `generate_events -o -n`. MadLAD integrates from scratch in every job; `GRID_MODE=fresh`
   does that.
3. **Shower and Delphes seeds set per task.** MadLAD changes only the matrix-element seed,
   which leaves every job on the same Pythia random-number stream.
4. **Tasks work on node-local disk.** The container image, gridpack and scripts are
   copied to the node first (image and gridpack into a per-node cache), because the
   shared filesystem stalls; a task touches it only at the start and for the final copy-back.
5. **Delphes run by the job script** on a stream (`gunzip -c | DelphesHepMC2`) with a
   slimmed output tree. All simulation modules still run; only the saved branches differ.
6. **Own container definition** on Ubuntu 20.04, because MG5 2.9.16 predates the Python
   in MadLAD's current base image.

## Conventions, checked against HyPER's `ttbar_reduced.root`

Settled by `validate/inspect_hyper.py` on the published file:

| Question | What HyPER's file shows | Setting |
|---|---|---|
| Jet acceptance | jets out to \|eta\| = 4.98, lowest pT 26.98 GeV: no \|eta\| < 2.5 cut, and the 25 GeV cut sits before the jet energy scale, as in MadLAD's card | `JET_ABS_ETA_MAX=5.0` |
| Which top gets labels 1-3 | the top; the stored list is always t, tbar, W+, W- | `TOP_ORDER=t-first` |
| Order of the two W quarks | quark first, antiquark second (the charm-like b-tag rate sits on labels 2 and 6) | `find_decay` |
| Stored particles | t, tbar, W+, W- | converter |
| One jet, two labels | never | `pylib/matching.py` |
| Jets sorted by pT; `nbTagged`; `allMatchedEvent` | as assumed | |
| Branch storage | `std::vector`; files written here use arrays with a counter, which uproot reads identically but ROOT-based code may not | `pylib/hyperio.py` |

Still assumed:

| Question | Assumed here | Where |
|---|---|---|
| Which copy of each quark is matched | the one produced at the top or W decay vertex, before it radiates | `find_decay` |
| Two quarks closest to one jet | the closer quark keeps the jet, the other is unmatched | `MATCH_CONFLICT`; changeable after production, see "Re-matching" |
| Whether >= 6 jets counts all jets or central jets | all jets | inspect script prints this |
| Jets beyond 20 | kept (HyPER's largest event has 20) | not capped |
| Negative-weight events | kept; sign stored in the extra branch `eventWeight` | converter |

**Open: label asymmetry.** In HyPER's file top 2 is fully matched more often than top 1
(61.8% against 57.8%), and each label of block 2 is present more often than its
counterpart in block 1. Nothing in the physics distinguishes t from tbar at this level,
so this is a property of HyPER's matching code that the converter may not reproduce.
`scripts/status.py` prints both rates for the new sample; compare them in the pilot,
and run `validate/rematch_scan.py` to see whether another matching rule reproduces them.

## Re-matching without regenerating

Per-task files in `$OUTDIR/raw` (in subfolders of 1,000 tasks: `raw/000`, `raw/001`, ...)
hold every event that passed the jet and b-tag cuts,
with the six decay quarks (b, q, qbar of the top, then of the antitop) alongside the
jets. The matched-W requirement is only a flag there. `scripts/merge.py` applies it and
drops the quarks, which gives HyPER's layout.

So the matching rule can change after production:

```
python3 validate/rematch_scan.py $OUTDIR/raw          # label rates under every rule, next to HyPER's
python3 scripts/merge.py $OUTDIR/raw $OUTDIR/merged_v2 --rematch --conflict last-wins --match-dr 0.4
```

Rules for two quarks closest to the same jet: `closest`, `first-wins`, `last-wins`,
`greedy` (the loser takes its next-closest free jet). Keep `$OUTDIR/raw` until the
matching convention is settled. What still needs new events: any change to the jet
or b-tag selection, or to which generator-record quarks are stored.

## Steps

1. **Inspect HyPER's file** (done once; rerun if the files change):
   `python3 validate/inspect_hyper.py ttbar_reduced.root`
2. **Edit `config.sh`.** Every line marked `CHECK` is a guess about HPC3.
3. **Get the container.** HPC3 does not allow container builds for ordinary accounts
   (no fakeroot mapping, and unprivileged builds fail at package installation), so the
   image is built by GitHub and pulled:
   - push this directory to a GitHub repository; the `build-container` workflow builds
     `container/Dockerfile` and publishes `ghcr.io/<owner>/madlad-hyper:latest`
   - on HPC3: `apptainer pull $WORK/madlad-hyper.sif docker://ghcr.io/<owner>/madlad-hyper:latest`

   `container/madlad-hyper.def` is the same stack for sites that do allow
   `singularity build --fakeroot`.
4. **Set up** (interactive session, 8 cores, about an hour):
   `scripts/00_setup.sh 8`
   The build runs on node-local scratch; only two tarballs and the Delphes card are
   copied to `$WORK`. Shared filesystems such as `/pub` are slow with many small files.
5. **Pilot** (interactive session): `scripts/job.sh 1 2000`, then `scripts/job.sh 2`
   for a full-size task. Then `python3 scripts/status.py $OUTDIR/raw`.
   Compare with HyPER: efficiency 14.67%, fully matched 26.7% of selected, top 1 / top 2
   matched 57.8% / 61.8%. Run `validate/rematch_scan.py $OUTDIR/raw`. Read off the
   time per event and scratch use, and adjust `WALLTIME`, `EVENTS_PER_JOB`, `CORES_PER_JOB`.
6. **First tranche = closure sample:** `scripts/submit.sh 1 2400` (60M generated, about
   8.8M selected, the size of HyPER itself). Run the validation below before going on.
7. **Production:** `scripts/submit.sh`. Finished tasks are skipped, so it is safe to rerun.
   `python3 scripts/status.py $OUTDIR/raw --tasks 21600 --missing missing.txt`
   then `scripts/submit.sh --list missing.txt` for the failures.
8. **Merge:** `python3 scripts/merge.py $OUTDIR/raw $OUTDIR/merged --holdout 440069`
   writes a held-out file the size of HyPER's test set and training shards the size of
   HyPER's training file.

## Validation

| Tier | Command | Reads as passing when |
|---|---|---|
| 0 | `scripts/status.py` | efficiency and fully-matched fraction agree with HyPER's |
| 1 | `validate/compare_spectra.py --ref ttbar_train.root --new 'merged/ttbar_new_train_*.root'` | no observable shows a deviation beyond its statistical error; check `spectra.pdf` |
| 2 | `validate/c2st.py --ref ... --new ...` (then again `--with-labels`) | new-vs-HyPER AUC lies within the HyPER-vs-HyPER spread |
| 3 | training code, not in this repository | see below |

Tier 3 decides. Train the same model and config on (a) HyPER's training file and (b) an
equal number of new events, several seeds each, and evaluate both on HyPER's test file
and on the new held-out file. The four numbers should agree within the seed spread,
inclusive and per jet multiplicity. Running the same 2x2 with the EveNet ttbar sample
shows what a real simulation shift looks like in the same units.

HyPER's test file stays the test set for every number quoted against published results.

## Not yet run

| Item | Risk |
|---|---|
| `container/madlad-hyper.def` | the Dockerfile with the same steps builds; the definition file itself has not completed a build |
| `generate_events -o -n -f` reusing a relocated process directory | standard gridpack practice, untested here; fall back to `GRID_MODE=fresh` |
| MadSpin and the shower running under `-f` | the job fails with a clear message if no showered file appears under `Events/*_decayed_*` |
| `rnd_seed` key in the shower card | checked by `00_setup.sh` before any job runs |
| Generator-record layout assumed by `find_decay` | the converter stops if more than 1% of events lack the t -> b W(-> q q) record |
| Reading real Delphes files (sub-branches `Jet/Jet.PT`) | tested only on a flat synthetic file with the same branch names |
| `cexec` flags (`--cleanenv --home`) | some sites restrict `--home` |
| Scratch use | unknown until the pilot; the HepMC file is written uncompressed before MG5 compresses it |

## Layout

```
config.sh                      all settings
container/Dockerfile           software stack, built by .github/workflows/build-container.yml
container/madlad-hyper.def     the same stack as a Singularity definition
process/ttbar-allhad-hyper.yaml   MadLAD process config
scripts/00_setup.sh            one-time: MadLAD, process directory, grids
scripts/job.sh                 one task
scripts/submit.sh              slurm arrays
scripts/status.py              bookkeeping, sizing, missing tasks
scripts/merge.py               shards in HyPER's layout
scripts/set_card.py, make_slim_card.py   helpers
convert/delphes_to_hyper.py    selection, truth matching, output
validate/                      inspect_hyper.py, rematch_scan.py, compare_spectra.py, c2st.py
pylib/hyperio.py               file layout, shared
pylib/matching.py              jet-quark matching rules, shared by converter and merge
tests/run_all.sh               tests that need no MadGraph
```
