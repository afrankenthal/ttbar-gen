#!/usr/bin/env python3
"""Tier 1 closure: compare distributions between the HyPER files and the new sample.

usage: compare_spectra.py --ref ttbar_train.root --new 'merged/ttbar_new_train_*.root'
                          [--entries N] [--out DIR]

Writes <out>/spectra.json (chi2/ndf and the largest bin-by-bin deviation for every
observable, with the statistical error on that deviation) and <out>/spectra.pdf (one page per observable, with a ratio panel).

Reading the numbers: with millions of events chi2/ndf detects differences far below
anything that matters, so look at "max dev" too: the largest relative difference in
any bin holding at least 1% of the events. Selection efficiency is not in these
files; scripts/status.py reports it from the job JSONs.
"""
import argparse
import glob
import json
import os
import sys

import awkward as ak
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pylib"))
import hyperio  # noqa: E402

C_REF, C_NEW = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"


def load(pattern, entries):
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit("compare: no files match %s" % pattern)
    chunks, n = [], 0
    for f in files:
        chunks.append(hyperio.read(f, entry_stop=None if entries is None else entries - n))
        n += len(chunks[-1]["jet_pt"])
        if entries is not None and n >= entries:
            break
    return hyperio.concatenate(chunks) if len(chunks) > 1 else chunks[0]


def mass_of(d, labels):
    """Invariant mass of the jets carrying the given labels; NaN where a label is missing."""
    lab = d["jet_truthmatch"]
    sel = lab == labels[0]
    for k in labels[1:]:
        sel = sel | (lab == k)
    ok = np.ones(len(lab), dtype=bool)
    for k in labels:
        ok &= np.asarray(ak.sum(lab == k, axis=1)) == 1
    pt, eta, phi, e = (d["jet_" + k][sel] for k in ["pt", "eta", "phi", "e"])
    px, py, pz, en = (np.asarray(ak.sum(x, axis=1)) for x in [pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta), e])
    m = np.sqrt(np.maximum(en ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))
    return m[ok]


def nth(arr, i):
    return np.asarray(ak.fill_none(ak.pad_none(arr, i + 1, axis=1)[:, i], np.nan))


def observables(d):
    lab, pt = d["jet_truthmatch"], d["jet_pt"]
    flat = lambda x: np.asarray(ak.flatten(x))  # noqa: E731
    tagged = d["jet_bTag"] > 0
    ibins = lambda lo, hi: np.arange(lo - 0.5, hi + 1.5)  # noqa: E731
    ptbins = np.geomspace(25, 1000, 41)
    return [
        ("jet multiplicity", np.asarray(ak.num(pt)), ibins(6, 20), False),
        ("b-tagged jet multiplicity", np.asarray(ak.sum(tagged, axis=1)), ibins(2, 8), False),
        ("matched quarks per event", np.asarray(ak.sum(lab > 0, axis=1)), ibins(2, 6), False),
        ("jet pT [GeV]", flat(pt), ptbins, True),
        ("leading jet pT [GeV]", nth(pt, 0), ptbins, True),
        ("6th jet pT [GeV]", nth(pt, 5), np.geomspace(25, 300, 31), True),
        ("HT [GeV]", np.asarray(ak.sum(pt, axis=1)), np.geomspace(200, 4000, 41), True),
        ("jet eta", flat(d["jet_eta"]), np.linspace(-2.5, 2.5, 51), False),
        ("jet phi", flat(d["jet_phi"]), np.linspace(-np.pi, np.pi, 33), False),
        ("jet mass [GeV]", flat(d["jet_m"]), np.linspace(0, 60, 41), False),
        ("b-tagged jet pT [GeV]", flat(pt[tagged]), ptbins, True),
        ("pT of jets matched to b quarks [GeV]", flat(pt[(lab == 1) | (lab == 4)]), ptbins, True),
        ("pT of jets matched to W quarks [GeV]", flat(pt[(lab == 2) | (lab == 3) | (lab == 5) | (lab == 6)]), ptbins, True),
        ("pT of unmatched jets [GeV]", flat(pt[lab == 0]), ptbins, True),
        ("eta of unmatched jets", flat(d["jet_eta"][lab == 0]), np.linspace(-2.5, 2.5, 51), False),
        ("W mass, top 1 (labels 2,3) [GeV]", mass_of(d, [2, 3]), np.linspace(40, 140, 51), False),
        ("W mass, top 2 (labels 5,6) [GeV]", mass_of(d, [5, 6]), np.linspace(40, 140, 51), False),
        ("top mass, top 1 (labels 1,2,3) [GeV]", mass_of(d, [1, 2, 3]), np.linspace(100, 260, 41), False),
        ("top mass, top 2 (labels 4,5,6) [GeV]", mass_of(d, [4, 5, 6]), np.linspace(100, 260, 41), False),
    ]


def rates(d):
    """Per-event and per-jet rates worth comparing as single numbers."""
    lab = d["jet_truthmatch"]
    has = [np.asarray(ak.any(lab == k, axis=1)) for k in range(1, 7)]
    t1, t2 = has[0] & has[1] & has[2], has[3] & has[4] & has[5]
    tagged = d["jet_bTag"] > 0
    f = lambda x: (float(np.sum(x)), float(len(x)))  # noqa: E731
    fl = lambda x: np.asarray(ak.flatten(x))  # noqa: E731
    return {
        "fully matched events": f(t1 & t2),
        "top 1 fully matched": f(t1),
        "top 2 fully matched": f(t2),
        "b-tag rate, jets matched to b quarks": f(fl(tagged[(lab == 1) | (lab == 4)])),
        "b-tag rate, jets matched to W quarks": f(fl(tagged[(lab == 2) | (lab == 3) | (lab == 5) | (lab == 6)])),
        "b-tag rate, unmatched jets": f(fl(tagged[lab == 0])),
    }


def shape_compare(a, b, bins):
    na, _ = np.histogram(a[np.isfinite(a)], bins)
    nb, _ = np.histogram(b[np.isfinite(b)], bins)
    ta, tb = float(na.sum()), float(nb.sum())
    use = (na + nb) > 0
    pa, pb = na / ta, nb / tb
    var = na / ta ** 2 + nb / tb ** 2
    chi2 = float(np.sum((pa[use] - pb[use]) ** 2 / var[use]))
    ndf = int(use.sum()) - 1
    big = np.flatnonzero((pa > 0.01) & (pb > 0))
    if big.size == 0:
        return na, nb, chi2, ndf, float("nan"), float("nan")
    dev = pb[big] / pa[big] - 1.0
    k = big[np.argmax(np.abs(dev))]
    maxdev = float(np.max(np.abs(dev)))
    err = float((pb[k] / pa[k]) * np.sqrt(1.0 / na[k] + 1.0 / nb[k]))     # statistical error on that ratio
    return na, nb, chi2, ndf, maxdev, err


def plot_page(pdf, plt, name, bins, na, nb, logx, chi2, ndf, maxdev, n_ref, n_new):
    ta, tb = na.sum(), nb.sum()
    width = np.diff(bins)
    pa, pb = na / ta / width, nb / tb / width
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True, facecolor=SURFACE,
                                 gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06})
    for axis in (ax, rx):
        axis.set_facecolor(SURFACE)
        axis.grid(True, axis="y", color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            axis.spines[side].set_color(AXIS)
        axis.tick_params(colors=MUTED, labelsize=9)
    ax.stairs(pa, bins, color=C_REF, linewidth=2, label="HyPER (%s events)" % format(n_ref, ","))
    ax.stairs(pb, bins, color=C_NEW, linewidth=2, linestyle=(0, (4, 2)), label="new sample (%s events)" % format(n_new, ","))
    ax.set_ylabel("normalised density", color=INK2, fontsize=10)
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    ax.set_title("%s\nchi2/ndf = %.1f/%d, largest deviation in populated bins %.1f%%"
                 % (name, chi2, ndf, 100 * maxdev), loc="left", color=INK, fontsize=11)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(pa > 0, pb / pa, np.nan)
        err = np.where((na > 0) & (nb > 0), ratio * np.sqrt(1.0 / np.maximum(na, 1) + 1.0 / np.maximum(nb, 1)), np.nan)
    centres = np.sqrt(bins[:-1] * bins[1:]) if logx else 0.5 * (bins[:-1] + bins[1:])
    rx.axhline(1.0, color=AXIS, linewidth=1)
    rx.errorbar(centres, ratio, yerr=err, fmt="o", markersize=3.5, color=C_NEW, ecolor=C_NEW, elinewidth=1)
    rx.set_ylim(0.9, 1.1)
    rx.set_ylabel("new / HyPER", color=INK2, fontsize=10)
    rx.set_xlabel(name, color=INK2, fontsize=10)
    if logx:
        rx.set_xscale("log")
    pdf.savefig(fig, facecolor=SURFACE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="HyPER file (glob allowed)")
    ap.add_argument("--new", required=True, help="new-sample files (glob allowed)")
    ap.add_argument("--entries", type=int, default=2000000, help="events read from each side")
    ap.add_argument("--out", default="closure_tier1")
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ref, new = load(a.ref, a.entries), load(a.new, a.entries)
    n_ref, n_new = len(ref["jet_pt"]), len(new["jet_pt"])
    print("compare: %d reference events, %d new events" % (n_ref, n_new))

    pdf = plt = None
    if not a.no_plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(os.path.join(a.out, "spectra.pdf"))

    results = []
    for (name, xa, bins, logx), (_, xb, _, _) in zip(observables(ref), observables(new)):
        na, nb, chi2, ndf, maxdev, deverr = shape_compare(xa, xb, bins)
        results.append({"observable": name, "chi2": chi2, "ndf": ndf,
                        "max_dev_populated_bins": maxdev, "max_dev_stat_error": deverr})
        if pdf is not None:
            plot_page(pdf, plt, name, bins, na, nb, logx, chi2, ndf, maxdev, n_ref, n_new)
    if pdf is not None:
        pdf.close()

    rate_rows = []
    ra, rb = rates(ref), rates(new)
    for name in ra:
        (ka, na_), (kb, nb_) = ra[name], rb[name]
        pa, pb = ka / na_, kb / nb_
        err = np.sqrt(pa * (1 - pa) / na_ + pb * (1 - pb) / nb_)
        rate_rows.append({"rate": name, "ref": pa, "new": pb, "difference": pb - pa,
                          "sigma": float((pb - pa) / err) if err > 0 else float("nan")})

    with open(os.path.join(a.out, "spectra.json"), "w") as f:
        json.dump({"n_ref": n_ref, "n_new": n_new, "shapes": results, "rates": rate_rows}, f, indent=1)

    print("\n%-42s %10s %18s" % ("observable", "chi2/ndf", "max dev (+- stat)"))
    for r in sorted(results, key=lambda r: -r["chi2"] / max(r["ndf"], 1)):
        print("%-42s %10.2f %9.1f%% +- %.1f%%" % (r["observable"], r["chi2"] / max(r["ndf"], 1),
                                                 100 * r["max_dev_populated_bins"], 100 * r["max_dev_stat_error"]))
    print("\n%-42s %9s %9s %9s %7s" % ("rate", "HyPER", "new", "diff", "sigma"))
    for r in rate_rows:
        print("%-42s %8.2f%% %8.2f%% %+8.2f%% %7.1f" % (r["rate"], 100 * r["ref"], 100 * r["new"], 100 * r["difference"], r["sigma"]))
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
