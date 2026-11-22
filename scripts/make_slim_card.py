#!/usr/bin/env python3
"""Write a copy of a Delphes card whose output tree keeps only what the converter reads.

Only the TreeWriter block changes. Every simulation and reconstruction module still
runs in the same order, so the jets are the same as with the full card; the output
file is just much smaller (no tracks, towers or energy-flow objects).

usage: make_slim_card.py <input card> <output card>
"""
import re
import sys

KEEP = [
    "  add Branch Delphes/allParticles Particle GenParticle",
    "  add Branch UniqueObjectFinder/jets Jet Jet",
    "  add Branch UniqueObjectFinder/electrons Electron Electron",
    "  add Branch UniqueObjectFinder/muons Muon Muon",
    "  add Branch MissingET/momentum MissingET MissingET",
    "  add Branch ScalarHT/energy ScalarHT ScalarHT",
]


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    with open(argv[1]) as f:
        text = f.read()
    pat = re.compile(r"module\s+TreeWriter\s+TreeWriter\s*\{.*?\n\}", re.DOTALL)
    blocks = pat.findall(text)
    if len(blocks) != 1:
        sys.exit("make_slim_card: expected one TreeWriter block, found %d" % len(blocks))
    for line in KEEP:
        if line.strip() not in [b.strip() for b in blocks[0].splitlines()]:
            sys.exit("make_slim_card: the input card has no line '%s'; "
                     "the card differs from the one this script was written for" % line.strip())
    new_block = "module TreeWriter TreeWriter {\n# slimmed by make_slim_card.py\n" + "\n".join(KEEP) + "\n}"
    with open(argv[2], "w") as f:
        f.write(pat.sub(lambda m: new_block, text))
    print("make_slim_card: wrote %s" % argv[2])


if __name__ == "__main__":
    main(sys.argv)
