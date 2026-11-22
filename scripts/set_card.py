#!/usr/bin/env python3
"""Set values in MG5_aMC cards, and fail loudly if a key is not there.

  run card     lines look like   "  25000 = nevents ! comment"   (value = key)
  shower card  lines look like   "rnd_seed = 0  # comment"       (key = value)

usage: set_card.py run    <card> key=value [key=value ...]
       set_card.py shower <card> key=value [key=value ...]
"""
import re
import sys


def set_run(lines, key, value):
    pat = re.compile(r"^(\s*)(\S+)(\s*=\s*)" + re.escape(key) + r"(\s*(!.*)?)$", re.IGNORECASE)
    hits = 0
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#"):
            continue
        m = pat.match(line.rstrip("\n"))
        if m:
            lines[i] = "%s%s%s%s%s\n" % (m.group(1), value, m.group(3), key, m.group(4))
            hits += 1
    return hits


def set_shower(lines, key, value):
    pat = re.compile(r"^(\s*)" + re.escape(key) + r"(\s*=\s*)(\S+)(.*)$", re.IGNORECASE)
    hits = 0
    for i, line in enumerate(lines):
        if line.lstrip().startswith("#"):
            continue
        m = pat.match(line.rstrip("\n"))
        if m:
            lines[i] = "%s%s%s%s%s\n" % (m.group(1), key, m.group(2), value, m.group(4))
            hits += 1
    return hits


def main(argv):
    if len(argv) < 4 or argv[1] not in ("run", "shower"):
        sys.exit(__doc__)
    kind, path, pairs = argv[1], argv[2], argv[3:]
    with open(path) as f:
        lines = f.readlines()
    setter = set_run if kind == "run" else set_shower
    for pair in pairs:
        key, _, value = pair.partition("=")
        hits = setter(lines, key, value)
        if hits != 1:
            sys.exit("set_card: expected exactly one '%s' line in %s, found %d" % (key, path, hits))
        print("set_card: %s = %s (%s)" % (key, value, path))
    with open(path, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    main(sys.argv)
