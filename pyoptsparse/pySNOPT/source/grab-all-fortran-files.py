#!/usr/bin/env python

# Standard Python modules
import os
from pathlib import Path
from argparse import ArgumentParser

CURDIR = os.path.abspath(os.path.dirname(__file__))


def main():
    parser = ArgumentParser()
    parser.add_argument("sourcedir", type=Path, default=None)
    args = parser.parse_args()
    if args.sourcedir is None:
        args.sourcedir = CURDIR

    TO_SKIP = [
        "sn27lu77.f",
        "sn27lu90.f",
        "snopth.f",
    ]

    for path in Path(args.sourcedir).glob("*.f"):
        if os.path.basename(path) not in TO_SKIP:
            print(path)


if __name__ == "__main__":
    main()
