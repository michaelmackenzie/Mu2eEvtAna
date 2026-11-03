#!/usr/bin/env python3
"""Make file lists and normalization counts for ntuple datasets.

For each ntuple dataset (e.g. nts.mu2e.CeEndpointMix1BB-KL.Run1B-013.root) this
  - writes file_lists/<dataset>.files with one /pnfs path per file, and
  - sums N(events) (rse.nevent) and N(generated) (gen.count) over the parent mcs files of the ntuple files,
then prints a datasets.C line to edit. For a partly produced dataset the counts cover only the ntuple files that
exist, and the totals of the full parent dataset are printed alongside.

Needs the metacat client and a metacat login (`muse setup ops`; `getToken`).

Usage:
  python3 Mu2eEvtAna/scripts/make_file_lists.py nts.mu2e.CeEndpointMix1BB-KL.Run1B-013.root [more datasets ...]
  python3 Mu2eEvtAna/scripts/make_file_lists.py --no-lists <datasets>   # counts only
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys

FILE_LIST_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "file_lists")
LOCATIONS = {"tape": "/pnfs/mu2e/tape", "disk": "/pnfs/mu2e/persistent/datasets"}
SEQUENCER = re.compile(r"\.(\d{6}_\d{8})\.\w+$")


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def query(mql, metadata=False):
    cmd = ["metacat", "query", "-j"] + (["-m", "all"] if metadata else []) + [mql]
    return json.loads(run(cmd))


def pnfs_path(name, location):
    """Standard Mu2e dCache path of a file: <root>/phy-<tier>/<tier>/<owner>/<desc>/<config>/<format>/xx/yy/<name>"""
    tier, owner, desc, config, _, fmt = name.split(".")
    h = hashlib.sha256(name.encode()).hexdigest()
    return f"{LOCATIONS[location]}/phy-{tier}/{tier}/{owner}/{desc}/{config}/{fmt}/{h[:2]}/{h[2:4]}/{name}"


def parent_dataset(file_name):
    """The mcs parent dataset of an ntuple file (without the .art extension)"""
    shown = run(["metacat", "file", "show", "-l", f"mu2e:{file_name}"])
    match = re.search(r"mu2e:(mcs\.\S+)\.\d{6}_\d{8}\.art", shown)
    if not match:
        raise RuntimeError(f"No mcs parent found for {file_name}")
    return match.group(1)


def process(dataset, args):
    names = sorted(f["name"] for f in query(f"files from mu2e:{dataset}"))
    if not names:
        print(f"{dataset}: no files", file=sys.stderr)
        return

    # Counts from the parent mcs files, matched by sequencer
    parent = parent_dataset(names[0])
    meta = {SEQUENCER.search(f["name"]).group(1): f["metadata"]
            for f in query(f"files from mu2e:{parent}.art", metadata=True)}
    seqs = [SEQUENCER.search(n).group(1) for n in names]
    missing = [s for s in seqs if s not in meta]
    if missing:
        print(f"{dataset}: WARNING: {len(missing)} files have no parent in {parent}.art (e.g. {missing[0]})", file=sys.stderr)
    n_events = sum(meta[s].get("rse.nevent", 0) for s in seqs if s in meta)
    n_gen = sum(meta[s].get("gen.count", 0) for s in seqs if s in meta)

    if not args.no_lists:
        location = args.location
        if location is None:  # use disk if the first file is there
            location = "disk" if os.path.exists(pnfs_path(names[0], "disk")) else "tape"
        out = os.path.join(args.out_dir, f"{dataset}.files")
        if os.path.exists(out) and not args.force:
            print(f"{dataset}: {out} exists, not overwritten (use --force)", file=sys.stderr)
        else:
            paths = [pnfs_path(n, location) for n in names]
            if not os.path.exists(paths[0]):
                print(f"{dataset}: WARNING: {paths[0]} not found", file=sys.stderr)
            with open(out, "w") as f:
                f.write("\n".join(paths) + "\n")
            print(f"{dataset}: wrote {out} ({location})")

    status = "complete" if len(names) == len(meta) else f"PARTIAL: {len(names)} of {len(meta)} parent files"
    print(f"{dataset}: {len(names)} files, parent {parent}.art, {status}")
    if len(names) != len(meta):
        print(f"  full parent dataset: N(events) = {sum(m.get('rse.nevent', 0) for m in meta.values())}, "
              f"N(gen) = {sum(m.get('gen.count', 0) for m in meta.values())}")
    print(f'  datasets.emplace_back(true , "<key>", "{dataset}", {n_events:>10}, {n_gen:>10});')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="+", help="ntuple dataset names (nts.<owner>.<desc>.<config>.root)")
    parser.add_argument("--out-dir", default=FILE_LIST_DIR, help="file list directory (default: Mu2eEvtAna/file_lists)")
    parser.add_argument("--location", choices=LOCATIONS, help="dCache location (default: disk if present, else tape)")
    parser.add_argument("--force", action="store_true", help="overwrite existing file lists")
    parser.add_argument("--no-lists", action="store_true", help="only print the counts")
    args = parser.parse_args()
    for dataset in args.datasets:
        process(dataset.removeprefix("mu2e:"), args)


if __name__ == "__main__":
    main()
