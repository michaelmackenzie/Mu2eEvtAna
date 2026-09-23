#!/usr/bin/env python3
"""Process an EventNtuple dataset with pyevtana and write a Mu2eEvtAna-style ROOT file.

The python counterpart of ``Mu2eEvtAna/scripts/make_histograms.C``: take a dataset key,
find its file list, run the selection sets over it with N workers, and write one ROOT file
holding a histogram book per selection set (``trk_<id>`` / ``evt_<id>``).

    # setup
    source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh
    pyenv rootana 2.5.0
    export PYTHONPATH=/exp/mu2e/app/users/mmackenz/main/pyevtana:$PYTHONPATH

    # list what is available
    python3 make_histograms.py --list

    # one dataset, 10 workers, first 20 files
    python3 make_histograms.py cele1b1s5r0100 --jobs 10 --max-files 20

    # write one ROOT file per worker and merge them with scripts/merge_hist_files.sh
    python3 make_histograms.py cele1b1s5r0100 --jobs 10 --split-files

Workers return numpy counts by default, which merge by addition, so no temporary files and
no hadd are needed; ``--split-files`` is there for when a per-worker file really is wanted.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import datasets as datasets_module
import histograms
import selections
from pyevtana import Dataset

#: the set whose survivors define "accepted" in the Norm tree
PRIMARY_SET = 75

MU2E_EVT_ANA = datasets_module.MU2E_EVT_ANA
MERGE_SCRIPT = os.path.join(MU2E_EVT_ANA, "scripts", "merge_hist_files.sh")


# --------------------------------------------------------------------------------------
# the per-partition worker
# --------------------------------------------------------------------------------------


def process(chunk):
    """Fill every selection set's book over one partition. Module level, so it pickles."""
    set_ids = sorted(selections.SELECTIONS)
    books = histograms.Books(set_ids)
    flows = {set_id: [0] * len(selections.SELECTIONS[set_id].cuts) for set_id in set_ids}
    n_events = 0
    n_accept = 0

    for entry in chunk:
        n_events += 1
        event = selections.Event(entry)
        for set_id in set_ids:
            selection = selections.SELECTIONS[set_id]
            selected, per_stage = selection.select(event)
            for i, n in enumerate(per_stage):          # per-event cumulative cut flow
                if n:
                    flows[set_id][i] += 1
            books[set_id].fill(event, selected)
            if set_id == PRIMARY_SET and selected:
                n_accept += 1

    books.flush()          # bin the accumulated values once, not once per event

    # every entry is read, so nseen == nntuple here; nneg stays 0 until weights are applied
    norm = histograms.Norm(ngen=0, nntuple=n_events, nseen=n_events,
                           naccept=n_accept, nneg=0)
    return {"books": books, "flows": flows, "n_events": n_events, "norm": norm}


def merge(results):
    results = [r for r in results if r]
    if not results:
        return None
    total = results[0]
    for other in results[1:]:
        total["books"] += other["books"]
        total["n_events"] += other["n_events"]
        total["norm"] += other["norm"]
        for set_id, flow in other["flows"].items():
            total["flows"][set_id] = [a + b for a, b in zip(total["flows"][set_id], flow)]
    return total


def process_and_write(chunk):
    """--split-files worker: write this partition's book to its own ROOT file."""
    result = process(chunk)
    path = os.path.join(chunk.dataset._split_dir,
                        f"part_{chunk.partition.index:04d}.root")
    # each part carries its own Norm entry; hadd concatenates them, so the merged tree has
    # one entry per part and the totals are the sums -- the same shape Mu2eEvtAna produces
    # when its jobs are merged
    norm = result["norm"]
    norm.ngen = int(round(chunk.dataset._gen_per_event * norm.nntuple))
    result["books"].write(path, norm)
    return {"path": path, "flows": result["flows"], "n_events": result["n_events"],
            "norm": norm}


def merge_paths(results):
    results = [r for r in results if r]
    if not results:
        return None
    total = {"paths": [], "flows": results[0]["flows"], "n_events": 0,
             "norm": histograms.Norm()}
    for i, r in enumerate(results):
        total["paths"].append(r["path"])
        total["n_events"] += r["n_events"]
        total["norm"] += r["norm"]
        if i:
            for set_id, flow in r["flows"].items():
                total["flows"][set_id] = [a + b for a, b in zip(total["flows"][set_id], flow)]
    return total


# --------------------------------------------------------------------------------------


def report(dataset, result, wall):
    n = result["n_events"]
    norm = result["norm"]
    print(f"\nnormalization (Ana/data/Norm)")
    print(f"  ngen     {norm.ngen:>12}   N(generated events)")
    print(f"  nntuple  {norm.nntuple:>12}   N(events in the input ntuple)")
    print(f"  nseen    {norm.nseen:>12}   N(processed events)")
    print(f"  naccept  {norm.naccept:>12}   N(accepted events, set {PRIMARY_SET})")
    print(f"  nneg     {norm.nneg:>12}   N(negative weight events)")
    for set_id in sorted(selections.SELECTIONS):
        selection = selections.SELECTIONS[set_id]
        flow = result["flows"][set_id]
        print(f"\nset {set_id} — {selection.name}   (events with >= 1 surviving track)")
        width = max(len(s) for s in selection.stages)
        previous = result["n_events"]
        for stage, count in zip(selection.stages, flow):
            rel = f"{100 * count / previous:6.1f}%" if previous else "   n/a"
            print(f"  {stage:<{width}}  {count:>9}  {rel}")
            previous = count


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset", nargs="?", help="dataset key, e.g. cele1b1s5r0100")
    parser.add_argument("--list", action="store_true", help="list known dataset keys")
    parser.add_argument("--jobs", type=int, default=1, help="worker processes")
    parser.add_argument("--max-files", type=int, default=0,
                        help="process only the first N files (0 = all)")
    parser.add_argument("--output", default=None, help="output ROOT file")
    parser.add_argument("--outdir", default=os.path.join(MU2E_EVT_ANA, "hist"),
                        help="directory for the output file")
    parser.add_argument("--split-files", action="store_true",
                        help="one ROOT file per worker, merged with merge_hist_files.sh")
    parser.add_argument("--keep-parts", action="store_true",
                        help="with --split-files, do not delete the per-worker files")
    args = parser.parse_args()

    if args.list or not args.dataset:
        print(f"{'key':<18}{'files':>7}  full name")
        for d in datasets_module.DATASETS:
            present = os.path.exists(d.file_list)
            n = len(d.files()) if present else 0
            flag = " " if present else "*"
            print(f"{d.key:<18}{n:>7}{flag} {d.full_name}")
        print("\n* no file list in file_lists/")
        return 0

    spec = datasets_module.get(args.dataset)
    files = spec.files()
    if args.max_files:
        files = files[:args.max_files]
    print(f"{spec.key}: {len(files)} file(s) from {os.path.basename(spec.file_list)}")
    print(f"selection sets: {', '.join(str(i) for i in sorted(selections.SELECTIONS))}")

    os.makedirs(args.outdir, exist_ok=True)
    output = args.output or os.path.join(args.outdir, f"{spec.key}.hist.root")

    # N(generated) comes from the dataset table rather than the ntuple; when only part of
    # the file list is processed it is scaled by the fraction of ntuple events seen, which
    # is an assumption worth seeing printed rather than buried
    gen_per_event = (spec.n_gen_events / spec.n_events) if spec.n_events else 0.0
    partial = len(files) < len(spec.files())

    dataset = Dataset(files, branches=selections.BRANCHES)
    dataset._gen_per_event = gen_per_event
    backend = "serial" if args.jobs == 1 else "process"

    if args.split_files:
        split_dir = os.path.join(args.outdir, f"{spec.key}_parts")
        os.makedirs(split_dir, exist_ok=True)
        dataset._split_dir = split_dir
        start = time.time()
        result = dataset.map(process_and_write, workers=args.jobs, backend=backend,
                             reduce=merge_paths)
        wall = time.time() - start
        listing = os.path.join(split_dir, "parts.txt")
        with open(listing, "w") as handle:
            handle.write("\n".join(result["paths"]) + "\n")
        command = [MERGE_SCRIPT, output, listing, "--force"]
        if not args.keep_parts:
            command.append("--clean")
        print(f"\nmerging {len(result['paths'])} files with merge_hist_files.sh")
        subprocess.run(command, check=True)
    else:
        start = time.time()
        result = dataset.map(process, workers=args.jobs, backend=backend, reduce=merge)
        wall = time.time() - start
        result["norm"].ngen = int(round(gen_per_event * result["norm"].nntuple))
        result["books"].write(output, result["norm"])

    if partial:
        print(f"\nnote: {len(files)} of {len(spec.files())} files processed, so ngen is "
              f"scaled from the dataset total ({spec.n_gen_events:,} for {spec.n_events:,} "
              f"ntuple events)")
    report(spec, result, wall)
    n = result["n_events"]
    rate = n / wall if wall else 0.0
    per_event = 1e6 * wall / n if n else 0.0
    print(f"\nprocessed {n} events in {wall:.1f} s")
    print(f"  {rate:,.0f} events/s  ({per_event:.1f} us/event, {args.jobs} worker(s))")
    if args.jobs > 1:
        print(f"  {rate/args.jobs:,.0f} events/s ({per_event * args.jobs:.1f} us/event) per worker")
    print(f"\nwrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
