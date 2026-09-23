"""Dataset key -> full dataset name, mirroring ``Mu2eEvtAna/scripts/datasets.C``.

A key such as ``cele1b1s5r0100`` maps to a full dataset name, whose file list lives at
``Mu2eEvtAna/file_lists/<full_name>.files``. Keep this table in step with ``datasets.C``;
the keys, names and event counts are copied from it.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

#: repository root, so file lists resolve wherever the script is run from
MU2E_EVT_ANA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILE_LISTS = os.path.join(MU2E_EVT_ANA, "file_lists")


@dataclass(frozen=True)
class Dataset:
    """One row of ``datasets.C``."""

    key: str
    full_name: str
    n_events: int
    n_gen_events: float
    process: bool = True

    @property
    def file_list(self) -> str:
        return os.path.join(FILE_LISTS, f"{self.full_name}.files")

    def files(self) -> list:
        """Paths in this dataset's file list."""
        if not os.path.exists(self.file_list):
            raise FileNotFoundError(
                f"no file list for {self.key!r} at {self.file_list}")
        with open(self.file_list) as handle:
            return [line.strip() for line in handle
                    if line.strip() and not line.startswith("#")]


_ELE, _POS = False, True          # which datasets each analysis enables, as in datasets.C

DATASETS = [
    Dataset("cele1b0s5r0100", "nts.mu2e.CeMLeadingLogOnSpill-reco-ntuple.MDC2025-002.root",     4009075,   10000000, False),
    Dataset("cry4ab0s5r0100", "nts.mu2e.CosmicSignalOnSpill-reco-ntuple.MDC2025-002.root",      4120241,    4437500, False),
    Dataset("dio00b0s5r0100", "nts.mu2e.DIOtail95OnSpill-reco-ntuple.MDC2025-002.root",         8780533,   25000000, False),

    Dataset("cele1b1s5r0100", "nts.mu2e.CeMLeadingLogMix1BB.MDC2025au_best_v1_1.root",          4141125,   10000000, _ELE),
    Dataset("cpos1b1s5r0100", "nts.mu2e.CePLeadingLogMix1BB.MDC2025au_best_v1_1.root",          3235878,   10000000, _POS),
    Dataset("cry4ab1s5r0100", "nts.mu2e.CosmicSignalMix1BB.MDC2025au_best_v1_1.root",           4155435,    4437713, True),
    Dataset("dio00b1s5r0100", "nts.mu2e.DIOtail95Mix1BB.MDC2025au_best_v1_1.root",              9368976,   25000000, _ELE),
    Dataset("rmce0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NExternalMix1BB.MDC2025au_best_v1_1.root", 4967393, 7000000000, True),
    Dataset("rmce1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NExternalMix1BB.MDC2025au_best_v1_1.root", 2974188, 7000000000, _POS),
    Dataset("rmci0b1s5r0100", "nts.mu2e.RMCPhaseSpace0NInternalMix1BB.MDC2025au_best_v1_1.root", 1229006,   50000000, _POS),
    Dataset("rmci1b1s5r0100", "nts.mu2e.RMCPhaseSpace1NInternalMix1BB.MDC2025au_best_v1_1.root",  517029,   50000000, _POS),
    Dataset("rpce1b1s5r0100", "nts.mu2e.RPCExternalPhysicalMix1BB.MDC2025au_best_v1_1.root",      458818, 5000000000, True),
    Dataset("rpci1b1s5r0100", "nts.mu2e.RPCInternalPhysicalMix1BB.MDC2025au_best_v1_1.root",     1899806,  125000000, True),
    Dataset("pbar1b1s5r0100", "nts.mu2e.PbarResamplingMix1BB.MDC2025au_best_v1_1.root",          6461314,   30000000, _POS),

    Dataset("mds3cb1s5r0100", "nts.mu2e.ensembleMDS3cMix1BB.MDC2025au_best_v1_1.root",          4926941,          1, False),
]

BY_KEY = {d.key: d for d in DATASETS}


def get(key: str) -> Dataset:
    if key not in BY_KEY:
        raise KeyError(f"unknown dataset key {key!r}; known keys: "
                       f"{', '.join(sorted(BY_KEY))}")
    return BY_KEY[key]


def enabled() -> list:
    """Datasets flagged for processing, as ``datasets.C``'s ``process_`` does."""
    return [d for d in DATASETS if d.process]


def available(key: Optional[str] = None) -> list:
    """Datasets whose file list actually exists on disk."""
    pool = [get(key)] if key else DATASETS
    return [d for d in pool if os.path.exists(d.file_list)]
