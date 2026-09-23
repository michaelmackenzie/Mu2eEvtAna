"""Histogram books, one per selection set ID.

Mirrors ``Mu2eEvtAna::BookHistograms``: a book per selection set, written to directories
``trk_<id>`` and ``evt_<id>``, with the names and binning ``Mu2eEvtAna`` uses so the
existing plotting macros (``mumep_ana/analysis/make_plots.C``) read the output unchanged.

A deliberately small subset to begin with -- the ones make_plots.C leans on most. Adding
one is a single row in :data:`TRACK_HISTS`.

Filling accumulates into numpy arrays rather than ROOT objects, because a worker process
has to send its results back: numpy histograms are picklable and merge by addition, so
N workers need no temporary ROOT files and no hadd. ROOT is touched only once, by the
parent, when the merged book is written out.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Norm:
    """Normalization counts, written as the ``Norm`` tree under ``Ana/data``.

    Branches and meanings follow ``Mu2eEvtAna.cc`` (``Norm_t.hh``), all ``Long64_t``:

    ``ngen``     N(generated events)
    ``nntuple``  N(events in the input ntuple)
    ``nseen``    N(processed events)
    ``naccept``  N(accepted events)
    ``nneg``     N(negative weight events)
    """

    ngen: int = 0
    nntuple: int = 0
    nseen: int = 0
    naccept: int = 0
    nneg: int = 0

    def __iadd__(self, other: "Norm") -> "Norm":
        self.ngen += other.ngen
        self.nntuple += other.nntuple
        self.nseen += other.nseen
        self.naccept += other.naccept
        self.nneg += other.nneg
        return self

#: (name, title, nbins, low, high, value) -- value pulls the quantity off a Track
TRACK_HISTS = [
    ("p",         "Track momentum",              300,    0.,  150., lambda t: t.p),
    ("p_2",       "Track momentum",              600,   80.,  110., lambda t: t.p),
    ("pt",        "track transverse momentum",   300,    0.,  300., lambda t: t.front_pt),
    ("t0",        "track t_{0}",                 400,    0., 2000., lambda t: t.t0),
    ("d0",        "track d0",                    200, -200.,  200., lambda t: t.d0),
    ("rMax",      "track rMax",                 2000,    0., 2000., lambda t: t.rmax),
    ("tanDip",    "track tanDip",                200,    0.,    2., lambda t: t.tandip),
    ("cosTheta",  "track cos(#theta)",           200,   -1.,    1., lambda t: t.front_cos),
    ("nActive",   "nHits used in fit",           150,    0.,  150., lambda t: t.nactive),
    ("trkQual",   "track MVA score",             200,   -1.,    1., lambda t: t.trkqual),
    ("trkpid",    "TrkPID MVA score",            200,   -1.,    1., lambda t: t.trkpid),
    ("clusterE",  "track's cluster energy",      600,    0.,  300., lambda t: t.cluster_e),
    ("ep",        "cluster E / track P",         200,    0.,    2., lambda t: t.ep),
]

#: event-level histograms, filled once per event that has >= 1 selected track
EVENT_HISTS = [
    ("nTracks", "nTracks", 50, 0., 50., lambda e, selected: len(selected)),
]


class Book:
    """The histograms for one selection set, as numpy counts.

    Values are appended to plain lists as events go by and binned in one pass by
    :meth:`flush`. Binning a handful of values per event instead would spend most of its
    time in numpy call overhead -- a few microseconds per call, times a histogram per
    selection set per event, is more than reading the data costs.
    """

    __slots__ = ("set_id", "track", "event", "n_events", "n_tracks", "_pending")

    def __init__(self, set_id: int):
        self.set_id = set_id
        self.track = {spec[0]: np.zeros(spec[2], dtype=np.float64) for spec in TRACK_HISTS}
        self.event = {spec[0]: np.zeros(spec[2], dtype=np.float64) for spec in EVENT_HISTS}
        self.n_events = 0
        self.n_tracks = 0
        self._pending = {spec[0]: [] for spec in TRACK_HISTS}
        self._pending.update({("evt", spec[0]): [] for spec in EVENT_HISTS})

    def fill(self, event, selected, weight: float = 1.0) -> None:
        if not selected:
            return
        self.n_events += 1
        self.n_tracks += len(selected)
        pending = self._pending
        for name, _, _, _, _, value in TRACK_HISTS:
            pending[name].extend(value(t) for t in selected)
        for name, _, _, _, _, value in EVENT_HISTS:
            pending[("evt", name)].append(value(event, selected))

    def flush(self, weight: float = 1.0) -> None:
        """Bin everything accumulated so far. Called once per partition."""
        for name, _, nbins, lo, hi, _ in TRACK_HISTS:
            _accumulate(self.track[name], self._pending[name], nbins, lo, hi, weight)
            self._pending[name] = []
        for name, _, nbins, lo, hi, _ in EVENT_HISTS:
            key = ("evt", name)
            _accumulate(self.event[name], self._pending[key], nbins, lo, hi, weight)
            self._pending[key] = []

    def __iadd__(self, other: "Book") -> "Book":
        for name in self.track:
            self.track[name] += other.track[name]
        for name in self.event:
            self.event[name] += other.event[name]
        self.n_events += other.n_events
        self.n_tracks += other.n_tracks
        return self


def _accumulate(counts, values, nbins, lo, hi, weight):
    """Bin values into an existing count array, clipping to the range like ROOT does."""
    if not values:
        return
    arr = np.asarray(values, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return
    idx = ((arr - lo) / (hi - lo) * nbins).astype(np.int64)
    idx = idx[(idx >= 0) & (idx < nbins)]        # drop under/overflow, as the plots do
    if idx.size:
        np.add.at(counts, idx, weight)


class Books:
    """All selection sets' books, keyed by set ID -- the parallel of Mu2eEvtAna's arrays."""

    __slots__ = ("books",)

    def __init__(self, set_ids):
        self.books = {set_id: Book(set_id) for set_id in set_ids}

    def __getitem__(self, set_id: int) -> Book:
        return self.books[set_id]

    def flush(self, weight: float = 1.0) -> None:
        for book in self.books.values():
            book.flush(weight)

    def __iadd__(self, other: "Books") -> "Books":
        for set_id, book in other.books.items():
            if set_id in self.books:
                self.books[set_id] += book
            else:
                self.books[set_id] = book
        return self

    def write(self, path: str, norm: "Norm" = None, dir_name: str = "Ana") -> None:
        """Write the books and the normalization tree in Mu2eEvtAna's layout.

        ``Ana/Hist/trk_<id>``, ``Ana/Hist/evt_<id>`` and ``Ana/data/Norm``, matching
        ``Mu2eEvtAna::InitializeOutput`` -- so ``make_plots.C``, which builds
        ``<dir>Hist/<type>_<set>/<hist>``, finds them where it expects.
        """
        import ROOT

        ROOT.gROOT.SetBatch(True)
        out = ROOT.TFile(path, "RECREATE")
        top = out.mkdir(dir_name)
        hist_dir = top.mkdir("Hist")
        data_dir = top.mkdir("data")

        for set_id in sorted(self.books):
            book = self.books[set_id]
            _write_dir(hist_dir, f"trk_{set_id}", TRACK_HISTS, book.track)
            _write_dir(hist_dir, f"evt_{set_id}", EVENT_HISTS, book.event)

        data_dir.cd()
        _write_norm(norm or Norm())
        out.Write()
        out.Close()


def _write_norm(norm: "Norm") -> None:
    """The Norm tree, with Mu2eEvtAna's branch names and Long64_t types."""
    import array

    import ROOT

    tree = ROOT.TTree("Norm", "Normalization information")
    buffers = {}
    for name in ("ngen", "nntuple", "nseen", "naccept", "nneg"):
        buffers[name] = array.array("l", [int(getattr(norm, name))])
        tree.Branch(name, buffers[name], f"{name}/L")
    tree.Fill()
    tree.Write()


def _write_dir(parent, folder, specs, counts):
    import ROOT

    subdir = parent.mkdir(folder)
    subdir.cd()
    for name, title, nbins, lo, hi, _ in specs:
        hist = ROOT.TH1F(name, f"{folder}: {title}", nbins, lo, hi)
        values = counts[name]
        for i in range(nbins):
            if values[i]:
                hist.SetBinContent(i + 1, values[i])
                hist.SetBinError(i + 1, values[i] ** 0.5)
        hist.SetEntries(float(values.sum()))
        hist.Write()
    parent.cd()
