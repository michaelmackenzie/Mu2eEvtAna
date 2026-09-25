"""Selection sets, keyed by the set ID that names the histogram directory.

``Mu2eEvtAna.cc`` books a histogram set per selection and writes it to ``trk_<id>`` /
``evt_<id>``; this does the same, with the set ID as the dictionary key. Adding a set is
one entry in :data:`SELECTIONS`.

Two sets to begin with:

``75``  the signal selection -- the same cuts ``Mu2eEvtAna`` applies, validated cut-flow
        stage by stage against pyfitter's implementation.
``0``   a basic ``p > 90 MeV/c`` set, for a look at the spectrum before the ID cuts.
"""

from __future__ import annotations

from statistics import fmean

from pyevtana import is_missing
from pyevtana.surfaces import surface_id

TT_FRONT = surface_id("TT_Front")
TT_MID = surface_id("TT_Mid")
ST_FOILS = surface_id("ST_Foils")
OPA = surface_id("OPA")
ST_BOUNDARY = {surface_id(n) for n in ("ST_Front", "ST_Back", "ST_Inner", "ST_Outer")}

#: branches every selection here needs
BRANCHES = ["evtinfo", "trk", "trkqual", "trkpid", "trksegs", "trksegpars_lh",
            "trkcalohit", "crvcoincs", "caloclusters"]

TRIGGERS = ("cpr_TrkDe_80m70p", "apr_TrkDe_80m70p")
MOM_LO, MOM_HI = 100.0, 115.0
T0_LO, T0_HI = 475.0, 1650.0


class Track:
    """One track, with everything the selections and histograms need read once.

    Reading a field off a proxy is cheap but not free, and several of these are needed by
    both the cuts and the histograms, so they are gathered in a single pass here.
    """

    __slots__ = ("index", "pdg", "nactive", "status", "goodfit", "trkqual", "trkpid",
                 "front_times", "front_p", "front_pt", "front_cos", "front_pz",
                 "mid_t0err", "tandip", "d0", "rmax", "t0", "cluster_e",
                 "n_st", "n_opa", "st_boundary", "is_good", "is_downstream",
                 "first_front_t", "mean_front_t")

    def __init__(self, track):
        self.index = track.index
        self.pdg = int(track.pdg)
        self.nactive = int(track.nactive)
        self.status = int(track.status)
        self.goodfit = int(track.goodfit)
        self.trkqual = float(track.qual().result)
        self.trkpid = float(track.pid().result)

        segs = track.segs()
        sids = segs._column("sid")
        times = segs._column("time")
        moms = segs._column("mom")

        # one pass over the segments: which are at the surfaces of interest, and the
        # counts. Five separate comprehensions over the same ~17 entries was most of the
        # per-track cost.
        front, mid = [], []
        n_st = n_opa = 0
        st_boundary = False
        for i, sid in enumerate(sids):
            if sid == TT_FRONT:
                front.append(i)
            elif sid == TT_MID:
                mid.append(i)
            if sid == ST_FOILS:
                n_st += 1
            elif sid == OPA:
                n_opa += 1
            if sid in ST_BOUNDARY:
                st_boundary = True
        self.n_st, self.n_opa, self.st_boundary = n_st, n_opa, st_boundary

        # build each front-segment vector once; .mag/.pt/.z off the same object
        front_vecs = [moms[i] for i in front]
        self.front_times = [times[i] for i in front]
        self.front_pz = [v.z for v in front_vecs]
        self.front_p = [v.mag for v in front_vecs]
        if front_vecs:
            lead = front_vecs[0]
            self.front_pt = lead.pt
            self.front_cos = lead.z / self.front_p[0] if self.front_p[0] > 0 else -2.0
        else:
            self.front_pt, self.front_cos = -1.0, -2.0

        pars = track.segpars()
        if is_missing(pars) or not front:
            self.mid_t0err, self.tandip, self.d0, self.rmax = [], -100.0, -1e4, -1.0
        else:
            first = front[0]
            t0err_col = pars._column("t0err")
            self.mid_t0err = [t0err_col[i] for i in mid]
            self.tandip = pars._column("tanDip")[first]
            self.d0 = pars._column("d0")[first]
            self.rmax = pars._column("maxr")[first]

        calohit = track.calohit()
        self.cluster_e = float(calohit.edep) if int(calohit.did) >= 0 else -1.0

        if front:
            earliest = min(self.front_times)
            self.is_downstream = any(pz > 0 for t, pz in zip(self.front_times, self.front_pz)
                                     if t == earliest)
            self.first_front_t = self.front_times[0]
            self.mean_front_t = fmean(self.front_times)
            self.t0 = self.front_times[0]
        else:
            self.is_downstream = False
            self.first_front_t = self.mean_front_t = self.t0 = 0.0

        self.is_good = self.status >= 0 and self.goodfit != 0

    @property
    def p(self) -> float:
        """Momentum at the tracker front, as ``Track_t::PFront()`` uses."""
        return self.front_p[0] if self.front_p else -1.0

    @property
    def ep(self) -> float:
        return self.cluster_e / self.p if self.cluster_e > 0 and self.p > 0 else -1.0


class Event:
    """Per-event quantities the cross-track and veto cuts need."""

    __slots__ = ("tracks", "trigger", "crv_times", "has_calo_energy", "run", "subrun",
                 "event")

    def __init__(self, event):
        self.tracks = [Track(t) for t in event.Tracks()]
        self.trigger = any(event.trigger(name) for name in TRIGGERS)
        self.crv_times = ([float(c.time) for c in event.CrvCoincs()]
                          if event.has("crvcoincs") else [])
        self.has_calo_energy = (any(c.energyDep_ > 0 for c in event.CaloClusters())
                                if event.has("caloclusters") else False)
        info = event.info
        self.run, self.subrun, self.event = int(info.run), int(info.subrun), int(info.event)


# --------------------------------------------------------------------------------------
# the cuts
# --------------------------------------------------------------------------------------


def _no_upstream_partner(t, e):
    if not (t.is_downstream and t.is_good):
        return True
    for other in e.tracks:
        if other.index == t.index or other.is_downstream or not other.is_good:
            continue
        if 40.0 <= t.first_front_t - other.first_front_t <= 110.0:
            return False
    return True


def _no_multi_track(t, e):
    if not (t.is_downstream and abs(t.pdg) == 11 and t.is_good):
        return True
    for other in e.tracks:
        if other.index == t.index:
            continue
        if not (other.is_downstream and abs(other.pdg) == 11 and other.is_good):
            continue
        if abs(t.mean_front_t - other.mean_front_t) < 150.0:
            return False
    return True


def _no_crv(t, e):
    return not any(0 < time - crv < 150
                   for time in t.front_times for crv in e.crv_times)


#: the signal selection, in the order Mu2eEvtAna applies it. Each entry is
#: ``(name, predicate)``; the cut flow is reported per event and cumulative.
SIGNAL_CUTS = [
    ("has_a_track",                  lambda t, e: bool(e.tracks)),
    ("is_good_track",                lambda t, e: t.status >= 0 and t.goodfit != 0),
    ("has_trk_front_seg",            lambda t, e: bool(t.front_times)),
    ("is_reco_electron_or_positron", lambda t, e: abs(t.pdg) == 11),
    ("has_downstream",               lambda t, e: t.is_downstream),
    ("charge_selection",             lambda t, e: t.pdg == 11),
    ("or_trigger",                   lambda t, e: e.trigger),
    ("upstream_veto",                _no_upstream_partner),
    ("no_multi_trk_veto",            _no_multi_track),
    ("good_trkpid",                  lambda t, e: t.trkpid > 0.54 and e.has_calo_energy),
    ("pz_over_pt",                   lambda t, e: 0.575 < t.tandip < 0.85),
    ("st_boundary",                  lambda t, e: t.st_boundary),
    ("has_st",                       lambda t, e: t.n_st > 0),
    ("no_opa",                       lambda t, e: t.n_opa == 0),
    ("good_trkqual",                 lambda t, e: t.trkqual > 0.155),
    ("has_hits",                     lambda t, e: t.nactive >= 20),
    ("within_t0err",                 lambda t, e: all(v < 0.85 for v in t.mid_t0err)),
    ("no_crv_veto",                  _no_crv),
    ("in_mom_range",                 lambda t, e: all(MOM_LO < p < MOM_HI for p in t.front_p)),
    ("within_t0_475",                lambda t, e: all(T0_LO < v < T0_HI
                                                      for v in t.front_times)),
]

#: a deliberately loose set: a downstream electron-hypothesis track above 90 MeV/c
BASIC_CUTS = [
    ("has_a_track",      lambda t, e: bool(e.tracks)),
    ("is_good_track",    lambda t, e: t.status >= 0 and t.goodfit != 0),
    ("is_reco_electron", lambda t, e: t.pdg == 11),
    ("has_downstream",   lambda t, e: t.is_downstream),
    ("p_gt_90",          lambda t, e: t.p > 90.0),
]


class SelectionSet:
    """A named, numbered list of cuts -- the key under which histograms are booked."""

    __slots__ = ("set_id", "name", "cuts")

    def __init__(self, set_id: int, name: str, cuts: list):
        self.set_id = set_id
        self.name = name
        self.cuts = cuts

    @property
    def stages(self) -> list:
        return [name for name, _ in self.cuts]

    def select(self, event: Event) -> tuple:
        """Return (surviving tracks, per-stage survivor counts) for one event."""
        surviving = event.tracks
        per_stage = []
        for _, cut in self.cuts:
            surviving = [t for t in surviving if cut(t, event)]
            per_stage.append(len(surviving))
        return surviving, per_stage


#: set ID -> selection. The ID is what names the output directories.
SELECTIONS = {
    75: SelectionSet(75, "signal selection", SIGNAL_CUTS),
    0:  SelectionSet(0, "p > 90 MeV/c", BASIC_CUTS),
}
