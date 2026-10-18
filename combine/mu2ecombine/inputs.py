"""Reading Mu2eEvtAna histogram files and normalizing MC samples. Histograms are returned as detached TH1D copies.

    f = Mu2eEvtAnaFile("ConvAna.cnv_ana.dio00b1s5r0101.m1.root")
    h = f.hist("trk_40/obs")                       # Ana/Hist/trk_40/obs
    up, down = f.variation("Scale", set=40)        # Ana/Hist/sys_40/obs_1, obs_2

    dio = Sample(f, weight=npot * rate / ngen, expected_events=9368976)
    h = dio.hist("trk_40/obs")                     # normalized to expected events

Histograms from any other file can be read with ROOT directly (ROOT.TFile(path).Get(name)).
"""

import os
import re
import warnings

from .hist import copy, is_th1


class InputError(RuntimeError):
    pass


_open_files = {}


def open_root_file(path):
    import ROOT
    f = _open_files.get(path)
    if f is None or not f.IsOpen():
        if "://" not in path and not os.path.exists(path):
            raise InputError(f"Input file not found: {path}")
        f = ROOT.TFile.Open(path, "READ")
        if not f or f.IsZombie():
            raise InputError(f"Unable to open {path}")
        _open_files[path] = f
    return f


#-----------------------------------------------------------------------------------------------------------------------
# Mu2eEvtAna systematic variations: names and histogram indices from Mu2eEvtAna/inc/Systematics.hh

def systematics_header():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "..", "inc", "Systematics.hh"))


_systematics = None


def mu2eevtana_systematics():
    """Mapping name -> (up index, down index) of the Mu2eEvtAna systematic histograms (sys_<set>/obs_<index>), read
    from the Sys_t enum of Systematics.hh (k<Name>Up = <index>, k<Name>Down = <index>)."""
    global _systematics
    if _systematics is None:
        path = systematics_header()
        if not os.path.exists(path):
            raise InputError(f"Mu2eEvtAna systematics header not found at {path}")
        with open(path) as f:
            text = f.read()
        found = {}
        for name, side, index in re.findall(r"\bk(\w+?)(Up|Down)\s*=\s*(\d+)", text):
            found.setdefault(name, {})[side] = int(index)
        _systematics = {name: (d["Up"], d["Down"]) for name, d in found.items() if "Up" in d and "Down" in d}
    return _systematics


class Mu2eEvtAnaFile:
    """A Mu2eEvtAna histogram file: histograms under <top>/Hist, the normalization tree <top>/data/Norm.

    hist("trk_40/obs") reads <top>/Hist/trk_40/obs; a path starting with '/' is taken from the file's top level.
    """

    def __init__(self, path, top="Ana"):
        self.path = str(path)
        self.top = top

    def hist(self, path):
        full = path[1:] if path.startswith("/") else f"{self.top}/Hist/{path}"
        h = open_root_file(self.path).Get(full)
        if not h:
            raise InputError(f"'{full}' not found in {self.path}")
        if not is_th1(h):
            raise InputError(f"'{full}' in {self.path} is not a 1D histogram")
        return copy(h, name=path.replace("/", "_"))

    def variation(self, name, set, hist="obs"):
        """(up, down) histograms of the Mu2eEvtAna systematic `name` (e.g. 'Scale') for selection `set`."""
        known = mu2eevtana_systematics()
        if name not in known:
            raise InputError(f"Unknown Mu2eEvtAna systematic '{name}' (defined: {', '.join(sorted(known))})")
        up, down = known[name]
        return self.hist(f"sys_{set}/{hist}_{up}"), self.hist(f"sys_{set}/{hist}_{down}")

    def events_seen(self):
        """Total events the analyzer processed (sum of nseen over the Norm tree entries), or None without the tree."""
        tree = open_root_file(self.path).Get(f"{self.top}/data/Norm")
        if not tree:
            return None
        return sum(int(entry.nseen) for entry in tree)

    def __repr__(self):
        return f"Mu2eEvtAnaFile('{self.path}')"


class Sample:
    """A Mu2eEvtAna file with a per-event normalization.

    weight:          per-MC-event weight (e.g. exposure x rate / N(generated)); 1 for data
    expected_events: events in the full dataset; when given, histograms are also scaled by expected / events seen
                     (from the file's Norm tree), correcting for a partially processed dataset
    file:            a path or a Mu2eEvtAnaFile
    """

    def __init__(self, file, weight=1., expected_events=None, name=None):
        self.file = file if isinstance(file, Mu2eEvtAnaFile) else Mu2eEvtAnaFile(file)
        self.name = name or os.path.basename(self.file.path)
        self.base_weight = float(weight)
        self.expected_events = expected_events
        self._completeness = None

    @property
    def completeness(self):
        if self.expected_events is None:
            return 1.
        if self._completeness is None:
            seen = self.file.events_seen()
            if seen is None:
                raise InputError(f"{self.name}: 'expected_events' is set but {self.file.path} has no Norm tree")
            if seen <= 0:
                raise InputError(f"{self.name}: the Norm tree reports no events seen")
            self._completeness = float(self.expected_events) / seen
            if self._completeness < 1. - 1e-9:
                warnings.warn(f"{self.name}: {seen} events seen but {self.expected_events:.0f} expected; the input "
                              "processed more events than the dataset has, check expected_events")
        return self._completeness

    @property
    def weight(self):
        """Weight applied to the histograms: weight x completeness."""
        return self.base_weight * self.completeness

    def hist(self, path):
        h = self.file.hist(path)
        h.Scale(self.weight)
        return h

    def variation(self, name, set, hist="obs"):
        """Normalized (up, down) Mu2eEvtAna systematic histograms (see Mu2eEvtAnaFile.variation)."""
        up, down = self.file.variation(name, set, hist)
        up.Scale(self.weight)
        down.Scale(self.weight)
        return up, down

    def __repr__(self):
        return f"Sample('{self.name}', weight={self.base_weight:.4g})"
