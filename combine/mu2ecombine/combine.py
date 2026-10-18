"""Small helpers around a Combine installation."""

COMBINE_LIBRARY = "libHiggsAnalysisCombinedLimit.so"


class CombineError(RuntimeError):
    pass


_loaded = False


def load_library():
    """Load the Combine library into ROOT (needed for RooMultiPdf envelopes)."""
    global _loaded
    if _loaded:
        return
    import ROOT
    if ROOT.gSystem.Load(COMBINE_LIBRARY) < 0:
        raise CombineError(f"Unable to load {COMBINE_LIBRARY}: set up Combine first (see README.md)")
    _loaded = True


def read_limits(path):
    """Rows of a Combine output's 'limit' tree as dicts (e.g. higgsCombine.<name>.AsymptoticLimits.mH120.root)."""
    import ROOT
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise CombineError(f"Unable to open {path}")
    t = f.Get("limit")
    if not t:
        raise CombineError(f"No limit tree in {path}")
    names = [b.GetName() for b in t.GetListOfBranches()]
    rows = [{n: getattr(entry, n) for n in names} for entry in t]
    f.Close()
    return rows


def asymptotic_limits(path):
    """{quantile: limit} from an AsymptoticLimits output (quantile -1 is the observed limit)."""
    return {round(float(r["quantileExpected"]), 3): float(r["limit"]) for r in read_limits(path)}
