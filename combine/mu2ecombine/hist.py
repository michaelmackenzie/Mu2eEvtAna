"""Operations on ROOT histograms used to build templates.

Every function takes TH1s and returns a new TH1D that is not attached to any file (SetDirectory(0)), with Sumw2 errors,
so MC statistical uncertainties follow every step. Inputs are never modified.

    h = rebin_to(h, edges)              # merge bins into edges that must all be edges of h
    up = transformed(h, shift=0.1)      # x -> x + 0.1 (or a stretch about the peak), contents redistributed
    t = convolved(truth, response)      # truth spectrum smeared by a response histogram (x_reco - x_true)
    h = replace_tail(h, 102., t)        # bins above 102 take the shape of t, keeping h's content there
    h = total([h1, h2])                 # sum
"""

import itertools

import numpy as np


class HistError(ValueError):
    pass


_counter = itertools.count()


def _unique(name):
    return f"{name or 'h'}_m2c{next(_counter)}"


def is_th1(obj):
    try:
        return bool(obj.InheritsFrom("TH1")) and obj.GetDimension() == 1
    except AttributeError:
        return False


def arrays(h):
    """(edges, contents, variances) of a 1D histogram as numpy arrays (no underflow/overflow)."""
    if not is_th1(h):
        raise HistError(f"expected a 1D ROOT histogram, got {type(h).__name__}")
    n = h.GetNbinsX()
    axis = h.GetXaxis()
    edges = np.array([axis.GetBinLowEdge(i) for i in range(1, n + 2)])
    values = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    variances = np.array([h.GetBinError(i) for i in range(1, n + 1)]) ** 2
    return edges, values, variances


def make(edges, values, variances=None, name="h", title=""):
    """A detached TH1D from bin edges, contents, and variances (default: no errors)."""
    import ROOT
    edges = np.ascontiguousarray(edges, dtype=float)
    h = ROOT.TH1D(_unique(name), title, len(edges) - 1, edges)
    h.SetDirectory(0)
    h.Sumw2()
    variances = np.zeros(len(values)) if variances is None else variances
    for i, (v, var) in enumerate(zip(values, variances), start=1):
        h.SetBinContent(i, float(v))
        h.SetBinError(i, float(np.sqrt(max(var, 0.))))
    return h


def copy(h, name=None):
    """A detached TH1D copy (also converts TH1F etc.)."""
    edges, values, variances = arrays(h)
    return make(edges, values, variances, name or h.GetName(), h.GetTitle())


def total(hists, name=None):
    """Sum of TH1s with the same binning, as a new TH1D."""
    hists = [hists] if is_th1(hists) else list(hists)
    if not hists:
        raise HistError("total: no histograms given")
    edges, values, variances = arrays(hists[0])
    for h in hists[1:]:
        e, v, var = arrays(h)
        if len(e) != len(edges) or not np.allclose(e, edges, rtol=0, atol=1e-9):
            raise HistError(f"total: '{h.GetName()}' and '{hists[0].GetName()}' have different binnings")
        values, variances = values + v, variances + var
    return make(edges, values, variances, name or hists[0].GetName())


def as_th1(h, what="histogram"):
    """A detached TH1D from a TH1 or a list of TH1s (summed)."""
    if is_th1(h):
        return copy(h)
    if isinstance(h, (list, tuple)) and h and all(is_th1(x) for x in h):
        return total(h)
    raise HistError(f"{what}: expected a 1D ROOT histogram or a list of them, got {type(h).__name__}")


def integral(h):
    """Sum of the bin contents (no underflow/overflow)."""
    return float(arrays(h)[1].sum())


def rel_error(h):
    """Relative MC statistical uncertainty of the integral."""
    _, v, var = arrays(h)
    return float(np.sqrt(var.sum()) / v.sum()) if v.sum() > 0 else 0.


def _edge_index(edges, value, tolerance=1e-6):
    match = np.where(np.abs(edges - value) <= tolerance * max(1., abs(value)))[0]
    return int(match[0]) if len(match) else None


def rebin_to(h, new_edges, name=None):
    """Merge bins into new_edges, which must all be edges of h. Contents outside the new range are dropped."""
    e, values, variances = arrays(h)
    new_edges = np.asarray(new_edges, dtype=float)
    index = []
    for x in new_edges:
        i = _edge_index(e, x)
        if i is None:
            raise HistError(f"'{h.GetName()}': bin edge {x:g} is not an edge of its binning "
                            f"({e[0]:g}-{e[-1]:g}, {len(e) - 1} bins)")
        index.append(i)
    return make(new_edges, [values[a:b].sum() for a, b in zip(index, index[1:])],
                [variances[a:b].sum() for a, b in zip(index, index[1:])], name or h.GetName())


def clip_negative(h):
    """Set negative bins to zero (in place). Returns h."""
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(i) < 0:
            h.SetBinContent(i, 0.)
    return h


def transformed(h, shift=0., stretch=0., about=None):
    """Template for x -> about + (x - about) * (1 + stretch) + shift, each bin's content spread uniformly over its
    transformed interval (an energy scale or resolution variation). `about` defaults to the peak bin."""
    e, values, variances = arrays(h)
    centers = 0.5 * (e[1:] + e[:-1])
    about = centers[np.argmax(values)] if about is None else about
    lo = about + (e[:-1] - about) * (1. + stretch) + shift
    hi = about + (e[1:] - about) * (1. + stretch) + shift
    out_v, out_var = np.zeros(len(values)), np.zeros(len(values))
    for i in range(len(values)):
        width = hi[i] - lo[i]
        if width <= 0:
            continue
        frac = np.clip(np.minimum(e[1:], hi[i]) - np.maximum(e[:-1], lo[i]), 0., None) / width
        out_v += frac * values[i]
        out_var += frac * variances[i]
    return make(e, out_v, out_var, h.GetName())


def convolved(h, response, new_edges=None, step=None):
    """Distribution of x + d, x from h and d from `response` (a histogram of x_reco - x_true), binned in new_edges
    (default: h's binning) and normalized to h's integral. The result has no errors."""
    e, values, _ = arrays(h)
    re, rv, _ = arrays(response)
    new_edges = e if new_edges is None else np.asarray(new_edges, dtype=float)
    step = step or min(np.diff(e).min(), np.diff(re).min()) / 4.
    xt = np.arange(e[0] + step / 2, e[-1], step)
    pt = np.interp(xt, 0.5 * (e[1:] + e[:-1]), values / np.diff(e), left=0., right=0.)
    d = np.arange(re[0] + step / 2, re[-1], step)
    pd = np.clip(np.interp(d, 0.5 * (re[1:] + re[:-1]), rv / np.diff(re), left=0., right=0.), 0., None)
    out, _ = np.histogram((xt[:, None] + d[None, :]).ravel(), bins=new_edges,
                          weights=(pt[:, None] * pd[None, :]).ravel())
    if out.sum() <= 0:
        raise HistError(f"'{h.GetName()}': the convolution is empty in the requested range")
    return make(new_edges, out * values.sum() / out.sum(), None, h.GetName())


def replace_tail(h, start, reference, norm_from=None, rel_error=None):
    """Copy of h with the bins above `start` replaced by the shape of `reference` (same binning), normalized to h's
    content above `norm_from` (default `start`). The replaced bins get errors of rel_error x content, or else h's
    relative uncertainty above norm_from."""
    e, values, variances = arrays(h)
    re, rv, _ = arrays(reference)
    if len(re) != len(e) or not np.allclose(re, e, rtol=0, atol=1e-9):
        raise HistError(f"replace_tail: '{reference.GetName()}' and '{h.GetName()}' have different binnings")
    centers = 0.5 * (e[1:] + e[:-1])
    norm_from = start if norm_from is None else norm_from
    above = centers > norm_from
    if rv[above].sum() <= 0:
        raise HistError(f"'{reference.GetName()}': no content above {norm_from:g}")
    content = values[above].sum()
    tail = centers > start
    values, variances = values.copy(), variances.copy()
    values[tail] = rv[tail] * content / rv[above].sum()
    if rel_error is None:
        rel_error = np.sqrt(variances[above].sum()) / max(content, 1e-300)
    variances[tail] = (rel_error * values[tail]) ** 2
    return make(e, values, variances, h.GetName())
