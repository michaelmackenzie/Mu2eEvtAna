"""Restrict a written model to narrower observable ranges, keeping its shapes: e.g. to scan the fit range.

    from mu2ecombine.narrow import narrow_model
    narrow_model("output/search", "output/search_102_105", {"sr": (102., 105.)})

The model is not refit. Each narrowed channel keeps its full-range shapes:
  - templates (and their variations, and the data) keep the bins inside the range, so inside the range the narrowed
    model is the full-range model bin by bin;
  - PDFs (fitted functions, user PDFs, an envelope's RooMultiPdf with all its functions) are renormalized over the
    range by Combine, with their yields scaled to the fraction inside the range: the card rate, or the value of a
    freely floating <pdf>_norm (an envelope). An envelope's fraction is the one of its current (default) function.
Rate uncertainties and MC statistical uncertainties stay as they were evaluated over the full range. Channels not
listed are copied unchanged. Range edges are moved to the nearest bin edges of the channel.
"""

import json
import os
import shutil

import numpy as np

from .hist import make


class NarrowError(ValueError):
    pass


def _read_card(path):
    with open(path) as f:
        return f.read().splitlines()


def _card_tables(lines):
    """Indices of the observation block (bin, observation) and of the process block (bin, process, process, rate)."""
    rows = [i for i, line in enumerate(lines) if line.split()[:1] in (["bin"], ["observation"], ["rate"])
            and line.split()[1:2] != ["group"]]   # not the 'rate group = ...' nuisance group
    obs = next(i for i in rows if lines[i].split()[0] == "observation")
    rate = next(i for i in rows if lines[i].split()[0] == "rate")
    return obs - 1, obs, rate - 3, rate


def _shapes(lines):
    """{(process, channel): (file, object)} of the card's shapes lines."""
    out = {}
    for line in lines:
        t = line.split()
        if len(t) >= 5 and t[0] == "shapes":
            out[(t[1], t[2])] = (t[3], t[4])
    return out


def _shape_of(shapes, process, channel):
    for key in ((process, channel), ("*", channel), (process, "*"), ("*", "*")):
        if key in shapes:
            file, obj = shapes[key]
            return file, obj.replace("$PROCESS", process).replace("$CHANNEL", channel)
    raise NarrowError(f"no shapes line for {channel}/{process}")


def _format_row(label, values, head, width):
    return label.ljust(head) + "".join(str(v).ljust(width) for v in values)


def _subrange(edges, lo, hi, channel):
    ilo, ihi = int(np.argmin(np.abs(edges - lo))), int(np.argmin(np.abs(edges - hi)))
    if ihi <= ilo:
        raise NarrowError(f"channel {channel}: the range {lo:g}-{hi:g} holds no bin")
    return ilo, ihi


def _datahist_arrays(dh):
    n = dh.numEntries()
    values, variances = np.zeros(n), np.zeros(n)
    for i in range(n):
        dh.get(i)
        values[i] = dh.weight()
        variances[i] = dh.weightSquared(i)
    return values, variances


def _narrow_workspace(src, dst, channel, info, ilo, ihi, procs):
    """Write the narrowed workspace of one channel; returns {process: card rate factor}, {process: fraction of the
    yield inside the range}, and the data integral."""
    import ROOT
    from .combine import load_library
    from .pdfs import bin_fractions, quiet_roofit
    load_library()
    quiet_roofit()
    f = ROOT.TFile.Open(src)
    w = f.Get("w")
    obs_name = info["observable"]["name"]
    obs = w.var(obs_name)
    edges = np.array(info["observable"]["edges"])
    sub = edges[ilo:ihi + 1]
    w2 = ROOT.RooWorkspace("w", "w")
    imp = getattr(w2, "import")
    obs2 = ROOT.RooRealVar(obs_name, obs.GetTitle(), float(sub[0]), float(sub[-1]))
    widths = np.diff(sub)
    if np.allclose(widths, widths[0]):
        obs2.setBins(len(widths))
    else:
        obs2.setBinning(ROOT.RooBinning(len(widths), np.ascontiguousarray(sub)))
    imp(obs2, ROOT.RooFit.Silence())

    factors, fractions, data_integral = {}, {}, None
    for dh in w.allData():
        values, variances = _datahist_arrays(dh)
        if len(values) != len(edges) - 1:
            raise NarrowError(f"{channel}: data set {dh.GetName()} has {len(values)} bins, the channel "
                              f"{len(edges) - 1}")
        h = make(sub, values[ilo:ihi], variances[ilo:ihi], dh.GetName())
        imp(ROOT.RooDataHist(dh.GetName(), dh.GetTitle(), ROOT.RooArgList(w2.var(obs_name)), h), ROOT.RooFit.Silence())
        if dh.GetName() == f"{channel}_data_obs":
            data_integral = float(values[ilo:ihi].sum())
    for p, name in procs.items():
        if w.data(name):
            values, _ = _datahist_arrays(w.data(name))
            factors[p] = fractions[p] = values[ilo:ihi].sum() / values.sum() if values.sum() > 0 else 0.
            continue
        pdf = w.pdf(name)
        if not pdf:
            raise NarrowError(f"{channel}/{p}: no data set or PDF '{name}' in the workspace")
        current = pdf.getCurrentPdf() if pdf.InheritsFrom("RooMultiPdf") else pdf
        fraction = fractions[p] = float(bin_fractions(current, obs, edges)[ilo:ihi].sum())
        imp(pdf, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
        norm = w.arg(f"{name}_norm")
        if norm and norm.InheritsFrom("RooRealVar") and not norm.isConstant():
            # a freely floating yield (an envelope): it starts at the yield inside the range
            imp(norm, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
            n2 = w2.var(f"{name}_norm")
            n2.setVal(n2.getVal() * fraction)
            factors[p] = 1.
        else:
            if norm:
                imp(norm, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
            factors[p] = fraction
    w2.writeToFile(dst)
    f.Close()
    return factors, fractions, data_integral


def _narrow_th1_file(src, dst, ilo, ihi):
    """Narrow the TH1 shapes of a binned channel; returns {name: (narrow integral, full integral)}."""
    import ROOT
    from .hist import arrays
    f = ROOT.TFile.Open(src)
    out = ROOT.TFile(dst, "RECREATE")
    sums = {}
    for key in f.GetListOfKeys():
        h = f.Get(key.GetName())
        if not h.InheritsFrom("TH1"):
            continue
        edges, values, variances = arrays(h)
        n = make(edges[ilo:ihi + 1], values[ilo:ihi], variances[ilo:ihi], key.GetName())
        out.WriteTObject(n, key.GetName())
        sums[key.GetName()] = (float(values[ilo:ihi].sum()), float(values.sum()))
    out.Close()
    f.Close()
    return sums


def narrow_model(model_dir, out_dir, ranges):
    """Write the model in model_dir to out_dir with the channels in ranges ({channel: (lo, hi)}) narrowed. Returns the
    {channel: (lo, hi)} ranges used (moved to bin edges)."""
    with open(os.path.join(model_dir, "model.json")) as f:
        summary = json.load(f)
    unknown = [ch for ch in ranges if ch not in summary["channels"]]
    if unknown:
        raise NarrowError(f"channels {', '.join(unknown)} are not in the model {model_dir}")
    if os.path.abspath(out_dir) == os.path.abspath(model_dir):
        raise NarrowError("the narrowed model needs its own directory")
    os.makedirs(os.path.join(out_dir, "shapes"), exist_ok=True)
    cards = sorted(f for f in os.listdir(model_dir) if f.endswith(".txt"))
    combined = _read_card(os.path.join(model_dir, "combined.txt"))
    shapes = _shapes(combined)

    used, factors, observations = {}, {}, {}
    for ch, info in summary["channels"].items():
        src = os.path.join(model_dir, "shapes", f"{ch}.root")
        dst = os.path.join(out_dir, "shapes", f"{ch}.root")
        if ch not in ranges:
            shutil.copyfile(src, dst)
            continue
        edges = np.array(info["observable"]["edges"])
        ilo, ihi = _subrange(edges, *ranges[ch], ch)
        used[ch] = (float(edges[ilo]), float(edges[ihi]))
        procs = list(info["processes"])
        if info["mode"] == "binned":
            sums = _narrow_th1_file(src, dst, ilo, ihi)
            factors[ch] = {p: (sums[p][0] / sums[p][1] if sums[p][1] > 0 else 0.) for p in procs}
            fractions = factors[ch]
            observations[ch] = sums["data_obs"][0]
        else:
            objects = {p: _shape_of(shapes, p, ch)[1].split(":", 1)[1] for p in procs}
            factors[ch], fractions, observations[ch] = _narrow_workspace(src, dst, ch, info, ilo, ihi, objects)
        info["observable"].update(range=list(used[ch]), bins=ihi - ilo, edges=[float(e) for e in edges[ilo:ihi + 1]])
        info["data"]["integral"] = observations[ch]
        for p in procs:
            info["processes"][p]["rate"] *= fractions[p]
        info["narrowed"] = {"from": os.path.abspath(model_dir), "range": list(used[ch])}

    for card in cards:
        lines = _read_card(os.path.join(model_dir, card))
        try:
            b, o, pb, r = _card_tables(lines)
        except StopIteration:
            shutil.copyfile(os.path.join(model_dir, card), os.path.join(out_dir, card))
            continue
        bins, obs_values = lines[b].split()[1:], lines[o].split()[1:]
        obs_values = [f"{observations[c]:.6g}" if c in observations else v for c, v in zip(bins, obs_values)]
        pbins, pnames = lines[pb].split()[1:], lines[pb + 1].split()[1:]
        rates = lines[r].split()[1:]
        new_rates = []
        for c, p, rate in zip(pbins, pnames, rates):
            new_rates.append(f"{float(rate) * factors[c][p]:.6g}" if c in factors else rate)
        width = max([len(x) for x in bins + pnames] + [12]) + 2
        lines[o] = _format_row("observation", obs_values, 16, width)
        lines[r] = _format_row("rate", new_rates, 28, width)
        with open(os.path.join(out_dir, card), "w") as f:
            f.write("\n".join(lines) + "\n")
    for extra in ("commands.sh",):
        if os.path.exists(os.path.join(model_dir, extra)):
            shutil.copyfile(os.path.join(model_dir, extra), os.path.join(out_dir, extra))
    with open(os.path.join(out_dir, "model.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return used
