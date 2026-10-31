"""Model figures: the expected (pre-fit) stack of each channel and the MC shape fits, and post-fit figures from Combine.

    model.write(outdir, plots=True)                                  # figures in <outdir>/figures
    python3 -m mu2ecombine.plots fitDiagnostics.<name>.root          # in <outdir>, after FitDiagnostics --saveShapes

Stacks are drawn as cumulative histograms rather than a THStack.
"""

import argparse
import json
import os
import re

import numpy as np

COLORS = [ "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581",
           "#92dadd"]


def _style():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetErrorX(0.)        # data points: vertical error bars only, without end caps
    ROOT.gStyle.SetEndErrorSize(0.)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning


def plot_model(results, figdir, label=None):
    """Validation figures of every channel, linear and log scale (<name>.png and <name>_log.png):
    <channel>_prefit             the expected stack (with the data unless it is Asimov), headed by label (see
                                 plot_channel)
    <channel>_inputs             each process's expected distribution over the data, with the pulls of the data
    <channel>_<process>_fit      each fit over the histogram it was fit to, with the pulls
    <channel>_<process>_template each binned template (MC) process, with its MC statistical uncertainties"""
    os.makedirs(figdir, exist_ok=True)
    _style()
    for name, r in results.items():
        for logy, suffix in ((False, ""), (True, "_log")):
            plot_channel(r, os.path.join(figdir, f"{name}_prefit{suffix}.png"), logy, label)
            plot_inputs(r, os.path.join(figdir, f"{name}_inputs{suffix}.png"), logy, label)
            for p, fit in r.fits.items():
                if p in r.processes:
                    plot_fit(r, p, fit, os.path.join(figdir, f"{name}_{p}_fit{suffix}.png"), logy, label)
            for p in r.processes:
                if is_template(r, p):
                    plot_template(r, p, os.path.join(figdir, f"{name}_{p}_template{suffix}.png"), logy, label)


def is_template(r, process):
    """Whether a process of a built channel is a binned (MC) template rather than a fit or an envelope."""
    return process in r.templates and process not in r.fits and process not in r.envelopes


def _color(r, process, index, signal_index=None):
    """The figure color of a process: its own (Process.color), else the palette (signals: red, magenta, green)."""
    import ROOT
    color = r.channel.processes[process].color
    if color is None:
        if signal_index is not None:
            return [ROOT.kRed, ROOT.kMagenta + 1, ROOT.kGreen + 2][signal_index % 3]
        color = COLORS[index % len(COLORS)]
    return ROOT.TColor.GetColor(color) if isinstance(color, str) else int(color)


FONT = 132   # the figure font (Times-like, as the stack figures)


def _axes(h, label_size, title_size, y_offset, x_offset=1.):
    """Font and sizes of a histogram's axes."""
    for axis in (h.GetXaxis(), h.GetYaxis()):
        axis.SetLabelFont(FONT)
        axis.SetTitleFont(FONT)
        axis.SetLabelSize(label_size)
        axis.SetTitleSize(title_size)
    h.GetYaxis().SetTitleOffset(y_offset)
    h.GetXaxis().SetTitleOffset(x_offset)


def _single_canvas():
    import ROOT
    c = ROOT.TCanvas("c", "", 1100, 825)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.13)
    c.SetTicks()
    return c


def _split_canvas():
    """A canvas with a main pad and a lower pad for ratios or pulls: (canvas, top, bottom). Draw the main pad's axes
    with _axes(h, 0.06, 0.065, 0.95) and its x labels hidden."""
    import ROOT
    c = ROOT.TCanvas("c", "", 1100, 1031)
    top = ROOT.TPad("top", "", 0, 0.3, 1, 1)
    bottom = ROOT.TPad("bottom", "", 0, 0, 1, 0.3)
    for pad in (top, bottom):
        pad.SetRightMargin(0.05)
        pad.SetLeftMargin(0.13)
        pad.SetTicks()
    top.SetTopMargin(0.05)
    top.SetBottomMargin(0.03)
    bottom.SetTopMargin(0.02)
    bottom.SetBottomMargin(0.3)
    top.Draw()
    bottom.Draw()
    return c, top, bottom


def _header(pad, label, keep, scale=1.):
    """The (left, right) TLatex header above a pad's frame (see plot_channel); scale: the text size relative to a
    single-pad canvas (the main pad of a split canvas is shorter)."""
    import ROOT
    if not label:
        return
    left, right = label
    text = ROOT.TLatex()
    text.SetNDC()
    y0 = 1. - pad.GetTopMargin() + 0.017 * scale
    if left:
        text.SetTextFont(42)
        text.SetTextSize(0.039 * scale)
        text.SetTextAlign(11)
        text.DrawLatex(pad.GetLeftMargin() + 0.015, y0, left)
    if right:
        text.SetTextFont(FONT)
        text.SetTextSize(0.03 * scale)
        text.SetTextAlign(31)
        text.DrawLatex(1. - pad.GetRightMargin(), y0, right)
    keep.append(text)


SPLIT_SCALE = 825. / (0.7 * 1031.)   # text in the main pad of a split canvas, relative to a single canvas


def _log_range(values, top):
    """(min, max) of a log-scale frame: from below the smallest positive value to well above top."""
    positive = values[values > 0]
    lo = 0.5 * positive.min() if positive.size else 1e-6 * top
    lo = max(lo, 1e-4 * top)
    return lo, lo * (top / lo) ** (1. / 0.8)


def _lower_frame(edges, x_title, y_title, lo, hi, keep):
    """Draw the axes of a lower pad, with a dashed line at 0."""
    import ROOT
    frame = _hist("_lower_frame", edges, np.zeros(len(edges) - 1))
    frame.SetMinimum(lo)
    frame.SetMaximum(hi)
    frame.GetXaxis().SetTitle(x_title)
    frame.GetYaxis().SetTitle(y_title)
    _axes(frame, 0.13, 0.15, 0.4, 0.9)
    frame.GetYaxis().SetTitleSize(0.1)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetNdivisions(505)
    frame.Draw("axis")
    line = ROOT.TLine(edges[0], 0., edges[-1], 0.)
    line.SetLineStyle(2)
    line.Draw()
    keep += [frame, line]


def _lower_points(name, edges, values, errors, keep):
    h = _hist(name, edges, values, errors)
    h.SetMarkerStyle(20)
    h.SetLineWidth(2)
    h.SetMarkerSize(0.8)
    h.SetLineColor(1)
    h.Draw("E same")
    keep.append(h)


def _bin_width_title(obs_title, edges):
    """'Events / <width> <unit>' for uniform bins, the unit taken from the "(unit)" or "[unit]" of the observable
    title; 'Events / bin' otherwise."""
    widths = np.diff(edges)
    if not np.allclose(widths, widths[0]):
        return "Events / bin"
    m = re.search(r"[\(\[]([^\)\]]+)[\)\]]\s*$", obs_title)
    return f"Events / {widths[0]:.2g} {m.group(1)}" if m else f"Events / {widths[0]:.2g}"


def plot_channel(r, path, logy=False, label=None):
    """The expected stack: the backgrounds stacked (smallest at the bottom), the signals as lines, and the data unless
    it is the Asimov data set (then it is the stack itself). label: (left, right) TLatex header above the frame, e.g.
    ("#bf{Experiment} #it{Simulation}", "<exposure>"); either may be None."""
    import ROOT
    from .hist import arrays
    ch = r.channel
    edges = ch.edges
    bkg = [p for p in r.processes if not ch.processes[p].signal]
    sig = [p for p in r.processes if ch.processes[p].signal]
    c = ROOT.TCanvas("c", "", 1100, 825)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.13)
    c.SetTicks()
    keep = []
    order = sorted(bkg, key=lambda p: r.expected[p].sum())
    cumulative = np.zeros(len(edges) - 1)
    stack = []
    for i, p in enumerate(order):
        cumulative = cumulative + r.expected[p]
        h = _hist(f"_stack_{p}", edges, cumulative)
        h.SetFillColor(_color(r, p, i))
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(1)
        stack.append((p, h))
    show_data = r.data_source != "asimov"
    data_vals = arrays(r.data)[1]
    data = _hist("_data", edges, data_vals, np.sqrt(np.clip(data_vals, 0., None)))
    data.SetMarkerStyle(20)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    sig_hists = []
    for i, p in enumerate(sig):
        h = _hist(f"_sig_{p}", edges, r.expected[p])
        h.SetLineColor(_color(r, p, i, signal_index=i))
        h.SetLineWidth(2)
        sig_hists.append((p, h))
    ymax = max([cumulative.max()] + ([data_vals.max()] if show_data else [])
               + [h.GetMaximum() for _, h in sig_hists] + [1e-9])
    frame = _hist("_frame", edges, np.zeros(len(edges) - 1))
    frame.SetTitle(f";{ch.obs.GetTitle()};{_bin_width_title(ch.obs.GetTitle(), edges)}")
    for axis in (frame.GetXaxis(), frame.GetYaxis()):
        axis.SetLabelFont(132)
        axis.SetTitleFont(132)
        axis.SetLabelSize(0.05)
        axis.SetTitleSize(0.05)
    if logy:   # the legend takes the top ~27% of the frame
        frame.SetMinimum(1e-4)
        frame.SetMaximum(1e-4 * (max(ymax, 1e-3) / 1e-4) ** (1. / 0.73))
    else:
        frame.SetMinimum(0.)
        frame.SetMaximum(1.35 * ymax)
    frame.Draw("axis")
    for _, h in reversed(stack):
        h.Draw("hist same")
    for _, h in sig_hists:
        h.Draw("hist same")
    if show_data:
        data.Draw("E same")
        if ch.blind:
            box = ROOT.TBox(ch.blind[0], frame.GetMinimum(), ch.blind[1], frame.GetMaximum())
            box.SetFillColorAlpha(ROOT.kGray, 0.35)
            box.Draw()
            keep.append(box)
    leg = ROOT.TLegend(0.03 + c.GetLeftMargin(), 0.75, 0.98 - c.GetRightMargin(), 0.97 - c.GetTopMargin())
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(132)
    leg.SetTextSize(0.043)
    if show_data:
        leg.AddEntry(data, "Data", "PE")
    for p, h in sig_hists:
        leg.AddEntry(h, ch.processes[p].title, "L")
    for p, h in reversed(stack):
        leg.AddEntry(h, ch.processes[p].title, "F")
    leg.Draw()
    _header(c, label, keep)
    c.SetLogy(logy)
    c.RedrawAxis()
    c.SaveAs(path)
    c.Close()


def plot_fit(r, process, fit, path, logy=False, label=None):
    """The fitted function (a smooth curve) over the histogram it was fit to, with the pulls of the histogram,
    (N - N_fit) / sigma. label: the header (see plot_channel)."""
    import ROOT
    from .hist import arrays
    from .pdfs import bin_fractions
    ch = r.channel
    edges = ch.edges
    h = fit.fit_hist(edges)
    _, values, variances = arrays(h)
    errors = np.sqrt(variances)
    centers = 0.5 * (edges[1:] + edges[:-1])
    frac = bin_fractions(fit.pdf, ch.obs, edges)
    mask = np.ones(len(centers), dtype=bool)
    if fit.fit_range is not None:
        mask = (centers > fit.fit_range[0]) & (centers < fit.fit_range[1])
    total = values[mask].sum() / max(frac[mask].sum(), 1e-300)
    predicted = frac * total
    # the curve: the PDF in sub-bins, in events per channel bin
    nsub = 10
    fine = np.concatenate([np.linspace(lo, hi, nsub, endpoint=False) for lo, hi in zip(edges[:-1], edges[1:])]
                          + [edges[-1:]])
    fine_frac = bin_fractions(fit.pdf, ch.obs, fine)
    width = np.repeat(np.diff(edges), nsub)
    curve_y = fine_frac * total * width / np.diff(fine)
    curve_x = 0.5 * (fine[1:] + fine[:-1])
    c, top, bottom = _split_canvas()
    keep = []
    top.cd()
    h.SetMarkerStyle(20)
    h.SetLineWidth(2)
    h.SetLineColor(ROOT.kBlack)
    h.SetTitle(f";;{_bin_width_title(ch.obs.GetTitle(), edges)}")
    _axes(h, 0.06, 0.065, 0.95)
    h.GetXaxis().SetLabelSize(0)
    ymax = max((values + errors).max(), curve_y.max(), 1e-9)
    if logy:
        lo, hi = _log_range(values, ymax)
        h.SetMinimum(lo)
        h.SetMaximum(hi)
    else:
        h.SetMinimum(0.)
        h.SetMaximum(1.2 * ymax)
    curve = ROOT.TGraph(len(curve_x), np.ascontiguousarray(curve_x), np.ascontiguousarray(curve_y))
    color = r.channel.processes[process].color
    curve.SetLineColor(ROOT.kBlue if color is None else _color(r, process, 0))
    curve.SetLineWidth(3)
    h.Draw("E")
    curve.Draw("L same")
    h.Draw("E same")
    keep.append(curve)
    _header(top, label, keep, SPLIT_SCALE)
    top.SetLogy(logy)
    top.RedrawAxis()
    # pulls in the bins with MC entries (only the fit range when the fit used one)
    bottom.cd()
    shown = mask & (values > 0) & (errors > 0)
    pull = np.where(shown, (values - predicted) / np.where(shown, errors, 1.), 0.)
    span = max(3.5, 1.15 * np.abs(pull[shown]).max()) if shown.any() else 3.5
    _lower_frame(edges, ch.obs.GetTitle(), "#frac{N_{MC} - N_{fit}}{#sigma_{MC}}", -span, span, keep)
    _lower_points("_fit_pull", edges, pull, np.where(shown, 1., 0.), keep)
    c.SaveAs(path)
    c.Close()


def plot_template(r, process, path, logy=False, label=None):
    """A binned template process at its expected yield, with the MC statistical uncertainties of its bins. label: the
    header (see plot_channel)."""
    import ROOT
    from .hist import arrays
    ch = r.channel
    edges = ch.edges
    _, values, variances = arrays(r.templates[process])
    errors = np.sqrt(variances)
    norm = r.expected[process].sum() / max(values.sum(), 1e-300)
    values, errors = values * norm, errors * norm
    c = _single_canvas()
    keep = []
    color = _color(r, process, 0, signal_index=0 if ch.processes[process].signal else None)
    line = _hist("_template_line", edges, values)
    line.SetLineColor(color)
    line.SetLineWidth(3)
    points = _hist("_template_points", edges, values, errors)
    points.SetMarkerStyle(20)
    points.SetLineWidth(2)
    points.SetLineColor(ROOT.kBlack)
    ymax = max((values + errors).max(), 1e-9)
    line.SetTitle(f";{ch.obs.GetTitle()};{_bin_width_title(ch.obs.GetTitle(), edges)}")
    _axes(line, 0.05, 0.05, 1.2)
    if logy:
        lo, hi = _log_range(values, ymax)
        line.SetMinimum(lo)
        line.SetMaximum(hi)
    else:
        line.SetMinimum(0.)
        line.SetMaximum(1.3 * ymax)
    line.Draw("hist")
    points.Draw("E same")
    leg = ROOT.TLegend(0.6, 0.78, 0.93, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(FONT)
    leg.SetTextSize(0.043)
    leg.AddEntry(line, f"{ch.processes[process].title} ({r.expected[process].sum():.3g})", "L")
    rel = r.mc_rel.get(process)
    leg.AddEntry(points, "MC" + (f" ({100. * rel:.2g}% stat.)" if rel is not None else ""), "PE")
    leg.Draw()
    _header(c, label, keep)
    c.SetLogy(logy)
    c.RedrawAxis()
    c.SaveAs(path)
    c.Close()


def plot_inputs(r, path, logy=False, label=None):
    """Each process's expected distribution as a line over the data, with the pulls of the data from the total
    expectation, (data - expected) / sqrt(data). label: the header (see plot_channel)."""
    import ROOT
    from .hist import arrays
    ch = r.channel
    edges = ch.edges
    c, top, bottom = _split_canvas()
    keep = []
    top.cd()
    sig = [p for p in r.processes if ch.processes[p].signal]
    lines = []
    for i, p in enumerate(r.processes):
        h = _hist(f"_input_{p}", edges, r.expected[p])
        h.SetLineColor(_color(r, p, i, signal_index=sig.index(p) if p in sig else None))
        h.SetLineWidth(2)
        lines.append((p, h))
    data_vals = arrays(r.data)[1]
    data_err = np.sqrt(np.clip(data_vals, 0., None))
    data = _hist("_inputs_data", edges, data_vals, data_err)
    data.SetMarkerStyle(20)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    expected = [r.expected[p] for p in r.processes]
    ymax = max([data_vals.max()] + [e.max() for e in expected] + [1e-9])
    positive = np.concatenate([e[e > 0] for e in expected] + [np.zeros(0)])
    ymin = max(0.3 * positive.min(), 1e-6 * ymax) if positive.size else 1e-6 * ymax
    frame = _hist("_inputs_frame", edges, np.zeros(len(edges) - 1))
    frame.SetTitle(f";;{_bin_width_title(ch.obs.GetTitle(), edges)}")
    _axes(frame, 0.06, 0.065, 0.95)
    frame.GetXaxis().SetLabelSize(0)
    frame.SetMaximum(ymin * (ymax / ymin) ** (1. / 0.7) if logy else 1.5 * ymax)
    frame.SetMinimum(ymin if logy else 0.)
    frame.Draw("axis")
    for _, h in lines:
        h.Draw("hist same")
    data.Draw("E same")
    leg = ROOT.TLegend(0.03 + top.GetLeftMargin(), 0.72, 0.98 - top.GetRightMargin(), 0.97 - top.GetTopMargin())
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(FONT)
    leg.SetTextSize(0.05)
    leg.AddEntry(data, f"{'Asimov d' if r.data_source == 'asimov' else 'D'}ata", "PE")
    for p, h in lines:
        leg.AddEntry(h, ch.processes[p].title, "L")
    leg.Draw()
    _header(top, label, keep, SPLIT_SCALE)
    top.SetLogy(logy)
    top.RedrawAxis()
    bottom.cd()
    total = np.sum(expected, axis=0)
    shown = data_err > 0
    pull = np.where(shown, (data_vals - total) / np.where(shown, data_err, 1.), 0.)
    span = max(3., 1.2 * np.abs(pull).max()) if shown.any() else 3.
    _lower_frame(edges, ch.obs.GetTitle(), "#frac{N_{data} - N_{exp}}{#sigma_{data}}", -span, span, keep)
    _lower_points("_inputs_pull", edges, pull, np.where(shown, 1., 0.), keep)
    c.SaveAs(path)
    c.Close()


FIT_DIRS = {"prefit": "Pre-fit", "fit_b": "Background-only fit", "fit_s": "Signal+background fit"}


def _graph_values(g, widths):
    n = g.GetN()
    y = np.array([g.GetY()[i] for i in range(n)]) * widths[:n]
    lo = np.array([g.GetErrorYlow(i) for i in range(n)]) * widths[:n]
    hi = np.array([g.GetErrorYhigh(i) for i in range(n)]) * widths[:n]
    return y, lo, hi


def read_fitdiagnostics(path, model_dir="."):
    """Shapes from a FitDiagnostics --saveShapes output, in events per bin on the channels' own binning.

    Returns {fit: {channel: {"edges", "processes": {p: (values, errors)}, "total", "total_background",
    "data": (values, err_lo, err_hi)}}} for the fits present ("prefit", "fit_b", "fit_s"). Combine stores workspace
    channels per unit of the observable and binned channels per bin index; both are converted with the stored axis.
    """
    import ROOT
    with open(os.path.join(model_dir, "model.json")) as f:
        model = json.load(f)
    tf = ROOT.TFile.Open(path)
    if not tf or tf.IsZombie():
        raise IOError(f"Unable to open {path}")
    out = {}
    for fit in FIT_DIRS:
        top = tf.Get(f"shapes_{fit}")
        if not top:
            continue
        out[fit] = {}
        for name, info in model["channels"].items():
            d = top.Get(name)
            if not d:
                continue
            edges = np.asarray(info["observable"]["edges"])
            ch = {"edges": edges, "processes": {}}
            widths = None
            for key in d.GetListOfKeys():
                obj = key.ReadObj()
                if not obj.InheritsFrom("TH1") or obj.InheritsFrom("TH2"):
                    continue
                n = obj.GetNbinsX()
                if n != len(edges) - 1:
                    raise ValueError(f"{name}/{key.GetName()}: {n} bins in {path}, {len(edges) - 1} in model.json")
                widths = np.array([obj.GetXaxis().GetBinWidth(b) for b in range(1, n + 1)])
                vals = np.array([obj.GetBinContent(b) for b in range(1, n + 1)]) * widths
                errs = np.array([obj.GetBinError(b) for b in range(1, n + 1)]) * widths
                if key.GetName() in ("total", "total_signal", "total_background"):
                    ch[key.GetName()] = (vals, errs)
                else:
                    ch["processes"][key.GetName()] = (vals, errs)
            g = d.Get("data")
            if g and widths is not None:
                y, lo, hi = _graph_values(g, widths)
                ch["asimov"] = bool(np.any(np.abs(y - np.round(y)) > 1e-6))
                if ch["asimov"]:      # Poisson intervals mean nothing for non-integer (Asimov) counts: use sqrt(N)
                    lo = hi = np.sqrt(np.clip(y, 0., None))
                ch["data"] = (y, lo, hi)
            out[fit][name] = ch
    tf.Close()
    return out, model


def _hist(name, edges, values, errors=None):
    from .hist import make
    return make(edges, values, None if errors is None else np.square(errors), name)


def plot_postfit_channel(ch, info, fit, path, logy=False):
    """Stacked backgrounds with the total uncertainty band, signals as lines, data, and a data / prediction ratio."""
    import ROOT
    edges = ch["edges"]
    procs = info["processes"]
    shown = [p for p in ch["processes"] if p in procs and ch["processes"][p][0].sum() > 0]
    bkg = [p for p in shown if not procs[p]["signal"]]
    sig = [p for p in shown if procs[p]["signal"]]
    sig_colors = [ROOT.kRed, ROOT.kMagenta + 1, ROOT.kGreen + 2, ROOT.kCyan + 1]
    total_key = "total" if fit == "fit_s" else "total_background"
    total, total_err = ch.get(total_key, ch.get("total"))
    c = ROOT.TCanvas("c", "", 800, 750)
    top = ROOT.TPad("top", "", 0, 0.3, 1, 1)
    bottom = ROOT.TPad("bottom", "", 0, 0, 1, 0.3)
    top.SetBottomMargin(0.02)
    bottom.SetTopMargin(0.03)
    bottom.SetBottomMargin(0.32)
    top.Draw()
    bottom.Draw()
    top.cd()
    keep = []
    order = sorted(bkg, key=lambda p: ch["processes"][p][0].sum())
    stacked = (order + sig) if fit == "fit_s" else order
    cumulative = np.zeros(len(edges) - 1)
    stack = []
    for i, p in enumerate(stacked):
        cumulative = cumulative + ch["processes"][p][0]
        h = _hist(f"_pf_{p}", edges, cumulative)
        color = sig_colors[sig.index(p) % 4] - 7 if p in sig else ROOT.TColor.GetColor(COLORS[i % len(COLORS)])
        h.SetFillColor(color)
        h.SetLineColor(color)
        stack.append((p, h))
    band = _hist("_pf_band", edges, total, total_err)
    band.SetFillColor(ROOT.kGray + 2)
    band.SetFillStyle(3354)
    band.SetMarkerSize(0)
    band.SetLineColor(ROOT.kGray + 2)
    lines = []
    if fit != "fit_s":
        for i, p in enumerate(sig):
            h = _hist(f"_pf_sig_{p}", edges, ch["processes"][p][0])
            h.SetLineColor(sig_colors[i % 4])
            h.SetLineWidth(3)
            lines.append((p, h))
    data = None
    ymax = max([total.max()] + [h.GetMaximum() for _, h in lines] + [1e-9])
    if "data" in ch:
        y, lo, hi = ch["data"]
        centers = 0.5 * (edges[1:] + edges[:-1])
        n = len(y)
        data = ROOT.TGraphAsymmErrors(n, np.ascontiguousarray(centers[:n]), np.ascontiguousarray(y),
                                      np.zeros(n), np.zeros(n), np.ascontiguousarray(lo), np.ascontiguousarray(hi))
        data.SetMarkerStyle(20)
        data.SetLineWidth(2)
        ymax = max(ymax, (y + hi).max())
    frame = _hist("_pf_frame", edges, np.zeros(len(edges) - 1))
    frame.GetYaxis().SetTitle("Events / bin")
    frame.GetXaxis().SetLabelSize(0)
    frame.SetMaximum(ymax * (30. if logy else 1.4))
    frame.SetMinimum(max(1e-4 * ymax, 1e-6) if logy else 0.)
    frame.Draw("axis")
    for _, h in reversed(stack):
        h.Draw("hist same")
    band.Draw("E2 same")
    for _, h in lines:
        h.Draw("hist same")
    if data:
        data.Draw("P same")
    if info.get("blind"):
        box = ROOT.TBox(info["blind"][0], frame.GetMinimum(), info["blind"][1], frame.GetMaximum())
        box.SetFillColorAlpha(ROOT.kGray, 0.35)
        box.Draw()
        keep.append(box)
    leg = ROOT.TLegend(0.6, 0.5, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetHeader(FIT_DIRS[fit])
    if data:
        leg.AddEntry(data, f"{'Asimov ' if ch['asimov'] else ''}data ({ch['data'][0].sum():.3g})",
                     "PE")
    for p, h in reversed(stack):
        leg.AddEntry(h, f"{procs[p]['title']} ({ch['processes'][p][0].sum():.3g})", "F")
    for p, h in lines:
        leg.AddEntry(h, f"{procs[p]['title']} ({ch['processes'][p][0].sum():.3g})", "L")
    leg.AddEntry(band, "Uncertainty", "F")
    leg.Draw()
    top.SetLogy(logy)
    top.RedrawAxis()
    bottom.cd()
    safe = np.where(total > 0, total, 1.)
    rband = _hist("_pf_rband", edges, np.where(total > 0, 1., 0.), np.where(total > 0, total_err / safe, 0.))
    rband.SetFillColor(ROOT.kGray + 2)
    rband.SetFillStyle(3354)
    rband.SetMarkerSize(0)
    rframe = _hist("_pf_rframe", edges, np.ones(len(edges) - 1))
    rframe.GetXaxis().SetTitle(info["observable"]["title"])
    rframe.GetYaxis().SetTitle("Data / pred.")
    for ax, size in ((rframe.GetXaxis(), 0.12), (rframe.GetYaxis(), 0.1)):
        ax.SetTitleSize(size)
        ax.SetLabelSize(size * 0.9)
    rframe.GetYaxis().SetTitleOffset(0.45)
    rframe.GetYaxis().SetNdivisions(505)
    rframe.SetMinimum(0.)
    rframe.SetMaximum(2.)
    rframe.SetLineColor(ROOT.kBlack)
    rframe.SetLineStyle(2)
    rframe.Draw("hist")
    rband.Draw("E2 same")
    if data:
        y, lo, hi = ch["data"]
        n = len(y)
        ratio = ROOT.TGraphAsymmErrors(n, np.ascontiguousarray(0.5 * (edges[1:] + edges[:-1])[:n]),
                                       np.ascontiguousarray(y / safe[:n]), np.zeros(n), np.zeros(n),
                                       np.ascontiguousarray(lo / safe[:n]), np.ascontiguousarray(hi / safe[:n]))
        ratio.SetMarkerStyle(20)
        ratio.SetLineWidth(2)
        ratio.Draw("P same")
        keep.append(ratio)
    c.SaveAs(path)
    c.Close()


def plot_fitdiagnostics(path, model_dir=".", figdir=None, logy=False):
    """Pre- and post-fit figures (<figdir>/<channel>_<fit>.png, default <model_dir>/figures/postfit) from a
    FitDiagnostics output made with --saveShapes (add --saveWithUncertainties for the bands). Also writes
    postfit_yields.json there. Returns {fit: {channel: {process: yield}}}."""
    shapes, model = read_fitdiagnostics(path, model_dir)
    if not shapes:
        raise IOError(f"No shapes_* directories in {path}: run FitDiagnostics with --saveShapes")
    figdir = figdir or os.path.join(model_dir, "figures", "postfit")
    os.makedirs(figdir, exist_ok=True)
    _style()
    yields = {}
    for fit, channels in shapes.items():
        yields[fit] = {}
        for name, ch in channels.items():
            plot_postfit_channel(ch, model["channels"][name], fit, os.path.join(figdir, f"{name}_{fit}.png"), logy)
            y = {p: float(v.sum()) for p, (v, _) in ch["processes"].items()}
            if "data" in ch:
                y["data"] = float(ch["data"][0].sum())
            yields[fit][name] = y
    with open(os.path.join(figdir, "postfit_yields.json"), "w") as f:
        json.dump(yields, f, indent=2)
    return yields


def main(argv=None):
    parser = argparse.ArgumentParser(description="Pre/post-fit figures from a Combine FitDiagnostics --saveShapes output")
    parser.add_argument("fitdiagnostics", help="fitDiagnostics.<name>.root")
    parser.add_argument("--model-dir", default=".", help="the mu2ecombine output directory (with model.json)")
    parser.add_argument("--figdir", default=None, help="default: <model-dir>/figures/postfit")
    parser.add_argument("--logy", action="store_true")
    args = parser.parse_args(argv)
    yields = plot_fitdiagnostics(args.fitdiagnostics, args.model_dir, args.figdir, args.logy)
    for fit, channels in yields.items():
        for name, y in channels.items():
            print(f"{fit:7s} {name:10s} " + "  ".join(f"{p}={v:.4g}" for p, v in y.items()))


if __name__ == "__main__":
    main()
