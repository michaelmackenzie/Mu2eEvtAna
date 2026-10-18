"""Model figures: the expected (pre-fit) stack of each channel and the MC shape fits, and post-fit figures from Combine.

    model.write(outdir, plots=True)                                  # figures in <outdir>/figures
    python3 -m mu2ecombine.plots fitDiagnostics.<name>.root          # in <outdir>, after FitDiagnostics --saveShapes

Stacks are drawn as cumulative histograms rather than a THStack.
"""

import argparse
import json
import os

import numpy as np

COLORS = [ "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581",
           "#92dadd"]


def _style():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning


def plot_model(results, figdir):
    os.makedirs(figdir, exist_ok=True)
    _style()
    for name, r in results.items():
        plot_channel(r, os.path.join(figdir, f"{name}_prefit.png"))
        for p, fit in r.fits.items():
            if p in r.processes:
                plot_fit(r, p, fit, os.path.join(figdir, f"{name}_{p}_fit.png"))


def plot_channel(r, path, logy=False):
    import ROOT
    from .hist import arrays
    ch = r.channel
    edges = ch.edges
    bkg = [p for p in r.processes if not ch.processes[p].signal]
    sig = [p for p in r.processes if ch.processes[p].signal]
    c = ROOT.TCanvas("c", "", 800, 600)
    keep = []
    order = sorted(bkg, key=lambda p: r.expected[p].sum())
    cumulative = np.zeros(len(edges) - 1)
    stack = []
    for i, p in enumerate(order):
        cumulative = cumulative + r.expected[p]
        h = _hist(f"_stack_{p}", edges, cumulative)
        color = ROOT.TColor.GetColor(COLORS[i % len(COLORS)])
        h.SetFillColor(color)
        h.SetLineColor(color)
        stack.append((p, h))
    data_vals = arrays(r.data)[1]
    data = _hist("_data", edges, data_vals, np.sqrt(np.clip(data_vals, 0., None)))
    data.SetMarkerStyle(20)
    data.SetLineColor(ROOT.kBlack)
    sig_hists = []
    for i, p in enumerate(sig):
        h = _hist(f"_sig_{p}", edges, r.expected[p])
        h.SetLineColor([ROOT.kRed, ROOT.kMagenta + 1, ROOT.kGreen + 2][i % 3])
        h.SetLineWidth(3)
        sig_hists.append((p, h))
    ymax = max([cumulative.max(), data_vals.max()] + [h.GetMaximum() for _, h in sig_hists] + [1e-9])
    frame = _hist("_frame", edges, np.zeros(len(edges) - 1))
    frame.GetXaxis().SetTitle(ch.obs.GetTitle())
    frame.GetYaxis().SetTitle("Events / bin")
    frame.SetMaximum(ymax * (30. if logy else 1.35))
    frame.SetMinimum(max(1e-4 * ymax, 1e-6) if logy else 0.)
    frame.Draw("axis")
    for _, h in reversed(stack):
        h.Draw("hist same")
    for _, h in sig_hists:
        h.Draw("hist same")
    data.Draw("E1 same")
    if ch.blind:
        box = ROOT.TBox(ch.blind[0], frame.GetMinimum(), ch.blind[1], frame.GetMaximum())
        box.SetFillColorAlpha(ROOT.kGray, 0.35)
        box.Draw()
        keep.append(box)
    leg = ROOT.TLegend(0.62, 0.55, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.AddEntry(data, f"{'Asimov ' if r.data_source == 'asimov' else ''}data ({data_vals.sum():.3g})", "PE")
    for p, h in reversed(stack):
        leg.AddEntry(h, f"{ch.processes[p].title} ({r.expected[p].sum():.3g})", "F")
    for p, h in sig_hists:
        leg.AddEntry(h, f"{ch.processes[p].title} ({r.expected[p].sum():.3g})", "L")
    leg.Draw()
    c.SetLogy(logy)
    c.RedrawAxis()
    c.SaveAs(path)
    c.Close()


def plot_fit(r, process, fit, path):
    """The fitted function over the histogram it was fit to (both in the channel binning)."""
    import ROOT
    from .hist import arrays
    from .pdfs import bin_fractions
    edges = r.channel.edges
    h = fit.fit_hist(edges)
    _, values, _ = arrays(h)
    centers = 0.5 * (edges[1:] + edges[:-1])
    frac = bin_fractions(fit.pdf, r.channel.obs, edges)
    if fit.fit_range is not None:
        mask = (centers > fit.fit_range[0]) & (centers < fit.fit_range[1])
        total = values[mask].sum() / max(frac[mask].sum(), 1e-300)
    else:
        total = values.sum()
    c = ROOT.TCanvas("c", "", 800, 600)
    h.SetMarkerStyle(20)
    h.SetLineColor(ROOT.kBlack)
    h.GetXaxis().SetTitle(r.channel.obs.GetTitle())
    h.GetYaxis().SetTitle("Events / bin")
    h.SetMaximum(1.35 * max(values.max(), (frac * total).max()))
    h.SetMinimum(0.)
    f = _hist("_fit_f", edges, frac * total)
    f.SetLineColor(ROOT.kBlue)
    f.SetLineWidth(2)
    h.Draw("E1")
    f.Draw("hist same")
    text = ROOT.TPaveText(0.55, 0.65, 0.89, 0.89, "NDC")
    text.SetFillStyle(0)
    text.SetBorderSize(0)
    text.SetTextAlign(12)
    text.AddText(f"{process}: {fit.info.get('function')} (status {fit.info.get('fit_status')})")
    for k, v in list(fit.info.get("params", {}).items())[:7]:
        text.AddText(f"{k.split(r.channel.name + '_' + process + '_')[-1]} = {v:.4g}")
    text.Draw()
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
