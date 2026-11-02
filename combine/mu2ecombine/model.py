"""Model, channels, processes, and nuisances -> Combine datacards and shape files.

    p = observable("p", 100, 110, bin_width=0.5, title="p [MeV/c]")   # a RooRealVar with its binning
    model = Model("bnv", poi_range=(0, 50))
    sr = model.channel("sr", p, blind=(102, 105))
    sr.add("signal", Fit("dscb", sig_th1), signal=True)
    sr.add("dio", dio_th1)                                       # a histogram is a binned template
    sr.add("cosmic", Template(cosmic_sr, shape=[cosmic_sr, cosmic_cr]))
    model.lnN("lumi", 1.1, processes="!cosmic")
    model.shape_shift("pscale", shift=0.1)                       # +-0.1 shift of every MC shape
    model.write("output/bnv")

A channel whose processes are all binned templates is written as TH1 shapes; a channel with any PDF is written as a
RooWorkspace, with its templates as RooDataHist on the channel observable. workspace=True (Model or channel) writes a
template channel as a workspace too, so that its data sets, e.g. toys generated from it, live on the channel observable
and can be fit by a PDF model of the same channel.

A channel on two observables, model.channel("sr", (p, t)), is a unrolled 1D channel: its templates are TH1Ds on the bin
index k = iy * nx + ix (x fastest, nx the bins of the first observable), made from TH2 inputs or a Product of two 1D
shapes (shapes.py). Only binned templates are supported in 2D channels (so far).

Process selections (processes=, channels=) are a name, a glob, or a list of them; '!pattern' excludes.
"""

import fnmatch
import json
import os
import shutil
import subprocess
import warnings

import numpy as np

from .datacard import Datacard, kappa_text
from .envelope import Envelope
from .hist import (HistError, arrays, as_th1, as_th2, clip_negative, integral, is_th1, is_th2, make, rebin2_to, rebin_to,
                   rel_error, transformed, transformed2, unroll)
from .shapes import Fit, Pdf, Product, ShapeError, Template, as_shape, in_channel, scaled


class ModelError(ValueError):
    pass


def matches(name, patterns):
    """True if name matches the selection (None or '*' select everything; '!pattern' excludes)."""
    if patterns is None:
        return True
    if isinstance(patterns, str):
        patterns = [patterns]
    positive = [p for p in patterns if not p.startswith("!")]
    negative = [p[1:] for p in patterns if p.startswith("!")]
    if positive and not any(fnmatch.fnmatchcase(name, p) for p in positive):
        return False
    return not any(fnmatch.fnmatchcase(name, p) for p in negative)


def observable(name, lo=None, hi=None, bins=None, bin_width=None, edges=None, title=None):
    """A RooRealVar with a binning: (lo, hi) with bins or bin_width, or explicit edges. Equivalent to building the
    RooRealVar yourself and calling setBins / setBinning."""
    import ROOT
    if edges is not None:
        edges = np.asarray(edges, dtype=float)
    else:
        if lo is None or hi is None or (bins is None) == (bin_width is None):
            raise ModelError(f"observable {name}: give lo, hi and one of bins / bin_width, or edges")
        n = int(bins) if bins is not None else int(round((hi - lo) / bin_width))
        if bin_width is not None and abs(n * bin_width - (hi - lo)) > 1e-6 * bin_width:
            raise ModelError(f"observable {name}: the range {lo}-{hi} is not a multiple of {bin_width}")
        edges = np.linspace(lo, hi, n + 1)
    x = ROOT.RooRealVar(name, title or name, float(edges[0]), float(edges[-1]))
    widths = np.diff(edges)
    if np.allclose(widths, widths[0]):
        x.setBins(len(widths))
    else:
        x.setBinning(ROOT.RooBinning(len(widths), np.ascontiguousarray(edges)))
    return x


def binning(obs):
    """Bin edges of a RooRealVar's default binning."""
    b = obs.getBinning()
    n = b.numBins()
    return np.array([b.binLow(i) for i in range(n)] + [b.binHigh(n - 1)])


class Process:
    def __init__(self, name, shape, signal=False, poi=None, title=None, color=None):
        self.name, self.shape, self.signal, self.poi = name, shape, signal, poi
        self.title = title or name
        self.color = color   # figure color: a ROOT color index or "#rrggbb" (default: the figure palette)


class Channel:
    """A channel on one RooRealVar, or on two (x, y): then `axes` are the two variables with their binnings in
    `axis_edges`, `edges` is the unrolled bin-index binning 0..nx*ny (global bin k = iy * nx + ix), and `obs` is the
    bin-index variable <name>_bin used in workspaces. For one observable, axes = (obs,) and edges = its binning."""

    def __init__(self, model, name, obs, data=None, blind=None, title=None, workspace=None):
        axes = tuple(obs) if isinstance(obs, (tuple, list)) else (obs,)
        for x in axes:
            try:
                ok = bool(x.InheritsFrom("RooRealVar"))
            except AttributeError:
                ok = False
            if not ok:
                raise ModelError(f"channel {name}: the observable must be a RooRealVar (see mu2ecombine.observable)")
        if len(axes) not in (1, 2):
            raise ModelError(f"channel {name}: give one RooRealVar, or a tuple of two for a 2D channel")
        if len(axes) == 2 and axes[0].GetName() == axes[1].GetName():
            raise ModelError(f"channel {name}: the two observables have the same name '{axes[0].GetName()}'")
        self.model, self.name, self.axes = model, name, axes
        self.axis_edges = tuple(binning(x) for x in axes)
        self.dim = len(axes)
        if self.dim == 2:
            import ROOT
            n = (len(self.axis_edges[0]) - 1) * (len(self.axis_edges[1]) - 1)
            self.obs = ROOT.RooRealVar(f"{name}_bin", f"{name} bin", 0., float(n))
            self.obs.setBins(n)
            self.edges = np.arange(n + 1, dtype=float)
        else:
            self.obs, self.edges = axes[0], self.axis_edges[0]
        self.data = self._hist(data, f"channel {name} data") if data is not None else None
        self.blind = tuple(blind) if blind is not None else None
        self.title = title or name
        self.workspace = workspace
        self.processes = {}

    def _hist(self, h, what):
        """A histogram on the channel binning as a detached TH1D: a TH1 (list) for 1D; a TH2 (list) rebinned to the
        axes and unrolled for 2D."""
        try:
            if self.dim == 1:
                if is_th2(h) or (isinstance(h, (list, tuple)) and h and is_th2(h[0])):
                    raise ModelError(f"{what}: a TH2 was given to the 1D channel {self.name}")
                return as_th1(h, what)
            if is_th1(h) or (isinstance(h, (list, tuple)) and h and is_th1(h[0])):
                raise ModelError(f"{what}: a TH1 was given to the 2D channel {self.name}; give a TH2")
            return unroll(rebin2_to(as_th2(h, what), *self.axis_edges))
        except HistError as err:
            raise ModelError(f"{what}: {err}") from None

    def add(self, name, shape, signal=False, poi=None, title=None, color=None):
        """Add a process. shape: a TH1 or list of TH1s (a binned template), Template, Fit, Pdf, or Envelope; in a 2D
        channel a TH2 (list), a Template of TH2s, or a Product of two 1D shapes. color: its figure color, a ROOT
        color index or "#rrggbb" (default: the figure palette)."""
        if name in self.processes:
            raise ModelError(f"channel {self.name}: process '{name}' already defined")
        if not name.replace("_", "").isalnum():
            raise ModelError(f"process name '{name}': use letters, digits, and underscores")
        if not isinstance(shape, (Template, Fit, Pdf, Envelope, Product)):
            try:
                shape = as_shape(shape, f"{self.name}/{name}")
            except (ShapeError, HistError) as err:
                raise ModelError(str(err)) from None
        if isinstance(shape, Envelope) and signal:
            raise ModelError(f"channel {self.name}: an envelope cannot be the signal")
        where = f"channel {self.name}/{name}"
        if self.dim == 2:
            if isinstance(shape, (Fit, Pdf, Envelope)):
                raise ModelError(f"{where}: a {type(shape).__name__} cannot be a process of a 2D channel by itself; "
                                 "use Product(x_shape, y_shape) with it as a component")
            if isinstance(shape, Template):
                if shape.dim != 2:
                    raise ModelError(f"{where}: a TH1 template was given to the 2D channel; give TH2s or a Product")
                try:
                    shape = shape.unrolled(*self.axis_edges)
                except (ShapeError, HistError) as err:
                    raise ModelError(f"{where}: {err}") from None
        else:
            if isinstance(shape, Product):
                raise ModelError(f"{where}: a Product needs a 2D channel (give model.channel a tuple of two observables)")
            if isinstance(shape, Template) and shape.dim == 2:
                raise ModelError(f"{where}: a TH2 template was given to the 1D channel")
        p = Process(name, shape, signal, poi, title, color)
        self.processes[name] = p
        return p

    def set_data(self, hist):
        self.data = self._hist(hist, f"channel {self.name} data")


class _Nuisance:
    def __init__(self, name, kind, processes="*", channels="*", group=None, **options):
        self.name, self.kind, self.processes, self.channels, self.group = name, kind, processes, channels, group
        self.options = options

    def applies(self, channel, process):
        return matches(channel, self.channels) and matches(process, self.processes)


class Model:
    def __init__(self, name="model", poi="r", poi_range=(0., 20.), workspace=False):
        self.name = name
        self.workspace = workspace   # write template-only channels as workspaces too
        self.default_poi = poi
        self.pois = {poi: tuple(poi_range)}
        self.channels = {}
        self.nuisances = []          # lnN and generated shape shifts, in order
        self.rate_params = []        # (name, processes, channels, init, range)
        self.params = {}             # name -> (mean, sigma)
        self.groups = {}             # name -> [nuisances]
        self.mcstats = {"enabled": True, "auto": False, "threshold": 10., "min_rel": 1e-4}
        self.free = set()            # unconstrained (flat-prior) nuisances
        self.figure_label = None     # (left, right) TLatex header of the stack figures, e.g. an experiment label and
                                     # the exposure; either may be None

    # -- definition ---------------------------------------------------------------------------------------------------
    def channel(self, name, obs, data=None, blind=None, title=None, workspace=None):
        """Add a channel on the RooRealVar obs (its binning is the channel binning), or on a tuple (x, y) of two for a
        2D channel (see Channel). data: a TH1, or a TH2 for a 2D channel (default: the Asimov background
        expectation); blind: (lo, hi) window excluded from envelope fits; workspace: write the channel as a
        RooWorkspace even if all its processes are templates (default: the Model's setting)."""
        if name in self.channels:
            raise ModelError(f"channel '{name}' already defined")
        ch = Channel(self, name, obs, data, blind, title, workspace)
        self.channels[name] = ch
        return ch

    def add_poi(self, name, range=(0., 20.)):
        """Another parameter of interest; signals select it with add(..., signal=True, poi=name)."""
        self.pois[name] = tuple(range)

    def lnN(self, name, value, processes="*", channels="*", group=None):
        """Log-normal rate uncertainty: value is kappa (e.g. 1.1) or (kappa_up, kappa_down)."""
        self.nuisances.append(_Nuisance(name, "lnN", processes, channels, group, value=value))

    def shape_shift(self, name, shift=0., stretch=0., about=None, processes="*", channels="*", group=None, axis=None):
        """A generated shape uncertainty: up/down are the MC shapes transformed by x -> about + (x - about)(1 +- stretch)
        +- shift (energy scale / resolution). Applies to Template and Fit processes (refit for fitted PDFs).
        axis: the name of the observable the shift acts on. Required in 2D channels, where it shifts the TH2 Templates
        (every slice) and the Product component along that axis; a 1D channel whose observable has another name is
        skipped, so one shift with axis="p" covers 1D p channels and 2D (p, t) channels."""
        self.nuisances.append(_Nuisance(name, "shift", processes, channels, group, shift=shift, stretch=stretch,
                                        about=about, axis=axis))

    def rate_param(self, name, processes, channels="*", init=1., range=(0., 10.)):
        """A free normalization factor shared by the selected processes (e.g. constrained by a control region)."""
        self.rate_params.append((name, processes, channels, float(init), tuple(range)))

    def param(self, name, mean=0., sigma=1.):
        """A Gaussian constraint on a parameter of your own PDFs (Combine 'param'). A Fit keeps such parameters at
        their values instead of fitting them."""
        self.params[name] = (float(mean), float(sigma))

    def unconstrained(self, *names):
        """Give nuisances a flat prior (no constraint term), e.g. a shape uncertainty whose size the data should
        determine: shape nuisances are written as 'shapeU', nuisances only inside workspaces as 'flatParam'. The
        up/down variations then only set the scale of the parameter (theta = +-1). Log-normal nuisances cannot be
        unconstrained."""
        self.free.update(names)

    def group(self, name, nuisances):
        self.groups.setdefault(name, []).extend([nuisances] if isinstance(nuisances, str) else nuisances)

    def mc_stats(self, enabled=True, auto=False, threshold=10., min_rel=1e-4):
        """MC statistical uncertainties: one log-normal per process (relative MC uncertainty of its yield, if above
        min_rel). auto=True uses Combine's autoMCStats (with threshold) in binned channels instead; it is off by default
        because it stalls fits in some Combine builds."""
        self.mcstats = {"enabled": enabled, "auto": auto, "threshold": float(threshold), "min_rel": float(min_rel)}

    def nuisance_names(self):
        """Every nuisance name defined in the model (lnN, shape shifts, params, and process variations/nuisances)."""
        names = {n.name for n in self.nuisances} | set(self.params)
        for ch in self.channels.values():
            for p in ch.processes.values():
                if isinstance(p.shape, Product):
                    names |= p.shape.nuisance_names()
                    continue
                names |= set(getattr(p.shape, "variations", {}))
                names |= set(getattr(p.shape, "nuisances", {}))
        return names

    # -- building ------------------------------------------------------------------------------------------------------
    def build(self):
        """Build every channel (fits, templates, variations). Returns {channel: ChannelResult}."""
        self._check_observables()
        return {name: _ChannelBuilder(self, ch).build() for name, ch in self.channels.items()}

    def _check_observables(self):
        seen = {}
        for ch in self.channels.values():
            for obs, edges in zip(ch.axes, ch.axis_edges):
                name = obs.GetName()
                if name in seen and (len(seen[name][1]) != len(edges) or not np.allclose(seen[name][1], edges)):
                    raise ModelError(f"channels {seen[name][0]} and {ch.name} use observables named '{name}' with "
                                     "different ranges or binnings: give each its own RooRealVar name")
                seen.setdefault(name, (ch.name, edges))

    def write(self, outdir, plots=False, text2workspace=False, shapes=True):
        """Write the Combine model to outdir: combined.txt and <channel>.txt cards (paths relative to outdir), shapes/,
        model.json, and commands.sh. shapes=False writes counting cards, without shapes/ or shapes lines: every channel
        must then have one bin and only log-normal uncertainties. Returns the summary dict."""
        results = self.build()
        if shapes:
            os.makedirs(os.path.join(outdir, "shapes"), exist_ok=True)
            for r in results.values():
                r.write_shapes(outdir)
        else:
            self._check_counting(results)
            os.makedirs(outdir, exist_ok=True)
        ids = self._process_ids(results)
        cards = {name: self._card([r], ids, shapes) for name, r in results.items()}
        cards["combined"] = self._card(list(results.values()), ids, shapes)
        for name, card in cards.items():
            card.write(os.path.join(outdir, f"{name}.txt"))
        summary = self._summary(results, ids)
        with open(os.path.join(outdir, "model.json"), "w") as f:
            json.dump(summary, f, indent=2)
        with open(os.path.join(outdir, "commands.sh"), "w") as f:
            f.write(self.commands(results))
        self.print_summary(summary)
        if plots:
            from . import plots as plotlib
            plotlib.plot_model(results, os.path.join(outdir, "figures"), self.figure_label)
        if text2workspace:
            self.text2workspace(outdir)
        print(f"Wrote {os.path.join(outdir, 'combined.txt')} (see commands.sh)")
        return summary

    def t2w_options(self):
        maps = []
        for ch in self.channels.values():
            for p in ch.processes.values():
                if p.signal:
                    poi = p.poi or self.default_poi
                    if poi not in self.pois:
                        raise ModelError(f"{ch.name}/{p.name}: unknown POI '{poi}' (add it with Model.add_poi)")
                    lo, hi = self.pois[poi]
                    maps.append(f"--PO 'map={ch.name}/{p.name}:{poi}[1,{lo:g},{hi:g}]'")
        return ["-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel"] + maps

    def text2workspace(self, outdir, card="combined.txt", output="model.root"):
        cmd = f"text2workspace.py {card} -o {output} " + " ".join(self.t2w_options())
        print(f"  (cd {outdir} && {cmd})")
        if shutil.which("text2workspace.py") is None:
            raise ModelError("text2workspace.py is not on PATH: set up Combine first")
        proc = subprocess.run(cmd, shell=True, cwd=outdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        with open(os.path.join(outdir, "text2workspace.log"), "w") as f:
            f.write(proc.stdout)
        if proc.returncode != 0:
            raise ModelError("text2workspace failed:\n" + "\n".join(proc.stdout.splitlines()[-20:]))

    def commands(self, results):
        envelope = any(r.envelopes for r in results.values())
        pois = list(self.pois)
        extra = (" --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams"
                 " --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints"
                 " --X-rtd MINIMIZER_multiMin_maskChannels=2") if envelope else ""

        def select(poi, value=None):
            """Options selecting one POI: the other POIs fixed to zero; value sets the POI for Asimov data."""
            others = [p for p in pois if p != poi]
            sets = ([f"{poi}={value:g}"] if value is not None else []) + [f"{p}=0" for p in others]
            out = f" --redefineSignalPOIs {poi}" if others else ""
            if sets:
                out += f" --setParameters {','.join(sets)}"
            if others:
                out += f" --freezeParameters {','.join(others)}"
            return out

        tag = lambda poi: f".{self.name}" + (f"_{poi}" if len(pois) > 1 else "")
        first = pois[0]
        lines = [
            "#!/bin/bash",
            f"# Combine commands for the model '{self.name}'. Run from this directory with Combine set up.",
            "# The cards are plain Combine datacards: edit or combine them as with any other.",
            "set -e",
            "",
            "# Build the workspace (signals mapped to their POIs)",
            f"text2workspace.py combined.txt -o model.root {' '.join(self.t2w_options())}",
            "",
            "# Expected upper limits (CLs, asymptotic); add --run blind for real data that must stay blind",
        ]
        if len(pois) > 1:
            lines.append("# one POI at a time, the other signals fixed to zero")
        for poi in pois:
            lines.append(f"combine -M AsymptoticLimits model.root -t -1{select(poi, 0.)} -n {tag(poi)}{extra}")
        lines += [
            "",
            f"# Expected significance for {first} = 1",
            f"combine -M Significance model.root -t -1{select(first, 1.)} -n {tag(first)}{extra}",
        ]
        if len(pois) > 1:
            lines += ["", f"# The commands below use {first}; change --redefineSignalPOIs/--setParameters for another POI"]
        lines += [
            "",
            "# Toy-based CLs limits (near-zero backgrounds), e.g. the median expected:",
            f"# combine -M HybridNew --LHCmode LHC-limits model.root -T 1000 --expectedFromGrid 0.5{select(first)}"
            f" -n {tag(first)}{extra}",
            "",
            "# Fit diagnostics with pre/post-fit shapes, then their figures (in figures/postfit)",
            f"# combine -M FitDiagnostics model.root -t -1{select(first, 1.)} --saveShapes --saveWithUncertainties"
            f" -n {tag(first)}{extra}",
            f"# python3 -m mu2ecombine.plots fitDiagnostics{tag(first)}.root",
            "",
            "# Likelihood scan of a POI",
            f"# combine -M MultiDimFit model.root -t -1{select(first, 1.)} --algo grid --points 50 -P {first}"
            f" -n {tag(first)}{extra}",
        ]
        groups = list(dict.fromkeys(list(self.groups) + [n.group for n in self.nuisances if n.group]))
        if groups:
            lines += ["", "# Freeze nuisance groups with, e.g., --freezeNuisanceGroups " + ",".join(groups)]
        return "\n".join(lines) + "\n"

    # -- cards ---------------------------------------------------------------------------------------------------------
    def _process_ids(self, results):
        ids, nsig, nbkg = {}, 0, 0
        for r in results.values():
            for p in r.processes:
                proc = r.channel.processes[p]
                if p in ids:
                    if (ids[p] <= 0) != proc.signal:
                        raise ModelError(f"process '{p}' is a signal in one channel and a background in another")
                    continue
                if proc.signal:
                    ids[p] = -nsig
                    nsig += 1
                else:
                    nbkg += 1
                    ids[p] = nbkg
        return ids

    def _check_counting(self, results):
        """A counting card has one bin per channel and only log-normal rate uncertainties."""
        if self.mcstats["enabled"] and self.mcstats["auto"]:
            raise ModelError("counting cards cannot use autoMCStats: use mc_stats(auto=False)")
        for r in results.values():
            ch = r.channel.name
            if len(r.channel.edges) != 2:
                raise ModelError(f"channel {ch} has {len(r.channel.edges) - 1} bins: a counting card needs one")
            if r.envelopes or r.pdfs:
                raise ModelError(f"channel {ch} has PDF or envelope processes: a counting card needs templates")
            shaped = sorted({name for name, kind, _, _ in r.effects if kind != "lnN"})
            if shaped:
                raise ModelError(f"channel {ch}: {', '.join(shaped)} are shape uncertainties; a counting card needs "
                                 "rate (lnN) uncertainties")

    def _card(self, results, ids, shapes=True):
        card = Datacard(f"Combine model '{self.name}' written by mu2ecombine; paths are relative to this file")
        effects = {}   # nuisance -> {"lnN": {(ch, p): text}, "shape": {(ch, p): text}, "ws": set}
        for r in results:
            ch = r.channel.name
            card.add_bin(ch, integral(r.data))
            if shapes:
                r.add_shapes(card)
            for p in r.processes:
                card.add_process(ch, p, ids[p], r.card_rates[p])
            for (name, kind, p, text) in r.effects:
                e = effects.setdefault(name, {"lnN": {}, "shape": {}, "ws": set()})
                if kind == "ws":
                    e["ws"].add((ch, p))
                else:
                    e[kind][(ch, p)] = text
        declared = set()
        for name, e in effects.items():
            if e["lnN"] and e["shape"]:
                raise ModelError(f"nuisance '{name}' is both a log-normal and a shape uncertainty: use different names")
            kind = "lnN" if e["lnN"] else "shape" if e["shape"] else None
            if kind == "lnN" and name in self.free:
                raise ModelError(f"nuisance '{name}' is a log-normal rate uncertainty: it cannot be unconstrained")
            if kind:
                declared.add(name)
                for (ch, p), text in e[kind].items():
                    card.add_systematic(name, "shapeU" if name in self.free else kind, ch, p, text)
        for name, e in effects.items():
            if name not in declared and e["ws"] and name not in self.params:
                card.extra.append(f"{name} flatParam" if name in self.free else f"{name} param 0 1")
        for name, procs, chans, init, (lo, hi) in self.rate_params:
            for r in results:
                for p in r.processes:
                    if matches(r.channel.name, chans) and matches(p, procs):
                        card.extra.append(f"{name} rateParam {r.channel.name} {p} {init:g} [{lo:g},{hi:g}]")
        params_in_card = set()
        for r in results:
            params_in_card |= r.free_params
        for name, (mean, sigma) in self.params.items():
            if name in params_in_card and name not in declared:
                card.extra.append(f"{name} param {mean:g} {sigma:g}")
        for r in results:
            for cat in r.envelopes.values():
                card.extra.append(f"{cat} discrete")
        groups = {}
        for n in self.nuisances:
            if n.group and n.name in effects:
                groups.setdefault(n.group, []).append(n.name)
        for g, names in self.groups.items():
            groups.setdefault(g, []).extend(x for x in names if x in effects or x in self.params)
        for g, names in groups.items():
            if names:
                card.extra.append(f"{g} group = {' '.join(dict.fromkeys(names))}")
        if self.mcstats["enabled"] and self.mcstats["auto"]:
            for r in results:
                if r.mode == "binned":
                    card.extra.append(f"{r.channel.name} autoMCStats {self.mcstats['threshold']:g}")
        if not results:
            raise ModelError("the model has no channels")
        return card

    # -- development -----------------------------------------------------------------------------------------------
    def fit_process(self, channel, process, figdir=None):
        """Build one process of one channel alone, to develop or debug its shape: a Fit is fit (and its variations
        refit), a Template is built, an Envelope is fit (after the channel's other processes, which it needs). Prints
        the result, draws the fit (or, for an envelope, the channel) to figdir if given, and returns the ChannelResult.
        Writes no model."""
        if channel not in self.channels:
            raise ModelError(f"no channel '{channel}' (channels: {', '.join(self.channels)})")
        ch = self.channels[channel]
        if process not in ch.processes:
            raise ModelError(f"channel {channel} has no process '{process}' (processes: {', '.join(ch.processes)})")
        r = _ChannelBuilder(self, ch, only=[process]).build()
        info = {"rate": r.card_rates.get(process), "mc_rel_error": r.mc_rel.get(process)}
        if process in r.fits:
            info["fit"] = r.fits[process].info
        if process in r.envelope_info:
            info["envelope"] = {k: v for k, v in r.envelope_info[process].items() if k != "tested"}
        print(f"{channel}/{process}: {json.dumps(info, indent=2, default=float)}")
        if figdir and ch.dim == 2:
            print("  (figures of 2D channels are not drawn by fit_process)")
        elif figdir:
            from . import plots as plotlib
            os.makedirs(figdir, exist_ok=True)
            plotlib._style()
            for logy, suffix in ((False, ""), (True, "_log")):
                if process in r.fits:
                    path = os.path.join(figdir, f"{channel}_{process}_fit{suffix}.png")
                    plotlib.plot_fit(r, process, r.fits[process], path, logy, self.figure_label)
                elif plotlib.is_template(r, process):
                    path = os.path.join(figdir, f"{channel}_{process}_template{suffix}.png")
                    plotlib.plot_template(r, process, path, logy, self.figure_label)
                else:
                    path = os.path.join(figdir, f"{channel}_{process}_prefit{suffix}.png")
                    plotlib.plot_channel(r, path, logy, self.figure_label)
                print(f"  figure: {path}")
        return r

    # -- summary -------------------------------------------------------------------------------------------------------
    def _summary(self, results, ids):
        out = {"name": self.name, "pois": {k: list(v) for k, v in self.pois.items()},
               "t2w_options": self.t2w_options(), "channels": {}}
        for name, r in results.items():
            out["channels"][name] = r.summary(ids)
        return out

    @staticmethod
    def print_summary(summary):
        for name, ch in summary["channels"].items():
            obs = ch["observable"]
            axes = obs.get("axes")
            where = (" x ".join(f"{a['name']} {a['range'][0]:g}-{a['range'][1]:g}, {a['bins']} bins" for a in axes)
                     + f", {obs['bins']} unrolled bins" if axes else
                     f"{obs['name']} {obs['range'][0]:g}-{obs['range'][1]:g}, {obs['bins']} bins")
            print(f"Channel {name} ({ch['mode']}, {where}), data: {ch['data']['source']} "
                  f"({ch['data']['integral']:.4g})")
            for p, info in ch["processes"].items():
                tag = " (signal)" if info["signal"] else ""
                unc = f" +- {info['mc_rel_error'] * 100:.1f}% (MC)" if info.get("mc_rel_error") else ""
                print(f"  {p:22s} {info['kind']:9s} {info['rate']:12.5g}{unc}{tag}")
            for p in ch.get("absorbed", []):
                print(f"  {p:22s} absorbed by the envelope")


def _shift_size(value, nominal):
    """A parameter shift per unit nuisance: a number, or 'x%' of the nominal value."""
    if isinstance(value, str) and value.strip().endswith("%"):
        return float(value.strip()[:-1]) / 100. * abs(nominal)
    return float(value)


class _ChannelBuilder:
    """Builds one channel: nominal shapes, variations, data, and the card entries."""

    def __init__(self, model, channel, only=None):
        self.model, self.channel = model, channel
        self.only = set(only) if only is not None else None   # build only these processes (see Model.fit_process)
        self.edges = channel.edges
        self.model_params = model.nuisance_names()

    def build(self):
        ch, model, edges = self.channel, self.model, self.edges
        r = ChannelResult(ch)
        if not ch.processes:
            raise ModelError(f"channel {ch.name} has no processes")
        envelopes = {n: p for n, p in ch.processes.items() if isinstance(p.shape, Envelope)
                     and (self.only is None or n in self.only)}
        # an envelope is fit with the channel's other backgrounds held fixed, so it needs all of them
        build_all = self.only is None or bool(envelopes)
        absorbed = set()
        for n, p in envelopes.items():
            absorbed |= {o for o in ch.processes if o not in envelopes and not ch.processes[o].signal
                         and any(matches(o, pat) for pat in p.shape.absorbs)}
        r.absorbed = sorted(absorbed)
        kinds = {n: p.shape.kind if not isinstance(p.shape, Envelope) else "envelope" for n, p in ch.processes.items()}
        workspace = model.workspace if ch.workspace is None else ch.workspace
        templates_only = all(kinds[n] == "template" for n in ch.processes if n not in absorbed)
        r.mode = "binned" if templates_only and not workspace else "workspace"
        shifts = [n for n in model.nuisances if n.kind == "shift"]
        self.thetas = {}   # nuisance name -> RooRealVar, one per channel workspace
        r.keep.append(self.thetas)
        expected = {}      # process -> expected events per bin (for data and envelope fits)

        for name, proc in ch.processes.items():
            shape = proc.shape
            if isinstance(shape, Envelope) or not (build_all or name in self.only):
                continue
            my_shifts = self._my_shifts(shifts, name)
            prefix = f"{ch.name}_{name}"
            try:
                if isinstance(shape, Product):
                    self._product(r, name, shape, my_shifts, prefix)
                elif isinstance(shape, Template):
                    self._template(r, name, shape, my_shifts)
                elif isinstance(shape, Fit) and shape.binned:
                    self._binned_fit(r, name, shape, my_shifts, prefix)
                elif isinstance(shape, Fit):
                    self._fit(r, name, shape, my_shifts, prefix)
                elif isinstance(shape, Pdf):
                    self._pdf(r, name, shape)
            except (ShapeError, HistError) as err:
                raise ModelError(f"{ch.name}/{name}: {err}") from None
            expected[name] = r.expected[name]

        # data: given, or the Asimov expectation of the backgrounds (absorbed processes included)
        bkg = [n for n in expected if not ch.processes[n].signal]
        if ch.data is not None:
            try:
                r.data = rebin_to(ch.data, edges, name="data_obs")
            except HistError as err:
                raise ModelError(f"{ch.name} data: {err}") from None
            r.data_source = "data"
        else:
            r.data = make(edges, sum((expected[n] for n in bkg), np.zeros(len(edges) - 1)), None, "data_obs")
            r.data_source = "asimov"

        for name, proc in envelopes.items():
            self._envelope(r, name, proc.shape, expected, absorbed)

        # processes kept in the card
        for name in ch.processes:
            if name in absorbed or name not in r.card_rates:
                continue
            if r.card_rates[name] <= 0 and kinds[name] != "envelope":
                warnings.warn(f"{ch.name}/{name}: no expected events in the channel range; dropped")
                continue
            r.processes.append(name)
        pdf_names = [r.pdfs[p].GetName() for p in r.processes if p in r.pdfs]
        if len(set(pdf_names)) != len(pdf_names):
            raise ModelError(f"channel {ch.name}: two processes use PDFs with the same name ({', '.join(pdf_names)})")

        # log-normal rate uncertainties
        for n in model.nuisances:
            if n.kind != "lnN":
                continue
            v = n.options["value"]
            text = kappa_text(v[1], v[0]) if isinstance(v, (tuple, list)) else kappa_text(float(v))
            for name in r.processes:
                if kinds[name] != "envelope" and n.applies(ch.name, name):
                    r.effects.append((n.name, "lnN", name, text))

        # MC statistics
        mc = model.mcstats
        if mc["enabled"]:
            for name in r.processes:
                rel = r.mc_rel.get(name, 0.)
                needs_lnn = not mc["auto"] or r.mode == "workspace" or name in r.smoothed
                if needs_lnn and rel > mc["min_rel"]:
                    r.effects.append((f"mcstat_{ch.name}_{name}", "lnN", name, kappa_text(1. + rel)))
        r.effects = [e for e in r.effects if e[2] in r.processes]
        return r

    # -- process kinds -------------------------------------------------------------------------------------------------
    def _my_shifts(self, shifts, name):
        """The shape shifts that act on a process of this channel (see Model.shape_shift: `axis`)."""
        ch, out = self.channel, []
        for s in shifts:
            if not s.applies(ch.name, name):
                continue
            axis = s.options.get("axis")
            if ch.dim == 2 and axis is None:
                raise ModelError(f"shape_shift '{s.name}' applies to {ch.name}/{name} in a 2D channel: give "
                                 f"axis= ({', '.join(a.GetName() for a in ch.axes)})")
            if axis is not None and axis not in [a.GetName() for a in ch.axes]:
                continue
            out.append(s)
        return out

    def _shift_source(self, shape, s, sign):
        """The shape's source histogram (unrolled TH1D for a 2D Template) transformed by a shape shift."""
        o = s.options
        raw = getattr(shape, "raw", None)
        if raw is None:
            return transformed(shape.source(), sign * o["shift"], sign * o["stretch"], o["about"])
        axis = [a.GetName() for a in self.channel.axes].index(o["axis"])
        return unroll(rebin2_to(transformed2(raw, axis, sign * o["shift"], sign * o["stretch"], o["about"]),
                                *self.channel.axis_edges))

    def _template_parts(self, name, shape, shifts, edges):
        """(nominal, mc_rel, [(nuisance, up, down)]) of a Template in `edges`."""
        nominal = shape.nominal(edges)
        if np.any(arrays(nominal)[1] < 0):
            warnings.warn(f"{self.channel.name}/{name}: negative bins set to zero")
            clip_negative(nominal)
        variations = [(vname, *shape.varied(vname, edges)) for vname in shape.variations]
        for s in shifts:
            varied = []
            for sign in (1., -1.):
                src = self._shift_source(shape, s, sign)
                varied.append(in_channel(src, edges, "Template") if shape.shape_hist is None
                              else shape.nominal(edges, source=src))
            variations.append((s.name, varied[0], varied[1]))
        return nominal, rel_error(in_channel(shape.hist, edges, "Template")), variations

    def _template(self, r, name, shape, shifts):
        nominal, mc_rel, variations = self._template_parts(name, shape, shifts, self.edges)
        r.add_template(name, name, nominal)
        r.expected[name] = arrays(nominal)[1]
        r.card_rates[name] = integral(nominal)
        r.mc_rel[name] = mc_rel
        for vname, up, down in variations:
            self._add_variation(r, name, vname, nominal, up, down)

    def _add_variation(self, r, name, vname, nominal, up, down):
        if integral(up) <= 0 or integral(down) <= 0:
            warnings.warn(f"{self.channel.name}/{name}: variation {vname} is empty; it has no effect on this process")
            up, down = scaled(nominal, 1.), scaled(nominal, 1.)
        r.add_template(name, f"{name}_{vname}Up", clip_negative(up))
        r.add_template(name, f"{name}_{vname}Down", clip_negative(down))
        r.effects.append((vname, "shape", name, "1"))

    def _binned_fit_parts(self, shape, obs, edges, prefix, shifts, name):
        """Fit `shape` on obs and use it binned: (nominal, rate, mc_rel, [(nuisance, up, down)]), all in `edges`."""
        from .pdfs import bin_fractions
        shape.build(obs, f"{prefix}_smooth", edges, self.model_params)
        rate = shape.rate(edges)
        nominal = make(edges, bin_fractions(shape.pdf, obs, edges) * rate, None, name)
        ref = integral(shape.fit_hist(edges))

        def template(values, k=1.):
            nominal_values = shape.values()
            shape.set_values(values)
            out = make(edges, bin_fractions(shape.pdf, obs, edges) * rate * k, None, name)
            shape.set_values(nominal_values)
            return out

        def from_hist(h, mode):
            values, _ = shape.refit(h, edges, obs)
            return template(values, integral(in_channel(h, edges, "Fit")) / ref if mode == "shape" else 1.)

        variations = []
        for vname, (up, down, mode) in shape.variations.items():
            if mode == "rate":
                variations.append((vname, scaled(nominal, _ratio(up, ref, edges)),
                                   scaled(nominal, _ratio(down, ref, edges))))
            else:
                variations.append((vname, from_hist(up, mode), from_hist(down, mode)))
        for s in shifts:
            varied = [from_hist(self._shift_source(shape, s, sign), "shape") for sign in (1., -1.)]
            variations.append((s.name, varied[0], varied[1]))
        for nuis, params in shape.nuisances.items():
            nominal_values = shape.values()
            varied = []
            for sign in (1., -1.):
                values = dict(nominal_values)
                for pname, size in params.items():
                    p = shape.param(pname, f"{prefix}_smooth")
                    values[p.GetName()] += sign * _shift_size(size, nominal_values[p.GetName()])
                varied.append(template(values))
            variations.append((nuis, varied[0], varied[1]))
        return nominal, rate, shape.rate_rel_error(edges), variations

    def _binned_fit(self, r, name, shape, shifts, prefix):
        nominal, rate, mc_rel, variations = self._binned_fit_parts(shape, self.channel.obs, self.edges, prefix, shifts,
                                                                   name)
        r.fits[name] = shape
        r.add_template(name, name, nominal)
        r.expected[name] = arrays(nominal)[1]
        r.card_rates[name] = rate
        r.mc_rel[name] = mc_rel
        r.smoothed.add(name)
        for vname, up, down in variations:
            self._add_variation(r, name, vname, nominal, up, down)

    def _product(self, r, name, shape, shifts, prefix):
        """A Product in a 2D channel: each component is built in 1D on its axis, the template is
        rate * outer(fx, fy) unrolled, and a variation of a component multiplies the other axis's nominal fractions
        (the variations of both components of a shared nuisance are applied together)."""
        ch = self.channel
        names = [a.GetName() for a in ch.axes]
        parts = []   # per axis: (nominal values, integral, mc_rel, {nuisance: (up values, down values)})
        for a, comp in enumerate(shape.components):
            obs, edges, other = ch.axes[a], ch.axis_edges[a], ch.axes[1 - a]
            cshifts = [s for s in shifts if s.options["axis"] == names[a]]
            if isinstance(comp, Template):
                nominal, mc_rel, variations = self._template_parts(name, comp, cshifts, edges)
            else:
                if not isinstance(comp.function, str) and comp.function.dependsOn(other):
                    raise ShapeError(f"the {'xy'[a]} component PDF '{comp.function.GetName()}' depends on the other "
                                     f"observable '{other.GetName()}'; a Fit component is on '{obs.GetName()}' only")
                nominal, _, mc_rel, variations = self._binned_fit_parts(comp, obs, edges, f"{prefix}_{'xy'[a]}",
                                                                        cshifts, name)
                r.components.setdefault(name, {})[names[a]] = comp
            total = integral(nominal)
            if total <= 0:
                raise ShapeError(f"the {'xy'[a]} component is empty in the channel range")
            parts.append((arrays(nominal)[1], total, mc_rel,
                          {v: (arrays(up)[1], arrays(dn)[1]) for v, up, dn in variations}))
        (vx, nx, mc_x, varx), (vy, ny, _, vary) = parts
        ex, ey = ch.axis_edges
        if shape.rate_source is None:
            rate, mc_rel = nx, mc_x
        elif isinstance(shape.rate_source, (int, float)):
            rate, mc_rel = float(shape.rate_source), 0.
        else:
            src = shape.rate_source
            h = unroll(rebin2_to(src, ex, ey)) if is_th2(src) else in_channel(src, ex, "Product rate")
            rate, mc_rel = integral(h), rel_error(h)

        def build2(fx, fy, label):
            return make(self.edges, (rate * np.outer(fx / nx, fy / ny)).T.ravel(), None, label)

        nominal = build2(vx, vy, name)
        r.add_template(name, name, nominal)
        r.expected[name] = arrays(nominal)[1]
        r.card_rates[name] = integral(nominal)
        r.mc_rel[name] = mc_rel
        r.smoothed.add(name)   # the template has no bin errors: the MC statistics are the log-normal
        for vname in dict.fromkeys(list(varx) + list(vary)):
            ux, dx = varx.get(vname, (vx, vx))
            uy, dy = vary.get(vname, (vy, vy))
            self._add_variation(r, name, vname, nominal, build2(ux, uy, name), build2(dx, dy, name))

    def _fit(self, r, name, shape, shifts, prefix):
        from .pdfs import add_nuisances, bin_fractions
        edges, obs = self.edges, self.channel.obs
        shape.build(obs, prefix, edges, self.model_params)
        r.fits[name] = shape
        rate = shape.rate(edges)
        nominal_values = shape.values()
        param_shifts = {}   # RooRealVar -> {nuisance: (up, down)}
        kappas = {}         # nuisance -> (down, up)
        ref = integral(shape.fit_hist(edges))

        def set_param_shifts(nuis, up_hist, down_hist):
            up_vals, s_up = shape.refit(up_hist, edges, obs)
            dn_vals, s_dn = shape.refit(down_hist, edges, obs)
            if s_up not in (0, 1) or s_dn not in (0, 1):
                warnings.warn(f"{prefix}: refit for {nuis} has status {s_up}/{s_dn}")
            for p in shape.floating:
                n = p.GetName()
                param_shifts.setdefault(p, {})[nuis] = (up_vals[n] - nominal_values[n], nominal_values[n] - dn_vals[n])

        for vname, (up, down, mode) in shape.variations.items():
            if mode != "rate":
                set_param_shifts(vname, up, down)
            if mode != "shape_only":
                kappas[vname] = (_ratio(down, ref, edges), _ratio(up, ref, edges))
        for s in shifts:
            hs = [transformed(shape.source(), sign * s.options["shift"], sign * s.options["stretch"],
                              s.options["about"]) for sign in (1., -1.)]
            set_param_shifts(s.name, hs[0], hs[1])
            if shape.rate_source is None and ref > 0:
                kappas[s.name] = (_ratio(hs[1], ref, edges), _ratio(hs[0], ref, edges))
        for nuis, params in shape.nuisances.items():
            for pname, size in params.items():
                p = shape.param(pname, prefix)
                d = _shift_size(size, nominal_values[p.GetName()])
                param_shifts.setdefault(p, {})[nuis] = (d, d)
        shape.freeze()
        if shape.pdf.GetName() == prefix:
            shape.pdf.SetName(f"{prefix}_nominal")
        pdf, keep = add_nuisances(shape.pdf, prefix, param_shifts, self.thetas)
        r.keep += keep
        r.pdfs[name] = pdf
        r.norms[name] = self._norm(r, pdf.GetName(), kappas)
        r.card_rates[name] = rate
        r.mc_rel[name] = shape.rate_rel_error(edges)
        r.expected[name] = bin_fractions(pdf, obs, edges) * rate
        nuisances = [n for s in param_shifts.values() for n in s] + list(kappas) + [p.GetName() for p in shape.hold]
        for nuis in dict.fromkeys(nuisances):
            r.effects.append((nuis, "ws", name, None))
        r.free_params |= {p.GetName() for p in shape.hold}

    def _pdf(self, r, name, shape):
        import ROOT
        from .pdfs import bin_fractions
        edges, obs = self.edges, self.channel.obs
        if not shape.pdf.dependsOn(obs):
            raise ModelError(f"{self.channel.name}/{name}: the PDF does not depend on the channel observable "
                             f"'{obs.GetName()}'")
        kappas = {vn: shape.kappas(vn, edges) for vn in shape.variations}
        r.pdfs[name] = shape.pdf
        r.norms[name] = self._norm(r, shape.pdf.GetName(), kappas, shape.norm)
        rate = shape.rate(edges)
        r.card_rates[name] = rate
        scale = shape.norm.getVal() if shape.norm is not None else 1.
        r.expected[name] = bin_fractions(shape.pdf, obs, edges) * rate * scale
        for vn in kappas:
            r.effects.append((vn, "ws", name, None))
        for p in shape.pdf.getParameters(ROOT.RooArgSet(obs)):
            r.free_params.add(p.GetName())
        if shape.norm is not None:
            for p in shape.norm.getVariables():
                r.free_params.add(p.GetName())

    def _norm(self, r, pdf_name, kappas, extra=None):
        """<pdf>_norm: the rate effects of the nuisances (asymmetric log-normal), times an optional extra factor."""
        import ROOT
        kappas = {k: v for k, v in kappas.items() if abs(v[0] - 1.) > 1e-6 or abs(v[1] - 1.) > 1e-6}
        if not kappas and extra is None:
            return None
        args = ROOT.RooArgList()
        terms = []
        for nuis, (down, up) in kappas.items():
            if down <= 0 or up <= 0:
                raise ModelError(f"{pdf_name}: non-positive rate factor for {nuis} ({down:g}/{up:g})")
            theta = self.thetas.get(nuis)
            if theta is None:
                theta = self.thetas[nuis] = ROOT.RooRealVar(nuis, nuis, 0., -7., 7.)
            i = args.getSize()
            args.add(theta)
            terms.append(f"((@{i}>0)*TMath::Power({up:.6g},@{i}) + (@{i}<=0)*TMath::Power({down:.6g},-@{i}))")
        if extra is not None:
            terms.append(f"@{args.getSize()}")
            args.add(extra)
        norm = ROOT.RooFormulaVar(f"{pdf_name}_norm", "", "*".join(terms), args)
        r.keep += [args, norm] + list(self.thetas.values())
        return norm

    def _envelope(self, r, name, env, expected, absorbed):
        import ROOT
        from .combine import load_library
        ch = self.channel
        load_library()
        others = [n for n in expected if n not in absorbed and not ch.processes[n].signal]
        fixed = sum((expected[n] for n in others), np.zeros(len(self.edges) - 1)) if others else None
        print(f"  {ch.name}/{name}: envelope fit to the {r.data_source} outside the blind window {ch.blind}")
        multi, cat, rate = env.build(f"{ch.name}_{name}", ch.obs, r.data, fixed, ch.blind)
        norm = ROOT.RooRealVar(f"{ch.name}_{name}_norm", f"{name} yield", rate, 0., max(10. * rate, 10.))
        r.keep += [env, norm]
        r.pdfs[name] = multi
        r.norms[name] = norm
        r.envelopes[name] = cat.GetName()
        r.envelope_info[name] = {"rate": rate, "best": env.best, "functions": [f.as_dict() for f in env.accepted],
                                 "tested": [f.as_dict() for f in env.fits]}
        r.card_rates[name] = 1.
        r.expected[name] = env.accepted[env.best].expected - (fixed if fixed is not None else 0.)


def _ratio(h, reference, edges):
    """Integral of h in the channel binning over a reference integral."""
    return integral(in_channel(h, edges, "variation")) / reference if reference > 0 else 1.


class ChannelResult:
    """A built channel: templates (TH1D), PDFs, data, and card entries."""

    def __init__(self, channel):
        self.channel = channel
        self.mode = "binned"
        self.templates = {}     # name -> TH1D (processes and <process>_<nuisance>Up/Down)
        self.owner = {}         # template name -> process
        self.pdfs = {}          # process -> RooAbsPdf
        self.norms = {}         # process -> RooAbsReal or None
        self.fits = {}          # process -> Fit
        self.components = {}    # Product process -> {axis observable name: Fit} of its fitted components (QA fits)
        self.expected = {}      # process -> expected events per bin (numpy)
        self.card_rates = {}
        self.mc_rel = {}
        self.smoothed = set()
        self.processes = []     # in the card
        self.absorbed = []
        self.effects = []       # (nuisance, kind 'lnN'|'shape'|'ws', process, text)
        self.envelopes = {}     # process -> category name
        self.envelope_info = {}
        self.free_params = set()
        self.data = None        # TH1D
        self.data_source = ""
        self.keep = []

    def add_template(self, process, key, hist):
        self.templates[key] = hist
        self.owner[key] = process

    @property
    def shapes_file(self):
        return f"shapes/{self.channel.name}.root"

    def write_shapes(self, outdir):
        import ROOT
        path = os.path.join(outdir, self.shapes_file)
        ch = self.channel.name
        if self.mode == "binned":
            f = ROOT.TFile(path, "RECREATE")
            for key, h in [("data_obs", self.data)] + list(self.templates.items()):
                if key == "data_obs" or self.owner[key] in self.processes:
                    h.SetName(key)
                    f.WriteTObject(h, key)
            f.Close()
            return
        obs = self.channel.obs
        ws = ROOT.RooWorkspace("w", "w")
        imp = getattr(ws, "import")
        imp(ROOT.RooDataHist(f"{ch}_data_obs", "", ROOT.RooArgList(obs), self.data), ROOT.RooFit.Silence())
        for key, h in self.templates.items():
            if self.owner[key] in self.processes:
                imp(ROOT.RooDataHist(f"{ch}_{key}", "", ROOT.RooArgList(obs), h), ROOT.RooFit.Silence())
        for name, pdf in self.pdfs.items():
            if name not in self.processes:
                continue
            imp(pdf, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
            norm = self.norms.get(name)
            if norm is not None:
                imp(norm, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
        ws.writeToFile(path)

    def add_shapes(self, card):
        ch, f = self.channel.name, self.shapes_file
        if self.mode == "binned":
            card.add_shapes("data_obs", ch, f, "data_obs")
            card.add_shapes("*", ch, f, "$PROCESS", "$PROCESS_$SYSTEMATIC")
            return
        card.add_shapes("data_obs", ch, f, f"w:{ch}_data_obs")
        for p in self.processes:
            if p in self.pdfs:
                card.add_shapes(p, ch, f, f"w:{self.pdfs[p].GetName()}")
            else:
                card.add_shapes(p, ch, f, f"w:{ch}_{p}", f"w:{ch}_{p}_$SYSTEMATIC")

    def summary(self, ids):
        ch = self.channel
        out = {"mode": self.mode, "title": ch.title, "blind": ch.blind,
               "observable": {"name": ch.obs.GetName(), "title": ch.obs.GetTitle(),
                              "range": [float(ch.edges[0]), float(ch.edges[-1])], "bins": len(ch.edges) - 1,
                              "edges": [float(e) for e in ch.edges]},
               "data": {"source": self.data_source, "integral": integral(self.data)},
               "processes": {}, "absorbed": self.absorbed}
        if ch.dim == 2:
            out["observable"]["axes"] = [{"name": x.GetName(), "title": x.GetTitle(),
                                          "range": [float(e[0]), float(e[-1])], "bins": len(e) - 1,
                                          "edges": [float(v) for v in e]} for x, e in zip(ch.axes, ch.axis_edges)]
            out["observable"]["unroll"] = "k = iy * nx + ix"
        for p in self.processes:
            proc = ch.processes[p]
            kind = "envelope" if p in self.envelopes else proc.shape.kind
            if isinstance(proc.shape, Fit):
                kind = "fit-binned" if proc.shape.binned else "fit"
            if isinstance(proc.shape, Product):
                kind = "product"
            info = {"kind": kind, "signal": proc.signal, "id": ids[p], "title": proc.title,
                    "rate": float(self.expected[p].sum()) if p in self.envelopes else self.card_rates[p],
                    "nuisances": sorted({e[0] for e in self.effects if e[2] == p})}
            if self.mc_rel.get(p):
                info["mc_rel_error"] = self.mc_rel[p]
            if p in self.fits:
                info["fit"] = self.fits[p].info
            if p in self.components:
                info["components"] = {axis: fit.info for axis, fit in self.components[p].items()}
            if p in self.envelope_info:
                info["envelope"] = self.envelope_info[p]
            out["processes"][p] = info
        return out
