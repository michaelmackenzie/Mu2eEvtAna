"""How a process gets its shape and rate. Histograms are ROOT TH1s (or lists of TH1s, summed); PDFs are RooFit PDFs.

  th1                                  binned template; the rate is the histogram integral in the channel range
  Template(th1, shape=[cr1, cr2])      rate from th1, shape from other histograms (control regions, merged sets)
  Fit("dscb", th1)                     a predefined function (pdfs.py) or any RooAbsPdf fit to an MC histogram, then
                                       frozen; rate from the histogram
  Fit(..., binned=True)                the fitted function, binned, used as a template (smoothing sparse MC)
  Pdf(roo_pdf, rate)                   a RooAbsPdf used as is; its free parameters float in the Combine fit
  Envelope(...)                        data-driven discrete profiling over function families (envelope.py)

  Product(x_shape, y_shape)            2D (e.g. p x t) shape from two 1D shapes (TH1, Template, or Fit), each on its own
                                       observable; for a channel on two observables (see Model.channel)

Histograms may come in any binning whose edges include the channel's bin edges; they are rebinned automatically.

Uncertainties from histogram variations, variations={name: (up, down)} or (up, down, mode):
  mode "shape"      (default) shape and rate change: Combine shape morphing for templates; for fitted PDFs, the
                    function is refit to up/down and its parameters move with the nuisance, plus the rate change
  mode "shape_only" the variations are normalized to the nominal rate
  mode "rate"       only the change in integral (an asymmetric log-normal rate uncertainty)
For a Pdf process, (up, down) may also be rate factors (kappa_up, kappa_down).
A nuisance with the same name in several processes/channels is one correlated parameter.
"""

from .hist import (arrays, as_hist, as_th1, integral, is_th1, is_th2, make, rebin2_to, rebin_to, rel_error,
                   unroll)

MODES = ("shape", "shape_only", "rate")


class ShapeError(ValueError):
    pass


def as_variations(variations, numbers=False):
    """{name: (up, down, mode)} with histograms converted to detached TH1Ds (or kept as numbers if allowed)."""
    out = {}
    for name, v in (variations or {}).items():
        if not isinstance(v, (tuple, list)) or len(v) not in (2, 3):
            raise ShapeError(f"variation '{name}': give (up, down) or (up, down, mode)")
        mode = v[2] if len(v) == 3 else "shape"
        if mode not in MODES:
            raise ShapeError(f"variation '{name}': mode '{mode}', use one of {', '.join(MODES)}")
        if numbers and all(isinstance(x, (int, float)) for x in v[:2]):
            out[name] = (float(v[0]), float(v[1]), mode)
        else:
            out[name] = (as_hist(v[0], f"variation {name} up"), as_hist(v[1], f"variation {name} down"), mode)
    return out


def in_channel(h, edges, what):
    try:
        return rebin_to(h, edges)
    except ValueError as err:
        raise ShapeError(f"{what}: {err}") from None


def _ratio(a, b):
    return a / b if b > 0 else 1.


def scaled(h, factor):
    out = make(*arrays(h), name=h.GetName())
    out.Scale(factor)
    return out


class Template:
    """A binned template: rate from `hist`, shape from `shape` (default `hist`). The histograms are TH1s, or TH2s for
    a channel on two observables (all of the same dimension). In a 2D channel the Template is converted to unrolled
    TH1Ds (see `unrolled`) and `raw` keeps the source TH2 for shape shifts along an axis."""

    kind = "template"

    def __init__(self, hist, shape=None, variations=None):
        self.hist = as_hist(hist, "Template")
        self.shape_hist = as_hist(shape, "Template shape") if shape is not None else None
        self.variations = as_variations(variations)
        self.raw = None
        dims = {2 if is_th2(h) else 1 for h in [self.hist, self.shape_hist] if h is not None}
        dims |= {2 if is_th2(h) else 1 for v in self.variations.values() for h in v[:2]}
        if len(dims) != 1:
            raise ShapeError("Template: the histograms must be all TH1 or all TH2")
        self.dim = dims.pop()

    def unrolled(self, xedges, yedges):
        """The 2D Template as a Template of unrolled TH1Ds in (xedges, yedges), with `raw` the source TH2."""
        if self.dim != 2:
            raise ShapeError("Template: expected 2D histograms for a 2D channel")

        def flat(h, what):
            try:
                return unroll(rebin2_to(h, xedges, yedges))
            except ValueError as err:
                raise ShapeError(f"{what}: {err}") from None
        variations = {n: (flat(up, f"variation {n} up"), flat(dn, f"variation {n} down"), mode)
                      for n, (up, dn, mode) in self.variations.items()}
        out = Template(flat(self.hist, "Template"), flat(self.shape_hist, "Template shape")
                       if self.shape_hist is not None else None, variations)
        out.raw = self.source()
        return out

    def source(self):
        """The histogram that defines the shape (before rebinning)."""
        return self.shape_hist if self.shape_hist is not None else self.hist

    def rate(self, edges):
        return integral(in_channel(self.hist, edges, "Template"))

    def nominal(self, edges, source=None):
        """The template in the channel binning; `source` replaces the shape histogram (e.g. a shifted copy)."""
        if self.shape_hist is None and source is None:
            return in_channel(self.hist, edges, "Template")
        shape = in_channel(self.source() if source is None else source, edges, "Template shape")
        total = integral(shape)
        if total <= 0:
            raise ShapeError("Template: the shape histogram is empty in the channel range")
        return scaled(shape, self.rate(edges) / total)   # keeps the shape's relative bin errors

    def varied(self, name, edges):
        """(up, down) templates of a histogram variation in the channel binning."""
        up, down, mode = self.variations[name]
        nominal = self.nominal(edges)
        up = in_channel(up, edges, f"variation {name} up")
        down = in_channel(down, edges, f"variation {name} down")
        reference = self.rate(edges)
        if mode == "rate" or self.shape_hist is not None:
            # with a separate shape source, the variations (of the rate histogram) give the rate effect only
            return scaled(nominal, _ratio(integral(up), reference)), scaled(nominal, _ratio(integral(down), reference))
        if mode == "shape_only":
            n = integral(nominal)
            return scaled(up, _ratio(n, integral(up))), scaled(down, _ratio(n, integral(down)))
        return up, down


class Fit:
    """A PDF fit to an MC histogram and frozen.

    function:  a predefined name ('dscb', 'exponential', ...; see pdfs.py) or a RooAbsPdf built on the channel
               observable. Its non-constant parameters are fit, except those declared with Model.param (nuisances
               already in your PDF), which stay at their values during the fit and float in Combine.
    hist:      histogram(s) the function is fit to (e.g. merged control regions); also the rate unless `rate` is given
    rate:      a number or a histogram (its integral in the channel range)
    order:     polynomial order for chebychev / bernstein / exp_poly / endpoint
    params:    {name: value (fixed) | (init, lo, hi)} for a predefined function
    fit_range: (lo, hi) to fit in, default the channel range
    binned:    use the fitted function as a binned template (no MC bin fluctuations); the rate's MC uncertainty is kept
               as a log-normal
    nuisances: {nuisance: {parameter: shift}} shifts of fitted parameters per unit nuisance; a shift is a number
               (observable units) or a string 'x%' (relative to the fitted value), e.g. {"pscale": {"mean": 0.05}}.
               Parameters are named as in pdfs.py for predefined functions, or by their RooRealVar names.
    """

    def __init__(self, function, hist, rate=None, order=None, params=None, fit_range=None, binned=False,
                 variations=None, nuisances=None):
        if not isinstance(function, str) and not _is_pdf(function):
            raise ShapeError(f"Fit: give a function name or a RooAbsPdf, not {type(function).__name__}")
        self.function = function
        self.hist = as_th1(hist, "Fit")
        self.rate_source = rate if rate is None or isinstance(rate, (int, float)) else as_th1(rate, "Fit rate")
        self.order, self.params, self.fit_range = order, params, fit_range
        self.binned = binned
        self.variations = as_variations(variations)
        self.nuisances = {n: dict(p) for n, p in (nuisances or {}).items()}
        self.kind = "template" if binned else "fit"
        self.pdf = None
        self.floating = []
        self.hold = []
        self.keep = []
        self.info = {}

    def source(self):
        return self.hist

    def rate(self, edges):
        if self.rate_source is None:
            return integral(in_channel(self.hist, edges, "Fit"))
        if isinstance(self.rate_source, (int, float)):
            return float(self.rate_source)
        return integral(in_channel(self.rate_source, edges, "Fit rate"))

    def rate_rel_error(self, edges):
        if isinstance(self.rate_source, (int, float)):
            return 0.
        return rel_error(in_channel(self.hist if self.rate_source is None else self.rate_source, edges, "Fit"))

    def fit_hist(self, edges, source=None):
        """The histogram the function is fit to, in the channel binning."""
        return in_channel(self.hist if source is None else source, edges, "Fit")

    def build(self, obs, name, edges, model_params=()):
        """Build (for a name) and fit the PDF on obs. Returns the RooAbsPdf."""
        from .pdfs import fit_pdf, floating_params, make_pdf, quiet_roofit
        quiet_roofit()
        h = self.fit_hist(edges)
        if integral(h) <= 0:
            raise ShapeError(f"{name}: the histogram to fit is empty in the channel range")
        if isinstance(self.function, str):
            self.pdf, self.keep = make_pdf(self.function, obs, name, h, self.order, self.params)
        else:
            if not self.function.dependsOn(obs):
                raise ShapeError(f"{name}: the PDF '{self.function.GetName()}' does not depend on the channel "
                                 f"observable '{obs.GetName()}'")
            if self.params:
                raise ShapeError(f"{name}: params= applies to predefined functions; set your RooRealVars directly")
            self.pdf = self.function
        self.floating = floating_params(self.pdf, obs, exclude=model_params)
        self.hold = [p for p in floating_params(self.pdf, obs) if p.GetName() in model_params]
        status = fit_pdf(self.pdf, obs, h, self.floating, self.fit_range, self.hold)
        self.info = {"function": self.function if isinstance(self.function, str) else self.function.GetName(),
                     "fit_status": status, "params": self.values()}
        if status not in (0, 1):
            print(f"  warning: {name}: fit status {status}")
        return self.pdf

    def values(self):
        return {p.GetName(): p.getVal() for p in self.floating}

    def set_values(self, values):
        for p in self.floating:
            p.setVal(values[p.GetName()])

    def param(self, name, prefix):
        """The fitted RooRealVar called `name` (or <prefix>_<name> for a predefined function)."""
        for p in self.floating:
            if p.GetName() in (name, f"{prefix}_{name}"):
                return p
        raise ShapeError(f"{prefix}: no fitted parameter '{name}' (fitted: "
                         f"{', '.join(p.GetName() for p in self.floating) or 'none'})")

    def refit(self, hist, edges, obs):
        """Parameter values from a fit to another histogram (starting from the nominal fit), restoring the nominal."""
        from .pdfs import fit_pdf
        nominal = self.values()
        status = fit_pdf(self.pdf, obs, self.fit_hist(edges, hist), self.floating, self.fit_range, self.hold)
        values = self.values()
        self.set_values(nominal)
        return values, status

    def freeze(self):
        for p in self.floating:
            p.setConstant(True)


class Pdf:
    """A RooAbsPdf you built on the channel observable, used as is.

    rate:        expected events (a number, or a histogram whose integral in the channel range is used)
    norm:        optional RooAbsReal multiplying the rate (e.g. a function of your own parameters)
    variations:  {name: (kappa_up, kappa_down)} rate factors, or (up, down) histograms (rate effect only)
    Free parameters of the PDF float in the fit; constrain them with Model.param(name, mean, sigma).
    """

    kind = "pdf"

    def __init__(self, pdf, rate, norm=None, variations=None):
        if not _is_pdf(pdf):
            raise ShapeError(f"Pdf: expected a RooAbsPdf, got {type(pdf).__name__}")
        self.pdf = pdf
        self.rate_source = rate if isinstance(rate, (int, float)) else as_th1(rate, "Pdf rate")
        self.norm = norm
        self.variations = as_variations(variations, numbers=True)

    def rate(self, edges):
        if isinstance(self.rate_source, (int, float)):
            return float(self.rate_source)
        return integral(in_channel(self.rate_source, edges, "Pdf rate"))

    def kappas(self, name, edges):
        """(kappa_down, kappa_up) of a variation."""
        up, down, _ = self.variations[name]
        if isinstance(up, float):
            return down, up
        rate = self.rate(edges)
        return (_ratio(integral(in_channel(down, edges, f"variation {name} down")), rate),
                _ratio(integral(in_channel(up, edges, f"variation {name} up")), rate))


class Product:
    """A factorized 2D shape for a channel on two observables (x, y): the template is rate x outer(fx, fy), with fx and
    fy the normalized bin fractions of two 1D shapes, each built in 1D on its own observable and binning, and unrolled.

    x_shape, y_shape:  a TH1 (or list of TH1s), Template, or Fit on that axis's RooRealVar. A Fit is fit in 1D and
                       used binned (a PDF component must depend on its own axis observable only); the fitted
                       components are kept in ChannelResult.components for the fit QA figures.
    rate:              a number, or a TH1/TH2 whose integral in the channel range is used; default the x component's
                       rate. The MC relative uncertainty is that of the rate source (the x component by default).
    Variations: those of each component (Template variations, Fit variations and `nuisances=`, and Model.shape_shift
    along the component's axis). A nuisance is varied in every component that has it, the other axis staying
    nominal; the rate effect of a variation (modes "shape" and "rate") stays with the component carrying it.
    """

    kind = "template"

    def __init__(self, x_shape, y_shape, rate=None):
        self.x, self.y = (self._component(c, a) for c, a in ((x_shape, "x"), (y_shape, "y")))
        if rate is None or isinstance(rate, (int, float)):
            self.rate_source = rate
        else:
            self.rate_source = as_hist(rate, "Product rate")

    @staticmethod
    def _component(obj, axis):
        if is_th1(obj) or (isinstance(obj, (list, tuple)) and obj and all(is_th1(o) for o in obj)):
            return Template(obj)
        if isinstance(obj, (Template, Fit)):
            if isinstance(obj, Template) and obj.dim != 1:
                raise ShapeError(f"Product: the {axis} component is a 2D Template; give 1D histograms")
            return obj
        raise ShapeError(f"Product: the {axis} component must be a TH1, a list of TH1s, a Template, or a Fit, got "
                         f"{type(obj).__name__}")

    @property
    def components(self):
        return (self.x, self.y)

    def nuisance_names(self):
        return {n for c in self.components for n in list(c.variations) + list(getattr(c, "nuisances", {}))}


def _is_pdf(obj):
    try:
        return bool(obj.InheritsFrom("RooAbsPdf"))
    except AttributeError:
        return False


def as_shape(obj, what):
    """A Template from a TH1 or TH2 or a list of them; other shapes pass through."""
    if is_th1(obj) or is_th2(obj) or isinstance(obj, (list, tuple)):
        return Template(obj)
    raise ShapeError(f"{what}: expected a TH1, a list of TH1s, Template, Fit, Pdf, Envelope, or Product, got "
                     f"{type(obj).__name__}")
