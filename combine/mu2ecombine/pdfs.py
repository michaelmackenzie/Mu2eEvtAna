"""Predefined RooFit PDFs, fits to histograms, and nuisance-dependent parameters.

Predefined functions (name: parameters), each built as an ordinary RooFit PDF:
  gauss:        mean, sigma                      (a RooCrystalBall with far tails: Combine treats a top-level
                                                  RooGaussian process PDF as a constraint term)
  cb:           mean, sigma, alpha, n            RooCrystalBall (low-side tail; negative alpha for a high-side tail)
  dscb:         mean, sigma, sigmaR, alphaL, nL, alphaR, nR    RooCrystalBall, double-sided
  exponential:  slope                            RooExponential
  power:        index                            RooPowerLaw (x^index, needs x > 0)
  chebychev:    c1..c<order>                     RooChebychev
  bernstein:    c0..c<order>                     RooBernstein (c0 fixed to 1)
  exp_poly:     c1..c<order>                     exp(sum c_i t^i), t the observable scaled to [-1, 1]
  endpoint:     e0, a, c1..c<order>              (e0 - x)^a (1 + sum c_i t^i), e0 above the range
  uniform:      (none)                           RooUniform

Parameter settings (Fit(params=...)): a number fixes the parameter, (init, lo, hi) floats it. Defaults are set from the
range and the histogram being fit, so most fits need none. Any RooAbsPdf can be used instead of a name.
"""

import numpy as np

from .hist import arrays

FUNCTIONS = ("gauss", "cb", "dscb", "exponential", "power", "chebychev", "bernstein", "exp_poly", "endpoint", "uniform")


class PdfError(ValueError):
    pass


def _setting(value):
    """(init, lo, hi, fixed) from a number or an (init, lo, hi) tuple."""
    if isinstance(value, (int, float)):
        return float(value), float(value), float(value), True
    v = tuple(float(x) for x in value)
    if len(v) == 3:
        return v[0], v[1], v[2], False
    if len(v) == 2:
        return 0.5 * (v[0] + v[1]), v[0], v[1], False
    raise PdfError(f"parameter setting {value!r}: use a number (fixed) or (init, lo, hi)")


def default_params(function, xmin, xmax, peak, width, order):
    """Default (init, lo, hi) of each parameter of a named function."""
    span = xmax - xmin
    if function in ("dscb", "cb", "gauss"):
        d = {"mean": (peak, xmin - span, xmax + span), "sigma": (width, 1e-3 * span, 2. * span)}
        if function == "dscb":
            d.update({"sigmaR": (width, 1e-3 * span, 2. * span),
                      "alphaL": (1.5, 0.1, 10.), "nL": (3., 1.01, 50.),
                      "alphaR": (1.5, 0.1, 10.), "nR": (3., 1.01, 50.)})
        elif function == "cb":
            d.update({"alpha": (1.5, 0.1, 10.), "n": (3., 1.01, 50.)})
        return d
    if function == "exponential":
        return {"slope": (-1. / span, -100. / span, 100. / span)}
    if function == "power":
        return {"index": (-2., -50., 50.)}
    if function == "chebychev":
        return {f"c{i}": (0., -1., 1.) for i in range(1, order + 1)}
    if function == "bernstein":
        return {f"c{i}": (1., 0., 100.) for i in range(order + 1)}
    if function == "exp_poly":
        return {f"c{i}": ((-1.5 if i == 1 else 0.), -20., 20.) for i in range(1, order + 1)}
    if function == "endpoint":
        d = {"e0": (xmax + 0.45 * span, xmax + 0.001 * span, xmax + 3. * span), "a": (3., 0., 20.)}
        d.update({f"c{i}": (0., -2., 2.) for i in range(1, order + 1)})
        return d
    if function == "uniform":
        return {}
    raise PdfError(f"unknown function '{function}' (available: {', '.join(FUNCTIONS)}, or pass a RooAbsPdf)")


_quiet = False


def quiet_roofit():
    """Show only RooFit warnings and errors."""
    global _quiet
    if not _quiet:
        import ROOT
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
        _quiet = True


def hist_moments(h):
    """(peak, width) of a histogram, used to seed the default parameters."""
    e, values, _ = arrays(h)
    centers = 0.5 * (e[1:] + e[:-1])
    w = np.clip(values, 0., None)
    if w.sum() <= 0:
        return float(centers.mean()), 0.1 * (e[-1] - e[0])
    peak = float(centers[np.argmax(w)])
    mean = float(np.average(centers, weights=w))
    width = float(np.sqrt(np.average((centers - mean) ** 2, weights=w)))
    return peak, max(width, 0.5 * float(np.diff(e).min()))


def make_pdf(function, obs, name, hist, order=None, params=None):
    """A predefined function as a RooFit PDF named `name` on obs, its parameters named <name>_<parameter> and started
    from the histogram's peak and width. Returns (pdf, objects to keep alive)."""
    import ROOT
    xmin, xmax = obs.getMin(), obs.getMax()
    peak, width = hist_moments(hist)
    order = (0 if function == "endpoint" else 2) if order is None else int(order)
    settings = {k: (*v, False) for k, v in default_params(function, xmin, xmax, peak, width, order).items()}
    if function == "bernstein":
        settings["c0"] = (1., 1., 1., True)  # the overall scale is not a shape parameter
    for pname, value in (params or {}).items():
        if pname not in settings:
            raise PdfError(f"{name}: '{function}' has no parameter '{pname}' (has {', '.join(settings) or 'none'})")
        settings[pname] = _setting(value)
    keep = []
    v = {}
    for pname, (init, lo, hi, fixed) in settings.items():
        var = ROOT.RooRealVar(f"{name}_{pname}", pname, init, min(lo, init), max(hi, init))
        var.setConstant(fixed)
        v[pname] = var
        keep.append(var)
    center, half = 0.5 * (xmin + xmax), 0.5 * (xmax - xmin)
    t = f"((@0 - {center:.10g})/{half:.10g})"
    if function == "dscb":
        pdf = ROOT.RooCrystalBall(name, name, obs, v["mean"], v["sigma"], v["sigmaR"], v["alphaL"], v["nL"],
                                  v["alphaR"], v["nR"])
    elif function == "cb":
        pdf = ROOT.RooCrystalBall(name, name, obs, v["mean"], v["sigma"], v["alpha"], v["n"])
    elif function == "gauss":
        far, n = ROOT.RooConstVar(f"{name}_far", "", 50.), ROOT.RooConstVar(f"{name}_n", "", 2.)
        keep += [far, n]
        pdf = ROOT.RooCrystalBall(name, name, obs, v["mean"], v["sigma"], far, n, far, n)
    elif function == "exponential":
        pdf = ROOT.RooExponential(name, name, obs, v["slope"])
    elif function == "power":
        if xmin <= 0:
            raise PdfError(f"{name}: 'power' needs a positive observable range")
        pdf = ROOT.RooPowerLaw(name, name, obs, v["index"])
    elif function in ("chebychev", "bernstein"):
        coeffs = ROOT.RooArgList()
        for var in v.values():
            coeffs.add(var)
        keep.append(coeffs)
        pdf = (ROOT.RooChebychev if function == "chebychev" else ROOT.RooBernstein)(name, name, obs, coeffs)
    elif function == "exp_poly":
        args = ROOT.RooArgList(obs)
        terms = []
        for i in range(1, order + 1):
            args.add(v[f"c{i}"])
            terms.append(f"@{i}*TMath::Power({t}, {i})")
        pdf = ROOT.RooGenericPdf(name, name, f"TMath::Exp({' + '.join(terms) or '0'})", args)
    elif function == "endpoint":
        args = ROOT.RooArgList(obs, v["e0"], v["a"])
        formula = "TMath::Power(TMath::Max(1.e-9, @1 - @0), @2)"
        if order > 0:
            formula += "*TMath::Max(1.e-9, 1" + "".join(f" + @{i + 2}*TMath::Power({t}, {i})"
                                                         for i in range(1, order + 1)) + ")"
            for i in range(1, order + 1):
                args.add(v[f"c{i}"])
        pdf = ROOT.RooGenericPdf(name, name, formula, args)
    elif function == "uniform":
        pdf = ROOT.RooUniform(name, name, ROOT.RooArgSet(obs))
    else:
        raise PdfError(f"unknown function '{function}'")
    keep.append(pdf)
    return pdf, keep


def floating_params(pdf, obs, exclude=()):
    """The non-constant RooRealVar parameters of pdf (other than obs), except those named in exclude."""
    import ROOT
    return [p for p in pdf.getParameters(ROOT.RooArgSet(obs))
            if isinstance(p, ROOT.RooRealVar) and not p.isConstant() and p.GetName() not in exclude]


def fit_pdf(pdf, obs, hist, params, fit_range=None, hold=()):
    """Weighted ML fit of the RooRealVars `params` to a histogram (SumW2 errors). The variables in `hold` (e.g.
    nuisances inside the PDF) are kept constant during the fit. Returns the fit status."""
    import ROOT
    if not params:
        return 0
    saved = [(v, v.isConstant()) for v in hold]
    for v, _ in saved:
        v.setConstant(True)
    dh = ROOT.RooDataHist(f"{pdf.GetName()}_fitdata", "", ROOT.RooArgList(obs), hist)
    opts = [ROOT.RooFit.Save(True), ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.Warnings(False), ROOT.RooFit.PrintEvalErrors(-1)]
    if fit_range is not None:
        obs.setRange("m2c_fit", float(fit_range[0]), float(fit_range[1]))
        opts.append(ROOT.RooFit.Range("m2c_fit"))
    result = pdf.fitTo(dh, *opts)
    for v, const in saved:
        v.setConstant(const)
    return result.status() if result else -1


def add_nuisances(pdf, name, shifts, thetas):
    """A copy of pdf (RooCustomizer) in which each parameter p becomes
    p + sum_k theta_k * (theta_k > 0 ? up_k : down_k).

    shifts: {RooRealVar p: {nuisance: (up, down)}}; thetas: {nuisance: RooRealVar}, created as needed (shared by the
    PDFs of one channel). Returns (new pdf named `name`, objects to keep alive), or (pdf, []) without shifts."""
    import ROOT
    shifts = {p: s for p, s in shifts.items() if any(abs(u) > 0 or abs(d) > 0 for u, d in s.values())}
    if not shifts:
        return pdf, []
    cust = ROOT.RooCustomizer(pdf, "m2c")
    keep = [cust]
    for p, by_nuisance in shifts.items():
        args = ROOT.RooArgList(p)
        terms = []
        for nuis, (up, down) in by_nuisance.items():
            theta = thetas.get(nuis)
            if theta is None:
                theta = thetas[nuis] = ROOT.RooRealVar(nuis, nuis, 0., -7., 7.)
            i = args.getSize()
            args.add(theta)
            terms.append(f"@{i}*((@{i}>0)*({up:.10g}) + (@{i}<=0)*({down:.10g}))")
        f = ROOT.RooFormulaVar(f"{p.GetName()}_eff", "", "@0 + " + " + ".join(terms), args)
        cust.replaceArg(p, f)
        keep += [args, f]
    new = cust.build(False)
    new.SetName(name)
    keep.append(new)
    return new, keep


def bin_fractions(pdf, obs, edges):
    """Fraction of a PDF (normalized over the observable range) in each bin, from RooFit integrals over the bins."""
    import ROOT
    nset = ROOT.RooArgSet(obs)
    out = np.zeros(len(edges) - 1)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        obs.setRange("m2c_bin", float(lo), float(hi))
        integral = pdf.createIntegral(nset, ROOT.RooFit.NormSet(nset), ROOT.RooFit.Range("m2c_bin"))
        ROOT.SetOwnership(integral, True)
        out[i] = integral.getVal()
    total = out.sum()
    return out / total if total > 0 else out
