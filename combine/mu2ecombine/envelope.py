"""Data-driven background envelope: function families fit to the data sidebands, packed into a RooMultiPdf that
Combine profiles over (discrete profiling).

    Envelope(families={"chebychev": (1, 4), "exp_poly": (1, 3)}, absorbs=["rmc*", "rpc*"])

Options:
  families:      {family: (min_order, max_order)}; default DEFAULT_FAMILIES (ALL_FAMILIES lists every family)
  chisq_p_min:   accept an order if p(chi2, ndof) exceeds this
  ftest_p_max:   stop a family once the next order does not improve with p(F) < this
  enforce_ftest: apply the F-test
  force_orders:  {family: order} keeps only these orders (no chi2/F-test selection)
  keep_all:      keep every tested order
  penalty:       RooMultiPdf NLL correction per floating parameter (Combine's default 0.5)
  starts:        fit starting points per function: the default one plus random ones

The channel's blind window is excluded from the fits. The other background processes of the channel are held fixed at
their expected shapes and rates during the fits and stay in the model with their own uncertainties. Parameters are set
relative to the observable range, so the families do not depend on the range.
"""

import math

import numpy as np
import ROOT
from scipy import stats  # chi2 / F-test p-values only

from .hist import arrays, make
from .pdfs import bin_fractions, quiet_roofit


class EnvelopeError(ValueError):
    pass

# family -> (orders scanned, label)
DEFAULT_FAMILIES = {
    "chebychev": (1, 4),
    "exponential": (1, 3),
    "power": (1, 3),
    "exp_poly": (1, 4),
    "endpoint": (0, 2),
    "gauss_poly": (0, 3),
    "gaussian": (1, 1),
}
ALL_FAMILIES = ("chebychev", "bernstein", "exponential", "power", "laurent", "inv_poly", "exp_poly", "endpoint",
                "gaussian", "gauss_poly", "gauss_expo", "gauss_power")

LABELS = {
    "chebychev": "Chebychev", "bernstein": "Bernstein", "exponential": "Exponential sum", "power": "Power-law sum",
    "laurent": "Laurent", "inv_poly": "Inverse polynomial", "exp_poly": "exp(polynomial)",
    "endpoint": "Endpoint power law", "gaussian": "Broad Gaussian", "gauss_poly": "Gaussian + polynomial",
    "gauss_expo": "Gaussian + exponential", "gauss_power": "Gaussian + power law",
}


class _Keeper(list):
    """Holds Python references to RooFit objects that must outlive their constructors' scope."""

    def var(self, name, init, lo, hi):
        v = ROOT.RooRealVar(name, name, init, lo, hi)
        self.append(v)
        return v

    def add(self, obj):
        self.append(obj)
        return obj


#----------------------------------------------------------------------------------------------------------------------
# Families. Each returns a RooAbsPdf (or None if the order is not available).

def make_family_pdf(family, obs, order, name, keep):
    xmin, xmax = obs.getMin(), obs.getMax()
    span = xmax - xmin
    center, half = 0.5 * (xmin + xmax), 0.5 * span
    t = f"((@0 - {center:.8g})/{half:.8g})"  # observable scaled to [-1, 1]
    base = f"{name}_{family}{order}"

    if family == "chebychev":
        if order < 1:
            return None
        coeffs = ROOT.RooArgList()
        for i in range(1, order + 1):
            coeffs.add(keep.var(f"{base}_c{i}", (-0.7 if i == 1 else 0.1 / 10 ** (i - 1)), -1., 1.))
        return keep.add(ROOT.RooChebychev(base, f"Chebychev, order {order}", obs, coeffs))

    if family == "bernstein":
        if order < 1:
            return None
        coeffs = ROOT.RooArgList()
        for i in range(order + 1):
            coeffs.add(keep.var(f"{base}_c{i}", 1. / (1 + i), 0., 10.))
        return keep.add(ROOT.RooBernstein(base, f"Bernstein, order {order}", obs, coeffs))

    if family == "exponential":
        if order < 1:
            return None
        pdfs, fracs = ROOT.RooArgList(), ROOT.RooArgList()
        for i in range(1, order + 1):
            c = keep.var(f"{base}_c{i}", -1. / span * i, -50. / span, 50. / span)
            pdfs.add(keep.add(ROOT.RooExponential(f"{base}_e{i}", "", obs, c)))
            if i < order:
                fracs.add(keep.var(f"{base}_f{i}", 0.5, 0., 1.))
        if order == 1:
            pdf = pdfs.at(0)
            pdf.SetName(base)
            pdf.SetTitle("Exponential, order 1")
            return pdf
        return keep.add(ROOT.RooAddPdf(base, f"Exponential sum, order {order}", pdfs, fracs, True))

    if family == "power":
        if order < 1 or xmin <= 0:
            return None
        pdfs, fracs = ROOT.RooArgList(), ROOT.RooArgList()
        for i in range(1, order + 1):
            c = keep.var(f"{base}_p{i}", -2. * i, -100., 20.)
            pdfs.add(keep.add(ROOT.RooPowerLaw(f"{base}_pw{i}", "", obs, c)))
            if i < order:
                fracs.add(keep.var(f"{base}_f{i}", 0.5, 0., 1.))
        if order == 1:
            pdf = pdfs.at(0)
            pdf.SetName(base)
            pdf.SetTitle("Power law, order 1")
            return pdf
        return keep.add(ROOT.RooAddPdf(base, f"Power-law sum, order {order}", pdfs, fracs, True))

    if family == "laurent":
        if order < 1 or order > 6 or xmin <= 0:
            return None
        # sum of x^-k, k = 2..(order + 1), scaled to the range
        args = ROOT.RooArgList(obs)
        terms = [f"TMath::Power(@0/{xmin:.8g}, -2)"]
        for i in range(1, order):
            args.add(keep.var(f"{base}_c{i}", 0.5, 0., 1e3))
            terms.append(f"@{i}*TMath::Power(@0/{xmin:.8g}, {-(i + 2)})")
        return keep.add(ROOT.RooGenericPdf(base, f"Laurent, order {order}", " + ".join(terms), args))

    if family == "inv_poly":
        if order < 1 or order > 4:
            return None
        args = ROOT.RooArgList(obs)
        terms = ["1"]
        for i in range(1, order + 1):
            args.add(keep.var(f"{base}_c{i}", 0.1, -0.99, 10.))
            terms.append(f"@{i}*TMath::Power({t}, {i})")
        return keep.add(ROOT.RooGenericPdf(base, f"Inverse polynomial, order {order}",
                                           f"1./TMath::Max(1.e-10, {' + '.join(terms)})", args))

    if family == "exp_poly":
        if order < 1 or order > 5:
            return None
        args = ROOT.RooArgList(obs)
        terms = []
        for i in range(1, order + 1):
            args.add(keep.var(f"{base}_c{i}", -1.5 if i == 1 else 0., -20., 20.))
            terms.append(f"@{i}*TMath::Power({t}, {i})")
        return keep.add(ROOT.RooGenericPdf(base, f"exp(polynomial), order {order}",
                                           f"TMath::Exp({' + '.join(terms)})", args))

    if family == "endpoint":
        if order < 0 or order > 3:
            return None
        args = ROOT.RooArgList(obs)
        args.add(keep.var(f"{base}_e0", xmax + 0.45 * span, xmax + 0.01 * span, xmax + 3. * span))
        args.add(keep.var(f"{base}_a", 3., 0., 20.))
        formula = "TMath::Power(TMath::Max(1.e-9, @1 - @0), @2)"
        if order > 0:
            formula += "*TMath::Max(1.e-9, 1"
            for i in range(1, order + 1):
                args.add(keep.var(f"{base}_c{i}", 0., -2., 2.))
                formula += f" + @{i + 2}*TMath::Power({t}, {i})"
            formula += ")"
        return keep.add(ROOT.RooGenericPdf(base, f"Endpoint power law, order {order}", formula, args))

    # Gaussian-based families: a broad Gaussian whose peak may lie below or inside the range
    def gauss_vars():
        mean = keep.var(f"{base}_mean", xmin - 1.2 * span, xmin - 6. * span, xmax)
        sigma = keep.var(f"{base}_sigma", 1.2 * span, 0.3 * span, 3. * span)
        return mean, sigma

    if family == "gaussian":
        if order != 1:
            return None
        mean, sigma = gauss_vars()
        # a Crystal Ball with far, fixed tails: Combine treats a top-level RooGaussian as a constraint term
        far = keep.add(ROOT.RooConstVar(f"{base}_far", "", 50.))
        n = keep.add(ROOT.RooConstVar(f"{base}_n", "", 2.))
        return keep.add(ROOT.RooCrystalBall(base, "Broad Gaussian", obs, mean, sigma, far, n, far, n))

    if family in ("gauss_poly", "gauss_expo", "gauss_power"):
        lo_order = 0 if family == "gauss_poly" else 1
        if order < lo_order or order > (3 if family == "gauss_poly" else 2):
            return None
        mean, sigma = gauss_vars()
        args = ROOT.RooArgList(obs, mean, sigma)
        formula = "TMath::Gaus(@0, @1, @2)"
        k = 3
        if family == "gauss_poly":
            for i in range(order + 1):
                args.add(keep.var(f"{base}_p{i}", 0., -1., 1.))
                formula += f" + @{k}*TMath::Power({t}, {i})"
                k += 1
            formula = f"TMath::Max(0., {formula})"
        elif family == "gauss_expo":
            for i in range(order):
                args.add(keep.var(f"{base}_n{i}", -1., -20., 5.))
                args.add(keep.var(f"{base}_s{i}", -1., -20., 20.))
                formula += f" + TMath::Exp(@{k} + @{k + 1}*{t})"
                k += 2
        else:
            if xmin <= 0:
                return None
            for i in range(order):
                args.add(keep.var(f"{base}_n{i}", 0.1, 0., 10.))
                args.add(keep.var(f"{base}_s{i}", -2., -50., 10.))
                formula += f" + @{k}*TMath::Power(@0/{xmin:.8g}, @{k + 1})"
                k += 2
        return keep.add(ROOT.RooGenericPdf(base, f"{LABELS[family]}, order {order}", formula, args))

    raise EnvelopeError(f"unknown envelope family '{family}' (available: {', '.join(ALL_FAMILIES)})")


def count_floating(pdf, obs):
    params = pdf.getParameters(ROOT.RooArgSet(obs))
    return sum(1 for p in params if isinstance(p, ROOT.RooRealVar) and not p.isConstant()
               and not p.GetName().endswith("_norm"))


#----------------------------------------------------------------------------------------------------------------------

class EnvelopeFit:
    def __init__(self, family, order, pdf, chi2, ndof, nparams, accepted, reason, pvalue, ftest):
        self.family, self.order, self.pdf = family, order, pdf
        self.chi2, self.ndof, self.nparams = chi2, ndof, nparams
        self.accepted, self.reason, self.pvalue, self.ftest = accepted, reason, pvalue, ftest
        self.expected = None  # expected events per bin (fit model incl. fixed processes)

    @property
    def label(self):
        return f"{LABELS[self.family]}, order {self.order}"

    def as_dict(self):
        return {"family": self.family, "order": self.order, "pdf": self.pdf.GetName(), "chi2": self.chi2,
                "ndof": self.ndof, "nparams": self.nparams, "p_chi2": self.pvalue, "p_ftest": self.ftest,
                "accepted": self.accepted, "reason": self.reason}


class Envelope:
    """A background shape from discrete profiling over function families (see the module docstring for options).

    build() arguments: obs, the channel observable; data, the TH1 (channel binning) fit outside the blind window; fixed,
    the expected events per bin of the explicitly modeled backgrounds (numpy array), or None.
    """

    def __init__(self, families=None, absorbs=(), chisq_p_min=0.001, ftest_p_max=0.2, enforce_ftest=True,
                 force_orders=None, keep_all=False, penalty=0.5, starts=4):
        fams = {k: tuple(v) for k, v in (families or DEFAULT_FAMILIES).items()}
        for fam, orders in fams.items():
            if fam not in ALL_FAMILIES:
                raise EnvelopeError(f"unknown envelope family '{fam}' (available: {', '.join(ALL_FAMILIES)})")
            if len(orders) != 2:
                raise EnvelopeError(f"envelope family {fam}: expected (min_order, max_order)")
        self.families = fams
        self.absorbs = [absorbs] if isinstance(absorbs, str) else list(absorbs)
        self.chisq_p_min = float(chisq_p_min)
        self.ftest_p_max = float(ftest_p_max)
        self.enforce_ftest = bool(enforce_ftest)
        self.force_orders = dict(force_orders or {})
        self.keep_all = bool(keep_all)
        self.penalty = float(penalty)
        self.starts = max(1, int(starts))
        self.where = "envelope"
        self.keep = _Keeper()
        self.fits = []

    def build(self, name, obs, data, fixed, blind, log=print):
        quiet_roofit()
        level = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = max(level, ROOT.kError)  # Hessian warnings from the random starting points
        try:
            return self._build(name, obs, data, fixed, blind, log)
        finally:
            ROOT.gErrorIgnoreLevel = level

    def _build(self, name, obs, data, fixed, blind, log):
        edges, data_vals, _ = arrays(data)
        centers = 0.5 * (edges[1:] + edges[:-1])
        if blind is None:
            raise EnvelopeError(f"{self.where}: the envelope needs the channel's blind window (it is fit to the sidebands)")
        side = ~((centers > blind[0]) & (centers < blind[1]))
        obs.setRange("m2c_low", float(edges[0]), float(blind[0]))
        obs.setRange("m2c_high", float(blind[1]), float(edges[-1]))
        n_data_side = float(data_vals[side].sum())
        if n_data_side <= 0:
            raise EnvelopeError(f"{self.where}: no data in the sidebands")
        fixed_vals = np.asarray(fixed) if fixed is not None else np.zeros(len(data_vals))
        n_fixed_side = float(fixed_vals[side].sum())
        if n_fixed_side >= n_data_side:
            raise EnvelopeError(f"{self.where}: the fixed backgrounds ({n_fixed_side:.3g}) exceed the sideband data "
                              f"({n_data_side:.3g})")

        accepted = []
        for family, (lo, hi) in self.families.items():
            prev = None
            forced = self.force_orders.get(family)
            for order in range(int(lo), int(hi) + 1):
                pdf = make_family_pdf(family, obs, order, name, self.keep)
                if pdf is None:
                    continue
                # nested families: also start from the previous order's best fit (new terms at their defaults)
                if prev is not None:
                    prev_vals = _suffix_values(prev.pdf, obs, f"{name}_{family}{prev.order}_")
                    seed = {}
                    for p_ in pdf.getParameters(ROOT.RooArgSet(obs)):
                        suffix = p_.GetName()[len(f"{name}_{family}{order}_"):]
                        if suffix in prev_vals:
                            seed[p_.GetName()] = prev_vals[suffix]
                else:
                    seed = None
                expected, chi2, nfit = self._fit(pdf, obs, edges, data_vals, fixed_vals, side, n_data_side,
                                                 n_fixed_side, seed)
                nparams = count_floating(pdf, obs)
                ndof = nfit - nparams - 1
                pvalue = float(stats.chi2.sf(chi2, ndof)) if ndof > 0 else 0.
                ftest = None
                if self.keep_all:
                    ok, reason = True, "keep_all"
                elif forced is not None:
                    ok, reason = order == forced, "forced order" if order == forced else f"not the forced order {forced}"
                else:
                    ok = pvalue > self.chisq_p_min
                    reason = f"p(chi2) = {pvalue:.3g} {'>' if ok else '<'} {self.chisq_p_min:g}"
                stop = False
                if prev is not None and self.enforce_ftest and not self.keep_all and forced is None:
                    dchi2 = prev.chi2 - chi2
                    dndof = prev.ndof - ndof
                    if dndof > 0:
                        ftest = float(stats.chi2.sf(max(dchi2, 0.), dndof)) if dchi2 > 0 else 1.
                    elif dndof == 0:
                        ftest = 0. if chi2 < prev.chi2 else 1.
                    if ftest is not None and ftest > self.ftest_p_max:
                        ok, stop = False, True
                        reason = f"no significant improvement over order {prev.order} (p(F) = {ftest:.3g})"
                fit = EnvelopeFit(family, order, pdf, chi2, ndof, nparams, ok, reason, pvalue, ftest)
                fit.expected = expected
                self.fits.append(fit)
                log(f"    {fit.label:34s} {nparams:2d} par  chi2/ndof = {chi2:8.2f}/{ndof:3d}  p = {pvalue:8.2e}"
                    f"{'' if ftest is None else f'  p(F) = {ftest:8.2e}'}  "
                    f"{'ACCEPTED' if ok else 'rejected'} ({reason})")
                if ok:
                    accepted.append(fit)
                prev = fit
                if stop:
                    break
        # A Gaussian is exp(polynomial) of order 2: keep only one of the two
        if any(f.family == "exp_poly" and f.order == 2 for f in accepted):
            dropped = [f for f in accepted if f.family == "gaussian"]
            for f in dropped:
                f.accepted, f.reason = False, "same shape as exp(polynomial), order 2"
                log(f"    {f.label}: dropped, the same shape as exp(polynomial), order 2")
            accepted = [f for f in accepted if f.family != "gaussian"]
        if not accepted:
            raise EnvelopeError(f"{self.where}: no function was accepted into the envelope")
        # starting function: the best chi2 + penalty, as Combine would choose on the fitted data
        scores = [f.chi2 + 2. * self.penalty * f.nparams for f in accepted]
        best = int(np.argmin(scores))
        cat = self.keep.add(ROOT.RooCategory(f"{name}_cat", f"{name} function index"))
        pdfs = ROOT.RooArgList()
        for f in accepted:
            pdfs.add(f.pdf)
            _freeze(f.pdf, obs, False)
        multi = self.keep.add(ROOT.RooMultiPdf(name, "Background envelope", cat, pdfs))
        multi.setCorrectionFactor(self.penalty)
        cat.setIndex(best)
        # starting yield over the full range from the best function's sideband fraction
        frac_side = float(bin_fractions(accepted[best].pdf, obs, edges)[side].sum())
        n_env = (n_data_side - n_fixed_side) / max(frac_side, 1e-12)
        self.accepted = accepted
        self.best = best
        self.rate = n_env
        self.multipdf = multi
        self.cat = cat
        return multi, cat, n_env

    def _fit(self, pdf, obs, edges, data_vals, fixed_vals, side, n_data_side, n_fixed_side, seed=None):
        """RooFit extended fit of pdf + the fixed backgrounds to the sideband data (the fixed part held at its
        expected sideband yield, the function's sideband yield floating), from several starting points. Returns the
        expected events per bin (all bins, function normalized to the fitted sideband yield), the Pearson chi2 over
        the sidebands, and the number of sideband bins."""
        params = [p for p in pdf.getParameters(ROOT.RooArgSet(obs))
                  if isinstance(p, ROOT.RooRealVar) and not p.isConstant()]
        data = self._roo_hist("data", obs, edges, data_vals)
        n_env = ROOT.RooRealVar(f"{pdf.GetName()}_sideband_yield", "", n_data_side - n_fixed_side, 0.,
                                10. * n_data_side)
        coefs, pdfs = ROOT.RooArgList(n_env), ROOT.RooArgList(pdf)
        keep = [n_env, data]
        if n_fixed_side > 0:
            fixed_dh = self._roo_hist("fixed", obs, edges, fixed_vals)
            fixed_pdf = ROOT.RooHistPdf(f"{pdf.GetName()}_fixed_pdf", "", ROOT.RooArgSet(obs), fixed_dh)
            n_fix = ROOT.RooRealVar(f"{pdf.GetName()}_fixed_yield", "", n_fixed_side)
            n_fix.setConstant(True)
            coefs.add(n_fix)
            pdfs.add(fixed_pdf)
            keep += [fixed_dh, fixed_pdf, n_fix]
        total = ROOT.RooAddPdf(f"{pdf.GetName()}_sideband_model", "", pdfs, coefs)
        keep.append(total)

        x0 = np.array([p.getVal() for p in params])
        starts = [x0]
        if seed:
            starts.append(np.array([seed.get(p.GetName(), v) for p, v in zip(params, x0)]))
        rng = np.random.default_rng(len(params) * 7919 + len(edges))
        for _ in range(self.starts - 1):
            starts.append(np.array([rng.uniform(p.getMin(), p.getMax()) for p in params]))
        best_nll, best_vals = None, None
        for start in starts:
            for p, v in zip(params, start):
                p.setVal(float(v))
            n_env.setVal(n_data_side - n_fixed_side)
            result = total.fitTo(data, ROOT.RooFit.Range("m2c_low,m2c_high"), ROOT.RooFit.Extended(True),
                                 ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Warnings(False),
                                 ROOT.RooFit.PrintEvalErrors(-1))
            nll = result.minNll() if result and result.status() in (0, 1) else None
            if nll is not None and np.isfinite(nll) and (best_nll is None or nll < best_nll):
                best_nll, best_vals = nll, [p.getVal() for p in params]
        if best_vals is None:
            return None, float("inf"), int(side.sum())
        for p, v in zip(params, best_vals):
            p.setVal(v)
        frac = bin_fractions(pdf, obs, edges)
        fs = frac[side].sum()
        if not np.isfinite(fs) or fs <= 0 or np.any(frac < 0):
            return None, float("inf"), int(side.sum())  # not a valid density over the range
        mu = fixed_vals + (n_data_side - n_fixed_side) * frac / fs
        m = np.clip(mu[side], 1e-12, None)
        chi2 = float(np.sum((data_vals[side] - m) ** 2 / m))
        return mu, chi2, int(side.sum())

    def _roo_hist(self, tag, obs, edges, values):
        h = make(edges, values, np.clip(values, 0, None), f"m2c_env_{tag}")
        dh = ROOT.RooDataHist(h.GetName() + "_dh", "", ROOT.RooArgList(obs), h)
        self.keep += [h, dh]
        return dh


def _suffix_values(pdf, obs, prefix):
    return {p.GetName()[len(prefix):]: p.getVal() for p in pdf.getParameters(ROOT.RooArgSet(obs))
            if isinstance(p, ROOT.RooRealVar) and p.GetName().startswith(prefix)}


def _freeze(pdf, obs, constant):
    for p in pdf.getParameters(ROOT.RooArgSet(obs)):
        if isinstance(p, ROOT.RooRealVar):
            p.setConstant(constant)
