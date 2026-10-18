"""mu2ecombine: build Combine models for Mu2e analyses in Python, from ROOT histograms and RooFit PDFs.

    import mu2ecombine as mc
    p = mc.observable("p", 100, 110, bin_width=0.5)       # a RooRealVar with its binning
    model = mc.Model("search", poi_range=(0, 50))
    sr = model.channel("sr", p)
    sr.add("signal", mc.Fit("dscb", signal_th1), signal=True)
    sr.add("dio", dio_th1)
    model.lnN("lumi", 1.1)
    model.write("output/search")

Combine does all of the statistics; the output directory is a plain Combine model (see README.md).
"""

from .hist import convolved, rebin_to, replace_tail, total, transformed
from .inputs import InputError, Mu2eEvtAnaFile, Sample, mu2eevtana_systematics
from .physics import mu2e_physics
from .shapes import Fit, Pdf, Template
from .envelope import Envelope, ALL_FAMILIES, DEFAULT_FAMILIES
from .model import Model, ModelError, observable
from .narrow import narrow_model

__version__ = "0.3.0"

__all__ = ["convolved", "rebin_to", "replace_tail", "total", "transformed", "InputError", "Mu2eEvtAnaFile", "Sample",
           "mu2eevtana_systematics", "mu2e_physics", "Fit", "Pdf", "Template", "Envelope", "ALL_FAMILIES",
           "DEFAULT_FAMILIES", "Model", "ModelError", "observable", "narrow_model"]
