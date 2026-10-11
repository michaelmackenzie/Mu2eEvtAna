//
// PhysicsWeights: per-event physics weights computed from MC truth
//   - RMC photon spectrum: Plestid phase-space model, for samples generated with a flat photon spectrum
//   - Beam intensity: re-weighting the simulated N(POT) per microbunch to a target intensity
// Michael MacKenzie (2026)

#ifndef MU2EEVTANA_PHYSICSWEIGHTS_HH
#define MU2EEVTANA_PHYSICSWEIGHTS_HH

// standard includes
#include <algorithm>
#include <cmath>
#include <regex>
#include <string>

namespace Mu2eEvtAna {

  //------------------------------------------------------------------------------------
  // Plestid phase-space RMC spectrum for a given number of nucleon knock-outs: p(x) ~ x (1 - x)^(2 + 1.5 n), x = E/kmax,
  // normalized to unit area over [0, kmax]
  inline double PlestidSpectrum(const double energy, const double kmax, const int knockout) {
    if(energy <= 0. || energy >= kmax || kmax <= 0. || knockout < 0) return 0.;
    const double power = 2. + 1.5*knockout;
    const double norm  = (power + 1.)*(power + 2.)/kmax;
    const double x = energy/kmax;
    return norm*x*std::pow(1. - x, power);
  }

  // Integral of PlestidSpectrum from k_1 to k_2
  inline double PlestidIntegral(double k_1, double k_2, const double kmax, const int knockout) {
    if(kmax <= 0. || knockout < 0) return 0.;
    k_1 = std::max(0., std::min(kmax, k_1));
    k_2 = std::max(0., std::min(kmax, k_2));
    if(k_1 >= k_2) return 0.;
    const double power = 2. + 1.5*knockout;
    auto primitive = [&](double x) { return (x - 1.)*std::pow(1. - x, power)*(power*x + x + 1.); };
    return primitive(k_2/kmax) - primitive(k_1/kmax);
  }

  // RMC photon spectrum: a mix of 0 and 1 knock-out channels, each normalized above a reference energy and weighted by its
  // branching fraction above it. The weight is a probability density (per MeV) that integrates to 1 above the reference
  // energy, so a flat sample normalized as R(RMC, E > E_ref) * (emax - emin) / N(gen) gives the physical spectrum.
  struct PlestidRMCWeight {
    double kmax_0;     // 0 knock-out kinematic endpoint (MeV)
    double kmax_1;     // 1 knock-out kinematic endpoint (MeV)
    double br_0;       // R(0 knock-out | E > ref) / R(RMC | E > ref)
    double br_1;       // R(1 knock-out | E > ref) / R(RMC | E > ref)
    double ref_energy; // reference energy (MeV)

    double Weight(const double energy) const {
      const double frac_0 = PlestidIntegral(ref_energy, kmax_0, kmax_0, 0);
      const double frac_1 = PlestidIntegral(ref_energy, kmax_1, kmax_1, 1);
      double weight = 0.;
      if(frac_0 > 0.) weight += br_0*PlestidSpectrum(energy, kmax_0, 0)/frac_0;
      if(frac_1 > 0.) weight += br_1*PlestidSpectrum(energy, kmax_1, 1)/frac_1;
      return weight;
    }

    // Aluminum (stopping target)
    static PlestidRMCWeight Aluminum() { return {101.866, 95.449, 0.099, 0.901, 57.}; }
    // Carbon (polyethylene target): kinematic endpoints from the 1992 carbon fit, approximate branching fractions
    static PlestidRMCWeight Carbon  () { return { 91.30 , 87.96 , 0.20 , 0.80 , 57.}; }
  };

  //------------------------------------------------------------------------------------
  // Beam intensity distributions: log-normal in N(POT) per microbunch, with mean mu and SDF = exp(-sigma^2)
  inline double LogNormalPDF(const double x, const double mu, const double sdf) {
    if(x <= 0. || sdf <= 0. || mu <= 0.) return 0.;
    const double sigma = std::sqrt(-std::log(sdf));
    const double mu0   = std::log(mu) - 0.5*sigma*sigma;
    return 1./(x*sigma*std::sqrt(2.*M_PI))*std::exp(-std::pow(std::log(x) - mu0, 2)/(2.*sigma*sigma));
  }

  inline double LogNormalCDF(const double x, const double mu, const double sdf) {
    if(x <= 0. || sdf <= 0. || mu <= 0.) return 0.;
    const double sigma = std::sqrt(-std::log(sdf));
    const double mu0   = std::log(mu) - 0.5*sigma*sigma;
    return 0.5*std::erfc(-(std::log(x) - mu0)/(sigma*std::sqrt(2.)));
  }

  // Intensity-weighted (x * log-normal) distribution: the intensity seen by a process whose rate scales with N(POT)
  inline double XLogNormalPDF(const double x, const double mu, const double sdf) {
    if(x <= 0. || sdf <= 0. || mu <= 0.) return 0.;
    const double sigma = std::sqrt(-std::log(sdf));
    const double mu0   = std::log(mu) - 0.5*sigma*sigma;
    const double d     = std::log(x) - (mu0 + sigma*sigma);
    return 1./(x*sigma*std::sqrt(2.*M_PI))*std::exp(-0.5*d*d/(sigma*sigma));
  }

  inline double XLogNormalCDF(const double x, const double mu, const double sdf) {
    if(x <= 0. || sdf <= 0. || sdf > 1. || mu <= 0.) return 0.;
    const double sigma = std::sqrt(-std::log(sdf));
    const double mu0   = std::log(mu) - 0.5*sigma*sigma;
    return 0.5*std::erfc(-(std::log(x) - (mu0 + sigma*sigma))/(sigma*std::sqrt(2.)));
  }

  // Weight from the simulated to the target intensity distribution, each normalized over the simulated range [0, max].
  // Processes whose rate scales with N(POT) (beam-induced signals) follow x * log-normal; pileup-only and cosmic-ray
  // samples follow the log-normal.
  struct BeamIntensityWeight {
    double mu_nominal  = 0.  ; // simulated mean N(POT) per microbunch
    double sdf_nominal = 0.8 ;
    double mu_goal     = 0.  ; // target mean N(POT) per microbunch
    double sdf_goal    = 0.8 ;
    double max         = 0.  ; // simulated cut-off

    bool Active() const { return mu_nominal > 0. && mu_goal > 0. && max > 0. && (mu_nominal != mu_goal || sdf_nominal != sdf_goal); }

    double Weight(const double npot, const bool scales_with_pot) const {
      if(npot <= 1. || !Active()) return 1.;
      double p_1, p_2, c_1, c_2;
      if(scales_with_pot) {
        p_1 = XLogNormalPDF(npot, mu_nominal, sdf_nominal); c_1 = XLogNormalCDF(max, mu_nominal, sdf_nominal);
        p_2 = XLogNormalPDF(npot, mu_goal   , sdf_goal   ); c_2 = XLogNormalCDF(max, mu_goal   , sdf_goal   );
      } else {
        p_1 = LogNormalPDF (npot, mu_nominal, sdf_nominal); c_1 = LogNormalCDF (max, mu_nominal, sdf_nominal);
        p_2 = LogNormalPDF (npot, mu_goal   , sdf_goal   ); c_2 = LogNormalCDF (max, mu_goal   , sdf_goal   );
      }
      if(p_1 <= 0. || c_1 <= 0. || c_2 <= 0.) return 0.;
      return (p_2/c_2)/(p_1/c_1);
    }
  };

  //------------------------------------------------------------------------------------
  // Dataset tag: <process (5 characters)>b<booster batches>s<stream>r<version>, e.g. cele0b1s51r0204, found anywhere in a name
  struct DatasetTag_t {
    std::string process ; // e.g. cele0
    int         batches = -1; // 0: no pileup, 1: 1BB, 2: 2BB
    int         version = -1; // e.g. 204 for r0204
  };

  inline bool ParseDatasetTag(const std::string& name, DatasetTag_t& tag) {
    static const std::regex pattern("([a-z0-9]{5})b([0-9])s[0-9]+r([0-9]{4})");
    std::smatch match;
    if(!std::regex_search(name, match, pattern)) return false;
    tag.process = match[1].str();
    tag.batches = std::stoi(match[2].str());
    tag.version = std::stoi(match[3].str());
    return true;
  }

  // Run 1B pileup samples (b1, r02XX), simulated at a reduced 1BB intensity
  inline bool IsRun1BSample(const DatasetTag_t& tag) { return tag.batches == 1 && tag.version / 100 == 2; }

  // Simulated beam intensity (production mixing PBISim.extendedMean, .cutMax) of a sample: false if unknown or without pileup
  //   Run 1B (b1r02XX): 5.92e6, 1BB (b1r00XX/b1r01XX): 1.58e7, 2BB (b2): 3.93e7, each cut off at 6x the mean
  inline bool NominalBeamIntensity(const DatasetTag_t& tag, double& mu, double& max) {
    mu = max = 0.;
    if     (IsRun1BSample(tag))                          { mu = 5.92e6; max = 35.5e6; }
    else if(tag.batches == 1 && tag.version / 100 <= 1) { mu = 1.58e7; max = 9.48e7; }
    else if(tag.batches == 2)                            { mu = 3.93e7; max = 2.36e8; }
    return mu > 0.;
  }

  // Whether a sample's rate scales with N(POT) (x * log-normal intensity): all but cosmic-ray (cry*, csms) and
  // no-primary pileup (mnbs) samples
  inline bool BeamRateScalesWithPOT(const DatasetTag_t& tag) {
    const std::string p = tag.process.substr(0, 4);
    return !(p.rfind("cry", 0) == 0 || p == "csms" || p == "mnbs");
  }
}

#endif
