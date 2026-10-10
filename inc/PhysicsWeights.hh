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
    double mu_nominal  = 5.92e6; // simulated mean N(POT) per microbunch (Run1Ban mixing)
    double sdf_nominal = 0.8   ;
    double mu_goal     = 6.14e6; // target mean N(POT) per microbunch (1.5 kW)
    double sdf_goal    = 0.8   ;
    double max         = 35.5e6; // simulated cut-off

    bool Active() const { return mu_nominal != mu_goal || sdf_nominal != sdf_goal; }

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
}

#endif
