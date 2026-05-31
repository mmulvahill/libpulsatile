#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <cmath>
#include <bp_mcmc/counter.h>
#include <bp_mcmc/proposalvariance.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>
#include <bp_datastructures/pulseestimates.h>
#include <bp_datastructures/utils.h>
#include <bpmod_singlesubject/ss_draw_tvarscale.h>
#include <testing/catch.h>


//
// SS_DrawTVarScale::posterior_function
//
// Regression test for the t-variance-scale (kappa) acceptance ratio. The
// scaled-normal representation has theta ~ N(mu, sigma^2 / kappa), so the
// kappa MH log-ratio includes a quadratic data term (theta - mu)^2 * 0.5 * kappa.
// A prior bug multiplied this term by an uninitialized 0, silently dropping it;
// the kappa samplers then ignored how far the random effect sat from its mean.
//
TEST_CASE( "tvarscale posterior_function includes the quadratic data term",
           "[draw_][tvarscale]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();

  // Known inputs. A small sigma makes the quadratic data term sizeable, so the
  // dropped-term bug is unmistakable.
  const double mu     = 3.5;   // patient mass_mean
  const double sigma  = 1.0;   // patient mass_sd  (sigma^2 = 1)
  const double theta  = 8.0;   // pulse mass (random effect), far from mu
  const double kappa  = 1.0;   // current tvarscale_mass
  const double prop   = 2.0;   // proposed tvarscale

  pat.estimates.mass_mean = mu;
  pat.estimates.mass_sd   = sigma;

  PulseEstimates pulse;          // empty pulse; posterior_function only reads
  pulse.mass           = theta;  // mass and tvarscale_mass (no mean_contrib)
  pulse.tvarscale_mass = kappa;

  // for_width = false -> sample the mass tvarscale
  SS_DrawTVarScale sampler(1.0, 500, 25000, 0.35, false, false, 5000);

  double got = sampler.posterior_function(&pulse, prop, &pat);

  // Analytic log posterior ratio (same closed form the sampler should compute).
  const double dev2        = (theta - mu) * (theta - mu);
  const double prior_ratio = std::log(Rf_dgamma(prop, 2, 0.5, 0)) -
                             std::log(Rf_dgamma(kappa, 2, 0.5, 0));
  const double stdold      = theta / (sigma / std::sqrt(kappa));
  const double stdnew      = theta / (sigma / std::sqrt(prop));
  const double quad        = (dev2 * 0.5 * kappa - dev2 * 0.5 * prop) / (sigma * sigma);
  const double re_ratio    = quad +
                             Rf_pnorm5(stdold, 0, 1, 1, 1) -
                             Rf_pnorm5(stdnew, 0, 1, 1, 1) -
                             0.5 * std::log(kappa) + 0.5 * std::log(prop);
  const double expected    = prior_ratio + re_ratio;

  SECTION( "matches the analytic scaled-normal log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "regression guard: the quadratic term is present and non-trivial" ) {
    // The quadratic term here is large (~ -10.1). Under the original bug it was
    // identically 0, so `got` would equal `expected - quad`. Assert it does not.
    REQUIRE( quad != Approx(0.0) );
    REQUIRE( got  != Approx(expected - quad) );
  }
}
