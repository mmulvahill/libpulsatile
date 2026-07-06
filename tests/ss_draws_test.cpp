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
#include <bpmod_singlesubject/ss_draw_randomeffects.h>
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


//
// SS_DrawRandomEffects::posterior_function
//
// Regression test for the pulse random-effect (mass/width) acceptance ratio.
// The random effect is a scale-mixture-of-normals: theta_k ~ N(mu, sigma^2/kappa_k),
// so the prior quadratic term must be weighted by the per-pulse t-scale kappa_k.
// A prior bug omitted kappa entirely, leaving prior variance sigma^2 and making
// the random-effect draw inconsistent with the tvarscale and SD-of-random-effects
// full conditionals (which do use kappa). The inconsistency hit the weakly
// identified pulse *width* hardest, letting its pulse-to-pulse SD space-fill.
//
// The likelihood term does not depend on kappa, so evaluating posterior_function
// at two different kappa values (holding the random effect, mean, sd, and proposal
// fixed) isolates the kappa-weighted prior term: the likelihood contribution
// cancels in the difference. Under the old bug this difference was identically 0.
//
TEST_CASE( "random-effects posterior_function weights the prior by the t-scale kappa",
           "[draw_][randomeffects]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();

  // The test patient is constructed with a single initial pulse; use it as the
  // random effect being sampled (it must live inside pat.pulses so the internal
  // likelihood() reflects the proposed value).
  REQUIRE( pat.get_pulsecount() >= 1 );
  PulseEstimates & pulse = pat.pulses.front();

  // Known inputs. Small sigma + a random effect far from its mean make the
  // kappa-weighted prior term sizeable, so the dropped-kappa bug is unmistakable.
  const double mu     = 3.5;   // patient mass_mean
  const double sigma  = 1.0;   // patient mass_sd (sigma^2 = 1)
  const double theta  = 8.0;   // pulse mass (random effect), far from mu
  const double prop   = 2.0;   // proposed pulse mass
  const double kappa1 = 1.0;   // first per-pulse t-scale
  const double kappa2 = 2.0;   // second per-pulse t-scale

  pat.estimates.mass_mean = mu;
  pat.estimates.mass_sd   = sigma;
  pulse.mass              = theta;

  // for_width = false -> sample the pulse mass random effect
  SS_DrawRandomEffects sampler(1.0, 500, 25000, 0.35, false, false, 5000);

  pulse.tvarscale_mass = kappa1;
  double r1 = sampler.posterior_function(&pulse, prop, &pat);

  pulse.tvarscale_mass = kappa2;
  double r2 = sampler.posterior_function(&pulse, prop, &pat);

  // Only the kappa-weighted prior term differs between the two calls:
  //   r(kappa) = kappa * [0.5(theta-mu)^2 - 0.5(prop-mu)^2] / sigma^2 + (loglik term)
  const double prior_unit = (0.5 * (theta - mu) * (theta - mu) -
                             0.5 * (prop  - mu) * (prop  - mu)) / (sigma * sigma);
  const double expected_diff = (kappa2 - kappa1) * prior_unit;

  SECTION( "kappa scales the prior term as predicted by the scaled-normal model" ) {
    REQUIRE( r2 - r1 == Approx(expected_diff) );
  }

  SECTION( "regression guard: the acceptance ratio actually depends on kappa" ) {
    // prior_unit (= 9.0 here) is non-trivial, so a kappa-dependent ratio must
    // differ between kappa1 and kappa2. Under the original bug r2 - r1 == 0.
    REQUIRE( prior_unit != Approx(0.0) );
    REQUIRE( r2 != Approx(r1) );
  }
}


//
// SS_DrawRandomEffects::posterior_function -- log-normal parameterization
//
// Under the papers' log-normal model, log(theta_k) ~ N(mu, sigma^2/kappa_k),
// where mu/sigma are the mean/SD of the LOG random effect. The pulse stores the
// natural-scale value, and the symmetric random-walk proposal is on that natural
// value, so the target of theta carries the log-normal's 1/theta Jacobian; the
// MH log-ratio therefore adds (log(theta) - log(theta')) to the kappa-weighted
// quadratic term. This test checks the full closed form exactly (the likelihood
// delta is computed the same way the sampler does, so what is verified is the
// prior + Jacobian ratio).
//
TEST_CASE( "random-effects posterior_function: log-normal prior + Jacobian",
           "[draw_][randomeffects][lognormal]" ) {

  DataStructuresUtils utils;
  Patient pat = utils.create_new_test_patient_obj();
  pat.lognormal_pulses = true;   // papers' parameterization

  REQUIRE( pat.get_pulsecount() >= 1 );
  PulseEstimates & pulse = pat.pulses.front();

  const double mu     = 3.7;   // mean of LOG width (e.g. log(~40))
  const double sigma  = 0.6;   // SD of LOG width
  const double theta  = 55.0;  // current width (natural scale, > 0)
  const double prop   = 30.0;  // proposed width (natural scale, > 0)
  const double kappa  = 1.3;   // per-pulse t-scale

  pat.estimates.width_mean = mu;
  pat.estimates.width_sd   = sigma;
  pulse.width              = theta;
  pulse.tvarscale_width    = kappa;

  // for_width = true -> sample the pulse width random effect
  SS_DrawRandomEffects sampler(1.0, 500, 25000, 0.35, true, false, 5000);
  double got = sampler.posterior_function(&pulse, prop, &pat);

  // Closed-form prior + Jacobian log-ratio.
  const double logold = std::log(theta);
  const double lognew = std::log(prop);
  const double prior_ratio =
      kappa * (0.5*(logold-mu)*(logold-mu) - 0.5*(lognew-mu)*(lognew-mu))
              / (sigma*sigma)
      + (logold - lognew);

  // Likelihood delta, evaluated exactly as the sampler does.
  const double clike = pat.likelihood(false);
  pulse.width = prop;
  const double plike = pat.likelihood(false);
  pulse.width = theta;  // reset

  const double expected = prior_ratio + (plike - clike);

  SECTION( "matches the analytic log-normal log-ratio" ) {
    REQUIRE( got == Approx(expected) );
  }

  SECTION( "differs from the natural-scale parameterization" ) {
    // The same inputs under the natural-scale prior give a different ratio,
    // confirming the lognormal branch is actually taken.
    pat.lognormal_pulses = false;
    double got_natural = sampler.posterior_function(&pulse, prop, &pat);
    REQUIRE( got != Approx(got_natural) );
  }
}
