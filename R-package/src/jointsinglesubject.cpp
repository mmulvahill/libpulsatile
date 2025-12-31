#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/bp_datastructures.h>
#include <bp_datastructures/associationparameters.h>
#include <bp_datastructures/chains.h>
#include <bpmod_joint/joint_mcmc_iteration.h>

using namespace Rcpp;


//
// jointsinglesubject.cpp
//   Rcpp export function for joint driver-response hormone model
//
// Author: Matt Mulvahill
// Created: 12/30/24
//
// This function implements the joint model for two coupled hormones:
//   - Driver hormone (e.g., LH): Standard pulsatile model
//   - Response hormone (e.g., FSH): Pulsatile model with coupling to driver
//
// The coupling is modeled through lambda (λ), which modulates the birth rate
// of response pulses based on driver pulse locations.
//


//' Joint Single Subject Model
//'
//' Fit a joint Bayesian model for two coupled pulsatile hormones measured
//' from a single subject. Driver pulses influence response pulse occurrence
//' through a coupling kernel parameterized by rho and nu.
//'
//' @param driver_concentration Numeric vector of driver hormone concentrations
//' @param driver_time Numeric vector of driver hormone measurement times
//' @param response_concentration Numeric vector of response hormone concentrations
//' @param response_time Numeric vector of response hormone measurement times
//' @param location_prior Character string: "strauss" (only option currently)
//' @param driver_priors List of prior parameters for driver hormone
//' @param response_priors List of prior parameters for response hormone
//' @param association_priors List of prior parameters for coupling (rho, nu)
//' @param proposalvars List of proposal variances for all parameters
//' @param driver_startingvals List of starting values for driver parameters
//' @param response_startingvals List of starting values for response parameters
//' @param association_startingvals List of starting values for coupling parameters
//' @param mcmc_iterations Integer: total MCMC iterations to run
//' @param thin Integer: thinning interval for saving samples
//' @param burnin Integer: number of initial iterations to discard
//' @param verbose Logical: print progress messages
//' @param pv_adjust_iter Integer: interval for proposal variance adaptation
//' @param pv_adjust_max_iter Integer: iteration to stop adaptation
//' @param bivariate_pv_target_ratio Numeric: target acceptance rate for bivariate proposals
//' @param univariate_pv_target_ratio Numeric: target acceptance rate for univariate proposals
//'
//' @return List containing MCMC chains for driver, response, and association parameters
//' @export
// [[Rcpp::export]]
Rcpp::List jointsinglesubject_(Rcpp::NumericVector driver_concentration,
                                Rcpp::NumericVector driver_time,
                                Rcpp::NumericVector response_concentration,
                                Rcpp::NumericVector response_time,
                                Rcpp::CharacterVector location_prior,
                                Rcpp::List driver_priors,
                                Rcpp::List response_priors,
                                Rcpp::List association_priors,
                                Rcpp::List proposalvars,
                                Rcpp::List driver_startingvals,
                                Rcpp::List response_startingvals,
                                Rcpp::List association_startingvals,
                                int mcmc_iterations,
                                int thin,
                                int burnin,
                                bool verbose,
                                int pv_adjust_iter,
                                int pv_adjust_max_iter,
                                double bivariate_pv_target_ratio,
                                double univariate_pv_target_ratio) {

  Rcpp::RNGScope rng_scope;

  //
  // Initialize driver hormone patient
  //
  PatientData driver_data(driver_time, driver_concentration);
  PatientPriors driver_patient_priors(
    driver_priors["baseline_mean"],
    driver_priors["baseline_variance"],
    driver_priors["halflife_mean"],
    driver_priors["halflife_variance"],
    driver_priors["mass_mean"],
    driver_priors["mass_variance"],
    driver_priors["width_mean"],
    driver_priors["width_variance"],
    driver_priors["mass_sd_param"],
    driver_priors["width_sd_param"],
    driver_priors["error_alpha"],
    driver_priors["error_beta"],
    driver_priors["pulse_count"],
    driver_priors["strauss_repulsion"],
    driver_priors["strauss_repulsion_range"]);
  PatientEstimates driver_estimates(
    driver_startingvals["baseline"],
    driver_startingvals["halflife"],
    driver_startingvals["errorsq"],
    driver_startingvals["mass_mean"],
    driver_startingvals["width_mean"],
    driver_startingvals["mass_sd"],
    driver_startingvals["width_sd"]);
  Patient driver_patient(driver_data, driver_patient_priors, driver_estimates);

  //
  // Initialize response hormone patient
  //
  PatientData response_data(response_time, response_concentration);
  PatientPriors response_patient_priors(
    response_priors["baseline_mean"],
    response_priors["baseline_variance"],
    response_priors["halflife_mean"],
    response_priors["halflife_variance"],
    response_priors["mass_mean"],
    response_priors["mass_variance"],
    response_priors["width_mean"],
    response_priors["width_variance"],
    response_priors["mass_sd_param"],
    response_priors["width_sd_param"],
    response_priors["error_alpha"],
    response_priors["error_beta"],
    response_priors["pulse_count"],
    response_priors["strauss_repulsion"],
    response_priors["strauss_repulsion_range"]);
  PatientEstimates response_estimates(
    response_startingvals["baseline"],
    response_startingvals["halflife"],
    response_startingvals["errorsq"],
    response_startingvals["mass_mean"],
    response_startingvals["width_mean"],
    response_startingvals["mass_sd"],
    response_startingvals["width_sd"]);
  Patient response_patient(response_data, response_patient_priors, response_estimates);

  //
  // Initialize association parameters
  //
  AssociationPriors assoc_priors(
    Rcpp::as<double>(association_priors["log_rho_mean"]),
    Rcpp::as<double>(association_priors["log_rho_var"]),
    Rcpp::as<double>(association_priors["log_nu_mean"]),
    Rcpp::as<double>(association_priors["log_nu_var"])
  );

  AssociationEstimates assoc_est(
    Rcpp::as<double>(association_startingvals["rho"]),
    Rcpp::as<double>(association_startingvals["nu"])
  );

  //
  // Initialize MCMC samplers
  //
  std::string loc_prior_str = Rcpp::as<std::string>(location_prior[0]);
  JointSamplers samplers(proposalvars, pv_adjust_iter, pv_adjust_max_iter,
                        bivariate_pv_target_ratio, univariate_pv_target_ratio,
                        verbose, mcmc_iterations, loc_prior_str);

  //
  // Initialize output chains
  //
  int num_saved = (mcmc_iterations - burnin) / thin;

  // Driver chains
  arma::mat driver_fixed_effects_chain(num_saved, 5, arma::fill::zeros);  // baseline, halflife, mass_mean, width_mean, errorsq
  MatrixVector driver_pulse_chains;

  // Response chains
  arma::mat response_fixed_effects_chain(num_saved, 5, arma::fill::zeros);  // baseline, halflife, mass_mean, width_mean, errorsq
  MatrixVector response_pulse_chains;

  // Association chains
  arma::mat association_chain(num_saved, 2, arma::fill::zeros);  // rho, nu

  //
  // Initialize lambda values for response pulses
  //
  update_lambda_values(driver_patient.pulses, response_patient.pulses, assoc_est);

  //
  // Main MCMC loop
  //
  if (verbose) Rcpp::Rcout << "Starting MCMC iterations..." << std::endl;

  int save_index = 0;

  for (int iter = 0; iter < mcmc_iterations; iter++) {

    // Run one MCMC iteration
    joint_mcmc_iteration(&driver_patient, &response_patient,
                        &assoc_est, &assoc_priors,
                        samplers, iter);

    // Save samples after burnin, applying thinning
    if (iter >= burnin && (iter - burnin) % thin == 0) {

      // Driver fixed effects (baseline, halflife, mass_mean, width_mean, errorsq)
      driver_fixed_effects_chain(save_index, 0) = driver_patient.estimates.baseline_halflife(0);
      driver_fixed_effects_chain(save_index, 1) = driver_patient.estimates.baseline_halflife(1);
      driver_fixed_effects_chain(save_index, 2) = driver_patient.estimates.mass_mean;
      driver_fixed_effects_chain(save_index, 3) = driver_patient.estimates.width_mean;
      driver_fixed_effects_chain(save_index, 4) = driver_patient.estimates.errorsq;

      // Driver pulses
      arma::mat driver_pulses_matrix(driver_patient.get_pulsecount(), 5);
      int pulse_idx = 0;
      for (const auto& pulse : driver_patient.pulses) {
        driver_pulses_matrix.row(pulse_idx++) = pulse.get_vector_of_values();
      }
      driver_pulse_chains.push_back(driver_pulses_matrix);

      // Response fixed effects (baseline, halflife, mass_mean, width_mean, errorsq)
      response_fixed_effects_chain(save_index, 0) = response_patient.estimates.baseline_halflife(0);
      response_fixed_effects_chain(save_index, 1) = response_patient.estimates.baseline_halflife(1);
      response_fixed_effects_chain(save_index, 2) = response_patient.estimates.mass_mean;
      response_fixed_effects_chain(save_index, 3) = response_patient.estimates.width_mean;
      response_fixed_effects_chain(save_index, 4) = response_patient.estimates.errorsq;

      // Response pulses (including lambda)
      arma::mat response_pulses_matrix(response_patient.get_pulsecount(), 6);
      pulse_idx = 0;
      for (const auto& pulse : response_patient.pulses) {
        response_pulses_matrix.row(pulse_idx++) = pulse.get_vector_of_values_with_lambda();
      }
      response_pulse_chains.push_back(response_pulses_matrix);

      // Association parameters
      association_chain(save_index, 0) = assoc_est.get_rho();
      association_chain(save_index, 1) = assoc_est.get_nu();

      save_index++;
    }

    // Progress reporting
    if (verbose && (iter + 1) % 1000 == 0) {
      Rcpp::Rcout << "Iteration " << (iter + 1) << " / " << mcmc_iterations << std::endl;
    }
  }

  if (verbose) Rcpp::Rcout << "MCMC complete!" << std::endl;

  //
  // Return results
  //
  return Rcpp::List::create(
    Rcpp::Named("driver_fixed_effects") = driver_fixed_effects_chain,
    Rcpp::Named("driver_pulses") = driver_pulse_chains,
    Rcpp::Named("response_fixed_effects") = response_fixed_effects_chain,
    Rcpp::Named("response_pulses") = response_pulse_chains,
    Rcpp::Named("association") = association_chain,
    Rcpp::Named("driver_colnames") = Rcpp::CharacterVector::create(
      "baseline", "halflife", "mass_mean", "width_mean", "errorsq"),
    Rcpp::Named("response_colnames") = Rcpp::CharacterVector::create(
      "baseline", "halflife", "mass_mean", "width_mean", "errorsq"),
    Rcpp::Named("driver_pulse_colnames") = Rcpp::CharacterVector::create(
      "time", "mass", "width", "tvarscale_mass", "tvarscale_width"),
    Rcpp::Named("response_pulse_colnames") = Rcpp::CharacterVector::create(
      "time", "mass", "width", "tvarscale_mass", "tvarscale_width", "lambda"),
    Rcpp::Named("association_colnames") = Rcpp::CharacterVector::create(
      "rho", "nu")
  );
}
