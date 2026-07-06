#-------------------------------------------------------------------------------
# test-joint.R - Tests for joint hormone model wrapper functions
#-------------------------------------------------------------------------------

test_that("joint_spec() creates valid specification", {
  spec <- joint_spec()

  # Check class
  expect_s3_class(spec, "joint_spec")

  # Check structure
  expect_true("location_prior" %in% names(spec))
  expect_true("driver_priors" %in% names(spec))
  expect_true("response_priors" %in% names(spec))
  expect_true("association_priors" %in% names(spec))
  expect_true("proposal_variances" %in% names(spec))
  expect_true("driver_starting_values" %in% names(spec))
  expect_true("response_starting_values" %in% names(spec))
  expect_true("association_starting_values" %in% names(spec))

  # Check driver priors -- default is now the lognormal (log-scale) column.
  expect_equal(spec$driver_priors$mass_mean, 1.2)
  expect_equal(spec$driver_priors$pulse_count, 12)

  # Check response priors
  expect_equal(spec$response_priors$mass_mean, 1.2)
  expect_equal(spec$response_priors$pulse_count, 12)

  # Check association priors
  expect_equal(spec$association_priors$log_rho_mean, 0.0)
  expect_equal(spec$association_priors$log_nu_mean, 3.0)

  # Check starting values (log-scale defaults)
  expect_equal(spec$driver_starting_values$mass_mean, 1.0)
  expect_equal(spec$response_starting_values$mass_mean, 1.0)
  expect_equal(spec$association_starting_values$rho, 1.0)
  expect_equal(spec$association_starting_values$nu, 20.0)
})


test_that("joint_spec() defaults to lognormal / uniform and converts SD starts", {
  spec <- joint_spec()
  expect_true(isTRUE(spec$driver_priors$lognormal_pulses))
  expect_true(isTRUE(spec$driver_priors$uniform_sd_prior))
  expect_true(isTRUE(spec$response_priors$lognormal_pulses))
  expect_true(isTRUE(spec$response_priors$uniform_sd_prior))
  # Log-scale prior means and starting values.
  expect_equal(spec$driver_priors$width_mean, 3.5)
  expect_equal(spec$driver_starting_values$width_mean, 3.0)
  # Starting SD converted to the log scale (0.7) so it sits inside the C++
  # Uniform(0, 10) support -- the natural-scale 35 would be outside it.
  expect_equal(spec$driver_starting_values$width_sd, 0.7)
  expect_equal(spec$response_starting_values$width_sd, 0.7)
  # Log-scale proposal variances.
  expect_true(spec$proposal_variances$driver_width_mean <= 1)
  expect_true(spec$proposal_variances$driver_pulse_width <= 1)
  expect_error(joint_spec(pulse_distribution = "banana"), "should be one of")
})


test_that("joint_spec() truncnorm/half_cauchy reproduces the pre-change spec", {
  spec <- joint_spec(pulse_distribution = "truncnorm",
                     sd_prior = "half_cauchy",
                     student_t_pulses = TRUE)
  expect_equal(spec$driver_priors$mass_mean, 3.5)
  expect_equal(spec$driver_priors$width_mean, 42)
  expect_equal(spec$driver_priors$width_variance, 1000)
  expect_equal(spec$driver_priors$mass_sd_param, 5)
  expect_equal(spec$driver_priors$width_sd_param, 5)
  expect_equal(spec$response_priors$mass_mean, 3.5)
  expect_equal(spec$response_priors$width_mean, 42)
  expect_equal(spec$driver_starting_values$mass_mean, 3.5)
  expect_equal(spec$driver_starting_values$width_mean, 42)
  expect_equal(spec$driver_starting_values$mass_sd, 1.6)
  expect_equal(spec$driver_starting_values$width_sd, 35)
  expect_equal(spec$response_starting_values$width_sd, 35)
  expect_equal(spec$proposal_variances$driver_width_mean, 3700)
  expect_equal(spec$proposal_variances$driver_pulse_width, 15000)
  expect_equal(spec$proposal_variances$driver_width_sd, 4000)
  expect_equal(spec$proposal_variances$response_width_mean, 3700)
})


test_that("joint_spec() validates input parameters", {
  # Invalid location prior
  expect_error(
    joint_spec(location_prior_type = "invalid"),
    "location_prior_type must be 'strauss'"
  )

  # Negative pulse count
  expect_error(
    joint_spec(prior_driver_pulse_count = -1),
    "prior pulse counts must be > 0"
  )

  # Invalid gamma range
  expect_error(
    joint_spec(prior_driver_location_gamma = 1.5),
    "prior_driver_location_gamma must be in \\[0,1\\]"
  )

  # Negative range
  expect_error(
    joint_spec(prior_response_location_range = -1),
    "prior_response_location_range must be >= 0"
  )

  # Invalid starting values
  expect_error(
    joint_spec(sv_rho = -1),
    "All starting value SD, variance, and association parameters must be > 0"
  )

  # student_t_pulses must be a single logical
  expect_error(
    joint_spec(student_t_pulses = "no"),
    "single logical"
  )
})


test_that("joint_spec() carries student_t_pulses into both hormones", {
  # Default is now Gaussian (student_t_pulses = FALSE).
  spec <- joint_spec()
  expect_false(isTRUE(spec$driver_priors$student_t_pulses))
  expect_false(isTRUE(spec$response_priors$student_t_pulses))

  spec_t <- joint_spec(student_t_pulses = TRUE)
  expect_true(isTRUE(spec_t$driver_priors$student_t_pulses))
  expect_true(isTRUE(spec_t$response_priors$student_t_pulses))
})


test_that("joint_spec() can be customized", {
  spec <- joint_spec(
    prior_driver_mass_mean = 4.0,
    prior_response_mass_mean = 2.5,
    prior_log_rho_mean = 0.5,
    prior_log_nu_mean = 2.5,
    sv_rho = 2.0,
    sv_nu = 15.0
  )

  expect_equal(spec$driver_priors$mass_mean, 4.0)
  expect_equal(spec$response_priors$mass_mean, 2.5)
  expect_equal(spec$association_priors$log_rho_mean, 0.5)
  expect_equal(spec$association_priors$log_nu_mean, 2.5)
  expect_equal(spec$association_starting_values$rho, 2.0)
  expect_equal(spec$association_starting_values$nu, 15.0)
})


test_that("print.joint_spec() works", {
  spec <- joint_spec()

  # Should not error
  expect_output(print(spec), "Bayesian Joint Hormone Model Specification")
  expect_output(print(spec), "Driver Hormone:")
  expect_output(print(spec), "Response Hormone:")
  expect_output(print(spec), "Association")
})


test_that("fit_pulse_joint() validates inputs", {
  spec <- joint_spec()

  # Invalid spec
  expect_error(
    fit_pulse_joint(
      driver_data = data.frame(time = 1:10, concentration = rnorm(10)),
      response_data = data.frame(time = 1:10, concentration = rnorm(10)),
      spec = list()
    ),
    "spec must be a joint_spec object"
  )

  # Burnin >= iters
  expect_error(
    fit_pulse_joint(
      driver_data = data.frame(time = 1:10, concentration = rnorm(10)),
      response_data = data.frame(time = 1:10, concentration = rnorm(10)),
      spec = spec,
      iters = 100,
      burnin = 100
    ),
    "burnin must be < iters"
  )

  # Empty data
  expect_error(
    fit_pulse_joint(
      driver_data = data.frame(time = numeric(0), concentration = numeric(0)),
      response_data = data.frame(time = 1:10, concentration = rnorm(10)),
      spec = spec
    ),
    "Driver and response data must have at least one observation"
  )
})


test_that("fit_pulse_joint() runs with simulated data (quick test)", {
  skip_on_cran()  # Skip on CRAN due to time constraints

  set.seed(2025)

  # Simulate simple driver and response data
  n_obs <- 12
  driver_data <- data.frame(
    time = seq(0, 1440, length.out = n_obs),
    concentration = abs(rnorm(n_obs, mean = 3, sd = 1))
  )

  response_data <- data.frame(
    time = seq(0, 1440, length.out = n_obs),
    concentration = abs(rnorm(n_obs, mean = 2.5, sd = 0.8))
  )

  # Create spec
  spec <- joint_spec(
    prior_driver_pulse_count = 4,
    prior_response_pulse_count = 4
  )

  # Fit model (very short run for testing)
  fit <- fit_pulse_joint(
    driver_data = driver_data,
    response_data = response_data,
    spec = spec,
    iters = 1000,
    thin = 10,
    burnin = 500,
    verbose = FALSE
  )

  # Check output structure
  expect_s3_class(fit, "joint_fit")
  expect_equal(fit$model, "joint")
  expect_true("driver_chain" %in% names(fit))
  expect_true("response_chain" %in% names(fit))
  expect_true("association_chain" %in% names(fit))
  expect_true("driver_pulse_chain" %in% names(fit))
  expect_true("response_pulse_chain" %in% names(fit))

  # Check chains have expected columns
  expect_true("num_pulses" %in% names(fit$driver_chain))
  expect_true("num_pulses" %in% names(fit$response_chain))
  expect_true("rho" %in% names(fit$association_chain))
  expect_true("nu" %in% names(fit$association_chain))

  # Check chain lengths
  expected_samples <- (1000 - 500) / 10
  expect_equal(nrow(fit$driver_chain), expected_samples)
  expect_equal(nrow(fit$response_chain), expected_samples)
  expect_equal(nrow(fit$association_chain), expected_samples)

  # Check parameters are numeric
  expect_true(is.numeric(fit$association_chain$rho))
  expect_true(is.numeric(fit$association_chain$nu))
  expect_true(all(fit$association_chain$rho > 0))
  expect_true(all(fit$association_chain$nu > 0))

  # Pulse-to-pulse SDs are now sampled (previously frozen) and surfaced in the
  # 7-column fixed-effects chains for both hormones.
  expect_true(all(c("mass_sd", "width_sd") %in% names(fit$driver_chain)))
  expect_true(all(c("mass_sd", "width_sd") %in% names(fit$response_chain)))
  # Both pulse-to-pulse SDs must actually move -- a frozen chain (the bug this
  # fixes) would be a constant, giving sd == 0.
  expect_true(stats::sd(fit$driver_chain$mass_sd) > 0)
  expect_true(stats::sd(fit$response_chain$mass_sd) > 0)
  expect_true(stats::sd(fit$driver_chain$width_sd) > 0)
  expect_true(stats::sd(fit$response_chain$width_sd) > 0)
})


test_that("print.joint_fit() works", {
  skip_on_cran()

  set.seed(2025)

  # Quick fit
  driver_data <- data.frame(
    time = seq(0, 1440, length.out = 12),
    concentration = abs(rnorm(12, mean = 3, sd = 1))
  )

  response_data <- data.frame(
    time = seq(0, 1440, length.out = 12),
    concentration = abs(rnorm(12, mean = 2.5, sd = 0.8))
  )

  spec <- joint_spec()

  fit <- fit_pulse_joint(
    driver_data = driver_data,
    response_data = response_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )

  # Should not error
  expect_output(print(fit), "Joint Hormone Model Fit")
  expect_output(print(fit), "Driver Hormone Parameters:")
  expect_output(print(fit), "Response Hormone Parameters:")
  expect_output(print(fit), "Association Parameters")
})


test_that("summary.joint_fit() works", {
  skip_on_cran()

  set.seed(2025)

  # Quick fit
  driver_data <- data.frame(
    time = seq(0, 1440, length.out = 12),
    concentration = abs(rnorm(12, mean = 3, sd = 1))
  )

  response_data <- data.frame(
    time = seq(0, 1440, length.out = 12),
    concentration = abs(rnorm(12, mean = 2.5, sd = 0.8))
  )

  spec <- joint_spec()

  fit <- fit_pulse_joint(
    driver_data = driver_data,
    response_data = response_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )

  # Should not error
  expect_output(summary(fit), "Joint Hormone Model Summary")
  expect_output(summary(fit), "Driver Hormone Parameter Estimates")
  expect_output(summary(fit), "Response Hormone Parameter Estimates")
  expect_output(summary(fit), "Association Parameter Estimates")
  expect_output(summary(fit), "Coupling Strength")
  expect_output(summary(fit), "Coupling Temporal Spread")
})


test_that("fit_pulse_joint() handles pulse_sim objects", {
  skip_on_cran()

  # Create pulse_sim objects
  set.seed(2025)
  driver_sim <- simulate_pulse(num_obs = 12, interval = 120)
  response_sim <- simulate_pulse(num_obs = 12, interval = 120)

  spec <- joint_spec()

  fit <- fit_pulse_joint(
    driver_data = driver_sim,
    response_data = response_sim,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )

  expect_s3_class(fit, "joint_fit")
  expect_equal(fit$model, "joint")
})


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
