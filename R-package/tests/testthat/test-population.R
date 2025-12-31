context("Population Model Tests")

test_that("population_spec() creates valid specification object", {
  spec <- population_spec()
  
  expect_s3_class(spec, "population_spec")
  expect_true(is.list(spec))
  expect_true("location_prior" %in% names(spec))
  expect_true("population_priors" %in% names(spec))
  expect_true("proposal_variances" %in% names(spec))
  expect_true("population_starting_values" %in% names(spec))
  
  # Check location prior is valid
  expect_true(spec$location_prior %in% c("strauss", "order-statistic"))
  
  # Check population priors has all required fields
  required_priors <- c("mass_mean_mean", "mass_mean_var", 
                       "width_mean_mean", "width_mean_var",
                       "baseline_mean_mean", "baseline_mean_var",
                       "halflife_mean_mean", "halflife_mean_var",
                       "mass_mean_sd_max", "width_mean_sd_max",
                       "baseline_sd_max", "halflife_sd_max",
                       "mass_sd_max", "width_sd_max",
                       "error_alpha", "error_beta", "pulse_count")
  expect_true(all(required_priors %in% names(spec$population_priors)))
  
  # Check all SD max values are positive
  sd_maxes <- c(spec$population_priors$mass_mean_sd_max,
                spec$population_priors$width_mean_sd_max,
                spec$population_priors$baseline_sd_max,
                spec$population_priors$halflife_sd_max,
                spec$population_priors$mass_sd_max,
                spec$population_priors$width_sd_max)
  expect_true(all(sd_maxes > 0))
})


test_that("population_spec() validates inputs correctly", {
  # Invalid location prior
  expect_error(
    population_spec(location_prior_type = "invalid"),
    "should be one of"
  )
  
  # Invalid pulse count
  expect_error(
    population_spec(prior_mean_pulse_count = -1),
    "must be > 0"
  )
  
  # Invalid SD max (negative)
  expect_error(
    population_spec(prior_mass_mean_sd_max = -1),
    "must be > 0"
  )
  
  # Strauss prior requires gamma and range
  expect_error(
    population_spec(location_prior_type = "strauss",
                    prior_location_gamma = NULL),
    "required"
  )
})


test_that("population_spec() print method works", {
  spec <- population_spec()
  expect_output(print(spec), "Bayesian Population Model")
  expect_output(print(spec), "Location prior")
  expect_output(print(spec), "Population Parameters")
})


test_that("fit_pulse_population() runs with simulated data", {
  skip_on_cran()
  skip_on_ci()
  
  # Create simulated data for 3 subjects
  set.seed(123)
  sim_data <- lapply(1:3, function(i) {
    simulate_pulse(num_obs = 24, interval = 10, 
                   ipi_mean = 12, mass_mean = 3.5, 
                   width_mean = 35)
  })
  
  # Create specification with small iterations for testing
  spec <- population_spec(
    prior_mass_mean_mean = 3.5,
    prior_width_mean_mean = 35,
    prior_baseline_mean_mean = 2.6,
    prior_halflife_mean_mean = 45,
    prior_mean_pulse_count = 12
  )
  
  # Fit model (short run for testing)
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 1000,
    thin = 10,
    burnin = 500,
    verbose = FALSE
  )
  
  # Check output structure
  expect_s3_class(fit, "population_fit")
  expect_equal(fit$model, "population")
  expect_equal(fit$num_subjects, 3)
  expect_true(is.data.frame(fit$population_chain))
  expect_true(is.list(fit$subject_chains))
  expect_equal(length(fit$subject_chains), 3)
  expect_true(is.list(fit$pulse_chains))
  expect_equal(length(fit$pulse_chains), 3)
  
  # Check population chain has correct columns
  expected_cols <- c("mass_mean", "width_mean", "baseline_mean", 
                     "halflife_mean", "mass_mean_sd", "width_mean_sd",
                     "baseline_sd", "halflife_sd", "mass_sd", "width_sd",
                     "errorsq")
  expect_true(all(expected_cols %in% names(fit$population_chain)))
  
  # Check number of saved samples is correct
  expected_samples <- (1000 - 500) / 10
  expect_equal(nrow(fit$population_chain), expected_samples)
  
  # Check subject chains exist and have data
  for (i in 1:3) {
    expect_true(is.data.frame(fit$subject_chains[[i]]))
    expect_equal(nrow(fit$subject_chains[[i]]), expected_samples)
    expect_true("baseline" %in% names(fit$subject_chains[[i]]))
    expect_true("halflife" %in% names(fit$subject_chains[[i]]))
    expect_true("mass_mean" %in% names(fit$subject_chains[[i]]))
    expect_true("width_mean" %in% names(fit$subject_chains[[i]]))
  }
  
  # Check all parameters are numeric and finite
  expect_true(all(is.finite(unlist(fit$population_chain))))
  expect_true(all(sapply(fit$subject_chains, function(sc) all(is.finite(unlist(sc))))))
})


test_that("fit_pulse_population() handles data frame with subject_id", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(456)
  
  # Create combined data frame
  sim_list <- lapply(1:2, function(i) {
    sim <- simulate_pulse(num_obs = 20, interval = 10)
    data.frame(
      subject_id = i,
      time = sim$data$time,
      concentration = sim$data$concentration
    )
  })
  combined_data <- do.call(rbind, sim_list)
  
  spec <- population_spec()
  
  # Fit using subject_id
  fit <- fit_pulse_population(
    data = combined_data,
    subject_id = "subject_id",
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "population_fit")
  expect_equal(fit$num_subjects, 2)
  expect_equal(length(fit$subject_chains), 2)
})


test_that("fit_pulse_population() validates inputs", {
  spec <- population_spec()
  
  # Invalid spec
  expect_error(
    fit_pulse_population(data = list(), spec = "not_a_spec"),
    "must be a population_spec"
  )
  
  # burnin >= iters
  expect_error(
    fit_pulse_population(data = list(), spec = spec, 
                         iters = 100, burnin = 150),
    "burnin must be < iters"
  )
  
  # No subjects
  expect_error(
    fit_pulse_population(data = list(), spec = spec),
    "No subjects found"
  )
})


test_that("fit_pulse_population() custom subject starting values work", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(789)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 20))
  
  # Custom starting values
  custom_starts <- list(
    list(baseline = 2.5, halflife = 40, mass_mean = 3.0, width_mean = 30),
    list(baseline = 2.8, halflife = 50, mass_mean = 4.0, width_mean = 40)
  )
  
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    subject_starts = custom_starts,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "population_fit")
  expect_equal(fit$num_subjects, 2)
})


test_that("print.population_fit method works", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(101)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 15))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  expect_output(print(fit), "Population Pulse Model Fit")
  expect_output(print(fit), "Number of subjects")
  expect_output(print(fit), "MCMC iterations")
})


test_that("summary.population_fit method works", {
  skip_on_cran()
  skip_on_ci()
  
  set.seed(102)
  sim_data <- lapply(1:2, function(i) simulate_pulse(num_obs = 15))
  spec <- population_spec()
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 500,
    thin = 5,
    burnin = 250,
    verbose = FALSE
  )
  
  expect_output(summary(fit), "Population Model Summary")
  expect_output(summary(fit), "Population-level parameter estimates")
  expect_output(summary(fit), "Subject-level parameter estimates")
})


test_that("fit_pulse_population() produces reasonable parameter estimates", {
  skip_on_cran()
  skip_on_ci()
  
  # Create data with known parameters
  set.seed(2024)
  true_mass_mean <- 3.5
  true_width_mean <- 35
  true_baseline <- 2.6
  true_halflife <- 45
  
  sim_data <- lapply(1:5, function(i) {
    simulate_pulse(
      num_obs = 48,
      interval = 10,
      mass_mean = true_mass_mean,
      width_mean = true_width_mean,
      mass_sd = 1.0,
      width_sd = 5,
      constant_baseline = true_baseline,
      constant_halflife = true_halflife,
      ipi_mean = 12
    )
  })
  
  spec <- population_spec(
    prior_mass_mean_mean = 3.5,
    prior_width_mean_mean = 35,
    prior_baseline_mean_mean = 2.6,
    prior_halflife_mean_mean = 45
  )
  
  fit <- fit_pulse_population(
    data = sim_data,
    spec = spec,
    iters = 5000,
    thin = 10,
    burnin = 2500,
    verbose = FALSE
  )
  
  # Check posterior means are reasonably close to true values
  # (allowing wide tolerance given short MCMC run)
  pop_means <- colMeans(fit$population_chain)
  
  expect_true(abs(pop_means["mass_mean"] - true_mass_mean) < 2.0,
              info = sprintf("mass_mean estimate %f vs true %f", 
                           pop_means["mass_mean"], true_mass_mean))
  
  expect_true(abs(pop_means["width_mean"] - true_width_mean) < 15,
              info = sprintf("width_mean estimate %f vs true %f",
                           pop_means["width_mean"], true_width_mean))
  
  expect_true(abs(pop_means["baseline_mean"] - true_baseline) < 1.0,
              info = sprintf("baseline_mean estimate %f vs true %f",
                           pop_means["baseline_mean"], true_baseline))
  
  expect_true(abs(pop_means["halflife_mean"] - true_halflife) < 20,
              info = sprintf("halflife_mean estimate %f vs true %f",
                           pop_means["halflife_mean"], true_halflife))
  
  # Check that SDs are positive
  expect_true(all(fit$population_chain$mass_mean_sd > 0))
  expect_true(all(fit$population_chain$width_mean_sd > 0))
  expect_true(all(fit$population_chain$mass_sd > 0))
  expect_true(all(fit$population_chain$width_sd > 0))
})
