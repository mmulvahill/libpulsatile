context("Single-subject fit")

test_that("student_t_pulses = FALSE holds every pulse t-scale (kappa) at 1", {
  skip_on_cran()

  set.seed(2025)
  n_obs <- 24
  dat <- data.frame(
    time = seq(0, 1440, length.out = n_obs),
    concentration = abs(rnorm(n_obs, mean = 3, sd = 1))
  )

  spec <- pulse_spec(student_t_pulses = FALSE, prior_mean_pulse_count = 4)

  fit <- fit_pulse(data = dat, spec = spec,
                   iters = 1000, thin = 10, burnin = 500, verbose = FALSE)

  # Under Gaussian random effects the per-pulse t-scale is fixed at 1 for every
  # pulse (never drawn; initialized to 1 in the birth-death). The single-subject
  # pulse chain reports these as eta_mass / eta_width.
  expect_true("eta_mass" %in% names(fit$pulse_chain))
  expect_true("eta_width" %in% names(fit$pulse_chain))
  expect_true(all(fit$pulse_chain$eta_mass == 1))
  expect_true(all(fit$pulse_chain$eta_width == 1))

  # width_sd is the PR's primary motivation: its draw was commented out and the
  # parameter sat frozen. It is sampled independently of the t-scale toggle, so
  # assert here that both pulse-to-pulse SDs actually move -- a frozen chain
  # (re-introduced `// draw_sd_widths...`) would be constant, giving sd == 0.
  expect_true(stats::sd(fit$patient_chain$width_sd) > 0)
  expect_true(stats::sd(fit$patient_chain$mass_sd) > 0)
})
