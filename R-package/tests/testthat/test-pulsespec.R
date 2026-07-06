context("pulse_spec")

test_that("pulse_spec() defaults to Student-t pulse random effects", {
  spec <- pulse_spec()
  expect_true("student_t_pulses" %in% names(spec$priors))
  expect_true(isTRUE(spec$priors$student_t_pulses))
})

test_that("pulse_spec() carries a Gaussian (kappa=1) request through the priors", {
  spec <- pulse_spec(student_t_pulses = FALSE)
  expect_false(isTRUE(spec$priors$student_t_pulses))
})

test_that("pulse_spec() validates student_t_pulses", {
  expect_error(pulse_spec(student_t_pulses = "yes"), "single logical")
  expect_error(pulse_spec(student_t_pulses = c(TRUE, FALSE)), "single logical")
  expect_error(pulse_spec(student_t_pulses = NA), "single logical")
})
