#-------------------------------------------------------------------------------
# joint_spec.R - Functions for creating a joint hormone model specification
#-------------------------------------------------------------------------------

#' joint_spec
#'
#' Generates a joint_spec object -- the specification object required for
#' fitting a joint hormone model with fit_pulse_joint().
#'
#' The joint model fits two coupled pulsatile hormones where driver hormone
#' pulses influence the intensity of response hormone pulse occurrence.
#'
#' @param location_prior_type Takes on value "strauss" (only option currently).
#'   "strauss" uses the Strauss interacting point-process as a prior and requires
#'   specification of "prior_mean_pulse_count", "prior_location_gamma", and
#'   "prior_location_range".
#'
#' @section Driver Hormone Priors:
#' @param prior_driver_mass_mean Prior mean for driver pulse mass
#' @param prior_driver_mass_var Prior variance for driver pulse mass
#' @param prior_driver_width_mean Prior mean for driver pulse width (on variance scale)
#' @param prior_driver_width_var Prior variance for driver pulse width (on variance scale)
#' @param prior_driver_baseline_mean Prior mean for driver baseline
#' @param prior_driver_baseline_var Prior variance for driver baseline
#' @param prior_driver_halflife_mean Prior mean for driver half-life
#' @param prior_driver_halflife_var Prior variance for driver half-life
#' @param prior_driver_error_alpha Gamma shape parameter for driver error
#' @param prior_driver_error_beta Gamma rate parameter for driver error
#' @param prior_driver_sd_mass Upper bound (uniform prior) for driver pulse-to-pulse SD of mass
#' @param prior_driver_sd_width Upper bound (uniform prior) for driver pulse-to-pulse SD of width
#' @param prior_driver_pulse_count Mean of Poisson prior on driver pulse count
#' @param prior_driver_location_gamma Strauss repulsion parameter for driver (0-1, 0=no repulsion)
#' @param prior_driver_location_range Strauss interaction range for driver (in time units)
#'
#' @section Response Hormone Priors:
#' @param prior_response_mass_mean Prior mean for response pulse mass
#' @param prior_response_mass_var Prior variance for response pulse mass
#' @param prior_response_width_mean Prior mean for response pulse width (on variance scale)
#' @param prior_response_width_var Prior variance for response pulse width (on variance scale)
#' @param prior_response_baseline_mean Prior mean for response baseline
#' @param prior_response_baseline_var Prior variance for response baseline
#' @param prior_response_halflife_mean Prior mean for response half-life
#' @param prior_response_halflife_var Prior variance for response half-life
#' @param prior_response_error_alpha Gamma shape parameter for response error
#' @param prior_response_error_beta Gamma rate parameter for response error
#' @param prior_response_sd_mass Upper bound (uniform prior) for response pulse-to-pulse SD of mass
#' @param prior_response_sd_width Upper bound (uniform prior) for response pulse-to-pulse SD of width
#' @param prior_response_pulse_count Mean of Poisson prior on response pulse count
#' @param prior_response_location_gamma Strauss repulsion parameter for response (0-1, 0=no repulsion)
#' @param prior_response_location_range Strauss interaction range for response (in time units)
#'
#' @section Association Priors:
#' @param prior_log_rho_mean Prior mean for log(rho) - coupling strength parameter
#' @param prior_log_rho_var Prior variance for log(rho)
#' @param prior_log_nu_mean Prior mean for log(nu) - coupling temporal spread parameter
#' @param prior_log_nu_var Prior variance for log(nu)
#'
#' @section Driver Starting Values:
#' @param sv_driver_mass_mean Starting value for driver mean pulse mass
#' @param sv_driver_width_mean Starting value for driver mean pulse width
#' @param sv_driver_baseline_mean Starting value for driver baseline
#' @param sv_driver_halflife_mean Starting value for driver half-life
#' @param sv_driver_error_var Starting value for driver error variance
#' @param sv_driver_mass_sd Starting value for driver pulse-to-pulse SD of mass
#' @param sv_driver_width_sd Starting value for driver pulse-to-pulse SD of width
#'
#' @section Response Starting Values:
#' @param sv_response_mass_mean Starting value for response mean pulse mass
#' @param sv_response_width_mean Starting value for response mean pulse width
#' @param sv_response_baseline_mean Starting value for response baseline
#' @param sv_response_halflife_mean Starting value for response half-life
#' @param sv_response_error_var Starting value for response error variance
#' @param sv_response_mass_sd Starting value for response pulse-to-pulse SD of mass
#' @param sv_response_width_sd Starting value for response pulse-to-pulse SD of width
#'
#' @section Association Starting Values:
#' @param sv_rho Starting value for rho (coupling strength)
#' @param sv_nu Starting value for nu (coupling temporal spread, is variance)
#'
#' @section Driver Proposal Variances:
#' @param pv_driver_baseline Proposal variance for driver baseline
#' @param pv_driver_halflife Proposal variance for driver half-life
#' @param pv_driver_mean_pulse_mass Proposal variance for driver mean pulse mass
#' @param pv_driver_mean_pulse_width Proposal variance for driver mean pulse width
#' @param pv_driver_indiv_pulse_mass Proposal variance for driver individual pulse masses
#' @param pv_driver_indiv_pulse_width Proposal variance for driver individual pulse widths
#' @param pv_driver_sd_pulse_mass Proposal variance for driver pulse mass SD
#' @param pv_driver_sd_pulse_width Proposal variance for driver pulse width SD
#' @param pv_driver_sdscale_pulse_mass Proposal variance for driver pulse mass t-distribution scales
#' @param pv_driver_sdscale_pulse_width Proposal variance for driver pulse width t-distribution scales
#' @param pv_driver_pulse_location Proposal variance for driver pulse locations
#'
#' @section Response Proposal Variances:
#' @param pv_response_baseline Proposal variance for response baseline
#' @param pv_response_halflife Proposal variance for response half-life
#' @param pv_response_mean_pulse_mass Proposal variance for response mean pulse mass
#' @param pv_response_mean_pulse_width Proposal variance for response mean pulse width
#' @param pv_response_indiv_pulse_mass Proposal variance for response individual pulse masses
#' @param pv_response_indiv_pulse_width Proposal variance for response individual pulse widths
#' @param pv_response_sd_pulse_mass Proposal variance for response pulse mass SD
#' @param pv_response_sd_pulse_width Proposal variance for response pulse width SD
#' @param pv_response_sdscale_pulse_mass Proposal variance for response pulse mass t-distribution scales
#' @param pv_response_sdscale_pulse_width Proposal variance for response pulse width t-distribution scales
#' @param pv_response_pulse_location Proposal variance for response pulse locations
#'
#' @section Association Proposal Variances:
#' @param pv_log_rho Proposal variance for log(rho)
#' @param pv_log_nu Proposal variance for log(nu)
#'
#' @return A list of class \code{joint_spec} containing:
#'   \item{location_prior}{Type of location prior ("strauss")}
#'   \item{driver_priors}{List of prior parameters for driver hormone}
#'   \item{response_priors}{List of prior parameters for response hormone}
#'   \item{association_priors}{List of prior parameters for coupling}
#'   \item{proposal_variances}{List of proposal variances for all parameters}
#'   \item{driver_starting_values}{List of starting values for driver parameters}
#'   \item{response_starting_values}{List of starting values for response parameters}
#'   \item{association_starting_values}{List of starting values for coupling parameters}
#'
#' @export
#' @keywords pulse joint
#' @examples
#' # Create a joint specification with default values
#' joint_spec <- joint_spec()
#'
#' # Customize for specific hormones (e.g., LH driving FSH)
#' custom_spec <- joint_spec(
#'   prior_driver_mass_mean = 4.0,
#'   prior_response_mass_mean = 2.5,
#'   prior_log_rho_mean = 0.0,  # log(1) = neutral coupling
#'   prior_log_nu_mean = 3.0    # log(20) for temporal spread
#' )
joint_spec <- function(
  # Location prior
  location_prior_type = "strauss",

  # Driver priors
  prior_driver_mass_mean = 3.5,
  prior_driver_mass_var = 100,
  prior_driver_width_mean = 42,
  prior_driver_width_var = 1000,
  prior_driver_baseline_mean = 2.6,
  prior_driver_baseline_var = 100,
  prior_driver_halflife_mean = 45,
  prior_driver_halflife_var = 100,
  prior_driver_error_alpha = 0.0001,
  prior_driver_error_beta = 0.0001,
  prior_driver_sd_mass = 5,
  prior_driver_sd_width = 5,
  prior_driver_pulse_count = 12,
  prior_driver_location_gamma = 0,
  prior_driver_location_range = 40,

  # Response priors
  prior_response_mass_mean = 3.5,
  prior_response_mass_var = 100,
  prior_response_width_mean = 42,
  prior_response_width_var = 1000,
  prior_response_baseline_mean = 2.6,
  prior_response_baseline_var = 100,
  prior_response_halflife_mean = 45,
  prior_response_halflife_var = 100,
  prior_response_error_alpha = 0.0001,
  prior_response_error_beta = 0.0001,
  prior_response_sd_mass = 5,
  prior_response_sd_width = 5,
  prior_response_pulse_count = 12,
  prior_response_location_gamma = 0,
  prior_response_location_range = 40,

  # Association priors (on log scale)
  prior_log_rho_mean = 0.0,    # log(1) = neutral coupling
  prior_log_rho_var = 1.0,
  prior_log_nu_mean = 3.0,     # log(20) for temporal spread
  prior_log_nu_var = 1.0,

  # Driver starting values
  sv_driver_mass_mean = 3.5,
  sv_driver_width_mean = 42,
  sv_driver_baseline_mean = 2.6,
  sv_driver_halflife_mean = 45,
  sv_driver_error_var = 0.005,
  sv_driver_mass_sd = 1.6,
  sv_driver_width_sd = 35,

  # Response starting values
  sv_response_mass_mean = 3.5,
  sv_response_width_mean = 42,
  sv_response_baseline_mean = 2.6,
  sv_response_halflife_mean = 45,
  sv_response_error_var = 0.005,
  sv_response_mass_sd = 1.6,
  sv_response_width_sd = 35,

  # Association starting values
  sv_rho = 1.0,
  sv_nu = 20.0,

  # Driver proposal variances
  pv_driver_baseline = 0.02,
  pv_driver_halflife = 1.5,
  pv_driver_mean_pulse_mass = 6,
  pv_driver_mean_pulse_width = 3700,
  pv_driver_indiv_pulse_mass = 1,
  pv_driver_indiv_pulse_width = 15000,
  pv_driver_sd_pulse_mass = 4.5,
  pv_driver_sd_pulse_width = 4000,
  pv_driver_sdscale_pulse_mass = 4,
  pv_driver_sdscale_pulse_width = 4,
  pv_driver_pulse_location = 65,

  # Response proposal variances
  pv_response_baseline = 0.02,
  pv_response_halflife = 1.5,
  pv_response_mean_pulse_mass = 6,
  pv_response_mean_pulse_width = 3700,
  pv_response_indiv_pulse_mass = 1,
  pv_response_indiv_pulse_width = 15000,
  pv_response_sd_pulse_mass = 4.5,
  pv_response_sd_pulse_width = 4000,
  pv_response_sdscale_pulse_mass = 4,
  pv_response_sdscale_pulse_width = 4,
  pv_response_pulse_location = 65,

  # Association proposal variances (on log scale)
  pv_log_rho = 0.1,
  pv_log_nu = 0.1
) {

  # Input validation
  if (location_prior_type != "strauss") {
    stop("location_prior_type must be 'strauss' (only option currently supported)")
  }

  if (prior_driver_pulse_count <= 0 || prior_response_pulse_count <= 0) {
    stop("prior pulse counts must be > 0")
  }

  # Validate Strauss prior parameters for driver
  if (prior_driver_location_gamma < 0 || prior_driver_location_gamma > 1) {
    stop("prior_driver_location_gamma must be in [0,1]")
  }
  if (prior_driver_location_range < 0) {
    stop("prior_driver_location_range must be >= 0")
  }

  # Validate Strauss prior parameters for response
  if (prior_response_location_gamma < 0 || prior_response_location_gamma > 1) {
    stop("prior_response_location_gamma must be in [0,1]")
  }
  if (prior_response_location_range < 0) {
    stop("prior_response_location_range must be >= 0")
  }

  # Validate positive parameters
  if (any(c(sv_driver_mass_sd, sv_driver_width_sd, sv_driver_error_var,
            sv_response_mass_sd, sv_response_width_sd, sv_response_error_var,
            sv_rho, sv_nu) <= 0)) {
    stop("All starting value SD, variance, and association parameters must be > 0")
  }

  # Create specification object
  spec_obj <- structure(
    list(
      location_prior = location_prior_type,

      # Driver priors
      driver_priors = list(
        baseline_mean = prior_driver_baseline_mean,
        baseline_variance = prior_driver_baseline_var,
        halflife_mean = prior_driver_halflife_mean,
        halflife_variance = prior_driver_halflife_var,
        mass_mean = prior_driver_mass_mean,
        mass_variance = prior_driver_mass_var,
        width_mean = prior_driver_width_mean,
        width_variance = prior_driver_width_var,
        mass_sd_param = prior_driver_sd_mass,
        width_sd_param = prior_driver_sd_width,
        error_alpha = prior_driver_error_alpha,
        error_beta = prior_driver_error_beta,
        pulse_count = prior_driver_pulse_count,
        strauss_repulsion = prior_driver_location_gamma,
        strauss_repulsion_range = prior_driver_location_range
      ),

      # Response priors
      response_priors = list(
        baseline_mean = prior_response_baseline_mean,
        baseline_variance = prior_response_baseline_var,
        halflife_mean = prior_response_halflife_mean,
        halflife_variance = prior_response_halflife_var,
        mass_mean = prior_response_mass_mean,
        mass_variance = prior_response_mass_var,
        width_mean = prior_response_width_mean,
        width_variance = prior_response_width_var,
        mass_sd_param = prior_response_sd_mass,
        width_sd_param = prior_response_sd_width,
        error_alpha = prior_response_error_alpha,
        error_beta = prior_response_error_beta,
        pulse_count = prior_response_pulse_count,
        strauss_repulsion = prior_response_location_gamma,
        strauss_repulsion_range = prior_response_location_range
      ),

      # Association priors
      association_priors = list(
        log_rho_mean = prior_log_rho_mean,
        log_rho_var = prior_log_rho_var,
        log_nu_mean = prior_log_nu_mean,
        log_nu_var = prior_log_nu_var
      ),

      # Proposal variances
      proposal_variances = list(
        # Driver
        driver_mass_mean = pv_driver_mean_pulse_mass,
        driver_width_mean = pv_driver_mean_pulse_width,
        driver_mass_sd = pv_driver_sd_pulse_mass,
        driver_width_sd = pv_driver_sd_pulse_width,
        driver_baseline = pv_driver_baseline,
        driver_halflife = pv_driver_halflife,
        driver_location = pv_driver_pulse_location,
        driver_pulse_mass = pv_driver_indiv_pulse_mass,
        driver_pulse_width = pv_driver_indiv_pulse_width,
        driver_sdscale_pulse_mass = pv_driver_sdscale_pulse_mass,
        driver_sdscale_pulse_width = pv_driver_sdscale_pulse_width,
        # Response
        response_mass_mean = pv_response_mean_pulse_mass,
        response_width_mean = pv_response_mean_pulse_width,
        response_mass_sd = pv_response_sd_pulse_mass,
        response_width_sd = pv_response_sd_pulse_width,
        response_baseline = pv_response_baseline,
        response_halflife = pv_response_halflife,
        response_location = pv_response_pulse_location,
        response_pulse_mass = pv_response_indiv_pulse_mass,
        response_pulse_width = pv_response_indiv_pulse_width,
        response_sdscale_pulse_mass = pv_response_sdscale_pulse_mass,
        response_sdscale_pulse_width = pv_response_sdscale_pulse_width,
        # Association (C++ code expects "rho" and "nu", not "log_rho" and "log_nu")
        rho = pv_log_rho,
        nu = pv_log_nu
      ),

      # Driver starting values
      driver_starting_values = list(
        mass_mean = sv_driver_mass_mean,
        width_mean = sv_driver_width_mean,
        baseline = sv_driver_baseline_mean,
        halflife = sv_driver_halflife_mean,
        mass_sd = sv_driver_mass_sd,
        width_sd = sv_driver_width_sd,
        errorsq = sv_driver_error_var
      ),

      # Response starting values
      response_starting_values = list(
        mass_mean = sv_response_mass_mean,
        width_mean = sv_response_width_mean,
        baseline = sv_response_baseline_mean,
        halflife = sv_response_halflife_mean,
        mass_sd = sv_response_mass_sd,
        width_sd = sv_response_width_sd,
        errorsq = sv_response_error_var
      ),

      # Association starting values
      association_starting_values = list(
        rho = sv_rho,
        nu = sv_nu
      )
    ),
    class = "joint_spec"
  )

  return(spec_obj)
}


#' @export
print.joint_spec <- function(x, ...) {
  cat("\nBayesian Joint Hormone Model Specification\n\n")
  cat("Location prior:", x$location_prior, "\n\n")

  cat("Driver Hormone:\n")
  cat("  Mean pulse count prior:", x$driver_priors$pulse_count, "\n")
  cat("  Pulse mass mean prior:", x$driver_priors$mass_mean,
      "(variance:", x$driver_priors$mass_variance, ")\n")
  cat("  Pulse width mean prior:", x$driver_priors$width_mean,
      "(variance:", x$driver_priors$width_variance, ")\n")
  cat("  Starting values: baseline =", x$driver_starting_values$baseline,
      ", halflife =", x$driver_starting_values$halflife, "\n\n")

  cat("Response Hormone:\n")
  cat("  Mean pulse count prior:", x$response_priors$pulse_count, "\n")
  cat("  Pulse mass mean prior:", x$response_priors$mass_mean,
      "(variance:", x$response_priors$mass_variance, ")\n")
  cat("  Pulse width mean prior:", x$response_priors$width_mean,
      "(variance:", x$response_priors$width_variance, ")\n")
  cat("  Starting values: baseline =", x$response_starting_values$baseline,
      ", halflife =", x$response_starting_values$halflife, "\n\n")

  cat("Association (Driver -> Response):\n")
  cat("  log(rho) prior: mean =", x$association_priors$log_rho_mean,
      ", variance =", x$association_priors$log_rho_var, "\n")
  cat("  log(nu) prior: mean =", x$association_priors$log_nu_mean,
      ", variance =", x$association_priors$log_nu_var, "\n")
  cat("  Starting values: rho =", x$association_starting_values$rho,
      ", nu =", x$association_starting_values$nu, "\n\n")

  invisible(x)
}


#------------------------------------------------------------------------------#
#    End of file                                                               #
#------------------------------------------------------------------------------#
