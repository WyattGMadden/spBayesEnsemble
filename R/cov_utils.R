###############################################################################
## covariance_utils.R
##
## Covariance kernels and kernel resolution. Kept separate from nngp_utils.R:
## nothing here is NNGP specific, and the NNGP routines take a kernel function
## as an argument rather than depending on this file.
###############################################################################

#' Exponential covariance kernel
#'
#' @param distance Distance, scalar / vector / matrix
#' @param theta Range parameter
#' @noRd
exponential_kernel <- function(distance, theta) {
    exp(-distance / theta)
}

#' Matern covariance kernel
#'
#' @param distance Distance, scalar / vector / matrix
#' @param theta Range parameter
#' @param nu Smoothness, one of 0.5, 1.5, 2.5
#' @noRd
matern_kernel <- function(distance, theta, nu = 1.5) {
    d <- distance / theta
    if (nu == 0.5) {
        exp(-d)
    } else if (nu == 1.5) {
        (1 + sqrt(3) * d) * exp(-sqrt(3) * d)
    } else if (nu == 2.5) {
        (1 + sqrt(5) * d + 5 / 3 * d^2) * exp(-sqrt(5) * d)
    } else {
        stop("'matern_nu' must be 0.5, 1.5, or 2.5.")
    }
}

#' Resolve a covariance kernel function from its name
#'
#' Returns a function of (distance, theta), which is the interface every
#' downstream routine assumes.
#'
#' @param covariance One of "exponential", "matern", "custom"
#' @param matern_nu Smoothness used when covariance = "matern"
#' @param covariance_kernel User supplied kernel used when covariance =
#'   "custom". Must be a function with "distance" and "theta" arguments.
#' @noRd
get_cov_kern <- function(covariance = "exponential",
                         matern_nu = 1.5,
                         covariance_kernel = NULL) {
    if (covariance == "exponential") {
        exponential_kernel
    } else if (covariance == "matern") {
        nu_val <- as.numeric(matern_nu)
        if (!(nu_val %in% c(0.5, 1.5, 2.5))) {
            stop("'matern_nu' must be 0.5, 1.5, or 2.5.")
        }
        function(distance, theta) matern_kernel(distance, theta, nu = nu_val)
    } else if (covariance == "custom") {
        if (!is.function(covariance_kernel)) {
            stop("If covariance = 'custom', 'covariance_kernel' must be a ",
                 "function with 'distance' and 'theta' arguments.")
        }
        covariance_kernel
    } else {
        stop("'covariance' must be one of 'exponential', 'matern', 'custom'.")
    }
}
