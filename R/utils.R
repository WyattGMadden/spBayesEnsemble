
logit <- function(x) log(x / (1 - x))

ilogit <- function(x) 1 / (1 + exp(-x))


#' Polya-Gamma draws that tolerate zero counts
#'
#' PG(0, z) is degenerate at 0, but BayesLogit::rpg.devroye() does not accept
#' h = 0. Zero counts arise whenever a location has no observations left for a
#' stick-breaking component, which is common once the weights concentrate.
#'
#' @param h Vector of PG shape parameters (counts)
#' @param z Vector of PG tilting parameters
#' @noRd
rpg_safe <- function(h, z) {
    out <- numeric(length(h))
    pos <- h > 0
    if (any(pos)) {
        out[pos] <- BayesLogit::rpg.devroye(num = sum(pos), h = h[pos], z = z[pos])
    }
    out
}

#' Pull the (M - 1) x S psi slice for one MCMC draw as an S x (M - 1) matrix
#'
#' The explicit matrix() guards the M = 2 case, where the array slice drops
#' to a vector.
#'
#' @param psi Array of dimension (M - 1) x S x n_iter
#' @param t Iteration index
#' @param Km1 Number of stick-breaking components, M - 1
#' @noRd
psi_slice <- function(psi, t, Km1) {
    t(matrix(psi[, , t], nrow = Km1))
}

#' Pull the (M - 1) x P beta slice for one MCMC draw as a P x (M - 1) matrix
#'
#' @param betas Array of dimension (M - 1) x P x n_iter
#' @param t Iteration index
#' @param Km1 Number of stick-breaking components, M - 1
#' @noRd
beta_slice <- function(betas, t, Km1) {
    t(matrix(betas[, , t], nrow = Km1))
}
