###############################################################################
## nngp_utils.R
##
## Covariance kernels live in covariance_utils.R; every routine here takes a
## kernel function of (distance, theta) as an argument.
##
## CONVENTIONS
## -----------
## * All routines are location level: length S vectors, or S x . matrices,
##   never observation level.
## * All routines expect their inputs in the NNGP ordering produced by
##   order_coords(). Callers permute in with coord_ordering and back out with
##   coord_reverse_ordering. Nothing here reorders silently.
## * Covariance is always parameterised as C = tau2 * K(theta, d).
##   get_F_and_B() is therefore always called with tau = 1, so B_s and F_s
##   refer to the unit kernel K. Multiply F_s by tau2 wherever the scaled
##   covariance is needed, and divide quadratic forms by tau2. This keeps the
##   tau2 Gibbs update, which needs the quadratic form under K rather than C,
##   free of rescaling.
###############################################################################


# ============================================================================
# 1. DISTANCE / ORDERING / NEIGHBOUR STRUCTURE
# ============================================================================

#' Euclidean distance
#' @noRd
euc_dist <- function(x1, y1, x2, y2) {
    sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

#' Order coordinates from upper-right to lower-left for NNGP
#'
#' Neighbours are always drawn from later entries of this ordering, so the
#' final location in the ordering is the root of the DAG.
#'
#' @param coords Matrix with columns "x" and "y" (S x 2), one row per location
#' @param space_id Location ID vector (S), aligned with rows of coords
#' @return List with ordered_coords, coord_ordering, coord_reverse_ordering
#' @noRd
order_coords <- function(coords, space_id) {
    coords <- as.matrix(coords)
    coords_dist <- euc_dist(coords[, "x"],
                            coords[, "y"],
                            max(coords[, "x"]),
                            max(coords[, "y"]))
    coord_ordering <- order(coords_dist)
    coords_ordered <- coords[coord_ordering, , drop = FALSE]
    space_id_ordered <- space_id[coord_ordering]
    coord_reverse_ordering <- match(space_id, space_id_ordered)

    list(ordered_coords = coords_ordered,
         coord_ordering = coord_ordering,
         coord_reverse_ordering = coord_reverse_ordering)
}

#' Get the m nearest neighbours, among successors in the ordering
#'
#' @param ordered_coords Matrix with columns "x","y" in NNGP order
#' @param m Number of neighbours
#' @noRd
get_neighbors <- function(ordered_coords, m) {
    ordered_coords <- as.matrix(ordered_coords)
    S <- nrow(ordered_coords)
    neighbors <- vector("list", S)
    for (j in 1:(S - 1)) {
        dists <- euc_dist(ordered_coords[j, "x"],
                          ordered_coords[j, "y"],
                          ordered_coords[(j + 1):S, "x"],
                          ordered_coords[(j + 1):S, "y"])
        neighbors[[j]] <- order(dists)[1:min(m, length(dists))] + j
    }
    neighbors[[S]] <- integer(0)
    neighbors
}

#' For each location, which other locations have it as a neighbour
#' @noRd
get_neighbors_inverse <- function(neighbors) {
    S <- length(neighbors)
    neighbors_inverse <- vector("list", S)
    for (i in seq_len(S)) {
        neighbors_inverse[[i]] <- which(vapply(neighbors,
                                               function(x) i %in% x,
                                               logical(1)))
    }
    neighbors_inverse
}

#' Pre-compute the position of each location within each neighbour vector
#'
#' Element s is a length S integer vector: 0 if s is not a neighbour of 
#' that location, otherwise its index within that location's neighbour vector.
#' @noRd
get_neighbor_positions <- function(neighbors) {
    S <- length(neighbors)
    pos_in_neighbors <- vector("list", S)
    for (s in seq_len(S)) {
        pos <- rep(0L, S)
        nb <- neighbors[[s]]
        if (length(nb) > 0) {
            pos[nb] <- seq_along(nb)
        }
        pos_in_neighbors[[s]] <- pos
    }
    pos_in_neighbors
}

#' Distance matrices for each location and its neighbour set
#'
#' Element i is the (1 + number of neighbours) square distance matrix with
#' location i first.
#' @noRd
get_dist_matrices <- function(coords, neighbors) {
    coords <- as.matrix(coords)
    S <- nrow(coords)
    dist_matrices <- vector("list", S)
    for (i in seq_len(S)) {
        sub <- coords[c(i, neighbors[[i]]), , drop = FALSE]
        dist_matrices[[i]] <- as.matrix(stats::dist(sub, upper = TRUE, diag = TRUE))
    }
    dist_matrices
}

#' Nearest reference (training) neighbours for each prediction location
#' @noRd
get_neighbors_ref <- function(ordered_coords, pred_coords, m) {
    ordered_coords <- as.matrix(ordered_coords)
    pred_coords <- as.matrix(pred_coords)
    neighbors_pred <- vector("list", nrow(pred_coords))
    for (i in seq_len(nrow(pred_coords))) {
        dists <- euc_dist(pred_coords[i, "x"],
                          pred_coords[i, "y"],
                          ordered_coords[, "x"],
                          ordered_coords[, "y"])
        neighbors_pred[[i]] <- order(dists)[1:min(m, nrow(ordered_coords))]
    }
    neighbors_pred
}

#' Distance matrices from prediction locations to their reference neighbours
#' @noRd
get_dist_matrices_ref <- function(ordered_coords, coords_pred, neighbors) {
    ordered_coords <- as.matrix(ordered_coords)
    coords_pred <- as.matrix(coords_pred)
    dist_matrices <- vector("list", nrow(coords_pred))
    for (i in seq_len(nrow(coords_pred))) {
        sub <- rbind(coords_pred[i, , drop = FALSE],
                     ordered_coords[neighbors[[i]], , drop = FALSE])
        dist_matrices[[i]] <- as.matrix(stats::dist(sub, upper = TRUE, diag = TRUE))
    }
    dist_matrices
}

#' Build the full NNGP structure for a set of unique locations
#'
#' Convenience wrapper so pg_ensemble() has a single setup call.
#'
#' @param coords Matrix with columns "x","y" (S x 2), one row per location
#' @param space_id Location IDs (S), aligned with rows of coords
#' @param m Number of neighbours
#' @noRd
nngp_setup <- function(coords, space_id, m) {
    coords <- as.matrix(coords)
    S <- nrow(coords)
    if (m < 1) {
        stop("'number_neighbors' must be at least 1.")
    }
    if (m >= S) {
        stop("'number_neighbors' must be smaller than the number of unique ",
             "spatial locations (", S, "). Use nngp = FALSE for small S.")
    }

    ord <- order_coords(coords = coords, space_id = space_id)
    neighbors <- get_neighbors(ord$ordered_coords, m = m)

    list(ordered_coords = ord$ordered_coords,
         coord_ordering = ord$coord_ordering,
         coord_reverse_ordering = ord$coord_reverse_ordering,
         neighbors = neighbors,
         neighbors_inverse = get_neighbors_inverse(neighbors),
         pos_in_neighbors = get_neighbor_positions(neighbors),
         dist_matrices = get_dist_matrices(ord$ordered_coords, neighbors),
         number_neighbors = m)
}


# ============================================================================
# 2. NNGP CONDITIONAL PARAMETERS (B, F)
# ============================================================================

#' Evaluate a kernel on every neighbour-set distance matrix
#'
#' Returns both the kernels and the inverses of the neighbour-only blocks,
#' which are what get_F_and_B() consumes. The last location has no neighbours,
#' so only the first S - 1 partial inverses are formed.
#'
#' @param dist_matrices List from get_dist_matrices()
#' @param theta Range parameter
#' @param cov_kern Kernel function of (distance, theta)
#' @noRd
get_nngp_kernels <- function(dist_matrices, theta, cov_kern) {
    kernels <- lapply(dist_matrices,
                      function(d) cov_kern(distance = d, theta = theta))
    kernels_partial_inv <- lapply(
        kernels[seq_len(length(kernels) - 1)],
        function(k) solve(k[-1, -1, drop = FALSE]))
    list(kernels = kernels, kernels_partial_inv = kernels_partial_inv)
}

#' NNGP conditional regression weights B_s and conditional variances F_s
#'
#' B_s = cross-cov %*% solve(neighbour-cov); F_s = conditional variance.
#' Call with tau = 1 to obtain unit-kernel quantities, per the file header.
#' @noRd
get_F_and_B <- function(kernels, tau, neighbors, kernels_partial_inverse) {
    S <- length(kernels)
    F_s <- numeric(S)
    B_s <- vector("list", S)

    for (s in seq_len(S)) {
        if (length(neighbors[[s]]) > 0) {
            B_s[[s]] <- as.numeric(kernels[[s]][1, -1, drop = FALSE] %*%
                                       kernels_partial_inverse[[s]])
            F_s[s] <- tau * (kernels[[s]][1, 1] - B_s[[s]] %*% kernels[[s]][-1, 1])
        } else {
            B_s[[s]] <- numeric(0)
            F_s[s] <- tau * kernels[[s]][1, 1]
        }
    }
    list(B = B_s, F = F_s)
}

#' One-call helper: distance matrices to (B, F) for a given range parameter
#' @noRd
get_nngp_B_and_F <- function(dist_matrices, theta, cov_kern, neighbors, tau = 1) {
    kerns <- get_nngp_kernels(dist_matrices, theta, cov_kern)
    get_F_and_B(kernels = kerns$kernels,
                tau = tau,
                neighbors = neighbors,
                kernels_partial_inverse = kerns$kernels_partial_inv)
}


# ============================================================================
# 3. SPARSE PRECISION ALGEBRA: C^{-1} = (I - B)' F^{-1} (I - B)
# ============================================================================

#' (I - B) %*% v, for v in NNGP order
#' @noRd
nngp_ImB <- function(v, B_s, neighbors) {
    out <- as.numeric(v)
    for (s in seq_along(out)) {
        nb <- neighbors[[s]]
        if (length(nb) > 0) {
            out[s] <- out[s] - sum(B_s[[s]] * v[nb])
        }
    }
    out
}

#' t(I - B) %*% v, for v in NNGP order
#' @noRd
nngp_ImBt <- function(v, B_s, neighbors) {
    out <- as.numeric(v)
    for (s in seq_along(v)) {
        nb <- neighbors[[s]]
        if (length(nb) > 0) {
            out[nb] <- out[nb] - B_s[[s]] * v[s]
        }
    }
    out
}

#' Quadratic form r' C^{-1} r under the NNGP
#'
#' Pass unit-kernel F_s to get the quadratic form under K(theta); divide the
#' result by tau2 for the form under C = tau2 * K(theta).
#' @noRd
nngp_quadform <- function(resid, B_s, F_s, neighbors) {
    r_tilde <- nngp_ImB(resid, B_s, neighbors)
    sum(r_tilde^2 / F_s)
}

#' C^{-1} %*% v under the NNGP, for v in NNGP order
#' @noRd
nngp_precision_multiply <- function(v, B_s, F_s, neighbors) {
    nngp_ImBt(nngp_ImB(v, B_s, neighbors) / F_s, B_s, neighbors)
}

#' NNGP log density of a mean-zero residual vector under C = tau * K
#'
#' F_s must be the unit-kernel conditional variances. Replaces looping dnngp()
#' calls: one pass, no per-site solve.
#' @noRd
nngp_loglik <- function(resid, B_s, F_s, neighbors, tau = 1) {
    S <- length(resid)
    quad <- nngp_quadform(resid, B_s, F_s, neighbors)
    -0.5 * (S * log(2 * pi) + S * log(tau) + sum(log(F_s)) + quad / tau)
}


# ============================================================================
# 4. SEQUENTIAL GIBBS UPDATE FOR AN NNGP SPATIAL FIELD
# ============================================================================

#' Sequential (site-by-site) Gibbs update of an NNGP field
#'
#' Generic scan covering both spatial fields in the ensemble model. The
#' likelihood contribution at site s must be expressible, up to a constant, as
#'
#'     exp( lik_lin[s] * u[s] - 0.5 * lik_prec[s] * u[s]^2 )
#'
#' and the prior is NNGP with mean mu:
#'
#'     u_s | u_nb ~ N( mu_s + B_s' (u_nb - mu_nb), F_s )
#'
#' Collecting the site's own prior term, the terms from every location that
#' has s as a neighbour (its "children"), and the likelihood gives
#'
#'     prec_s = lik_prec[s] + 1/F_s + sum_tt B_tts^2 / F_tt
#'     m_s = prec_s^{-1} * ( lik_lin[s]
#'                           + (1/F_s) * (mu_s + B_s'(u_nb - mu_nb))
#'                           + sum_tt B_tts r_tt / F_tt
#'                           + mu_s * sum_tt B_tts^2 / F_tt )
#'
#' where r_tt = (u_tt - mu_tt) - sum_{l in nb(tt), l != s} B_ttl (u_l - mu_l).
#' The final term is what the zero-mean GRM version does not need.
#'
#' Usage:
#'   psi_k: lik_prec = omega_k (Polya-Gamma), lik_lin = kappa_k, mu = X beta_k
#'   delta: lik_prec = b_s, lik_lin = e_s, mu = 0
#'
#' All arguments in NNGP order. F_s must be the scaled conditional variances
#' (tau2 * unit F), since the prior here is on the actual field.
#'
#' NOTE: this mirrors the argument order of the GRM's
#' mcmc_draw_spatial_nngp_cpp() so it can be swapped for an Rcpp version
#' without touching the callers.
#' @noRd
mcmc_draw_nngp_effect <- function(effect,
                                  mu,
                                  lik_prec,
                                  lik_lin,
                                  neighbors,
                                  neighbors_inverse,
                                  pos_in_neighbors,
                                  B_s,
                                  F_s) {

    S <- length(effect)

    for (s in seq_len(S)) {

        # contributions from locations that have s as a neighbour
        sum_B_F_inv_B <- 0
        sum_B_F_inv_a <- 0

        for (tt in neighbors_inverse[[s]]) {
            B_tt <- B_s[[tt]]
            B_tts <- B_tt[pos_in_neighbors[[tt]][s]]
            F_tt <- F_s[tt]

            r_tt <- effect[tt] - mu[tt]
            for (l in neighbors[[tt]]) {
                if (l != s) {
                    r_tt <- r_tt - B_tt[pos_in_neighbors[[tt]][l]] * (effect[l] - mu[l])
                }
            }

            sum_B_F_inv_B <- sum_B_F_inv_B + B_tts^2 / F_tt
            sum_B_F_inv_a <- sum_B_F_inv_a + B_tts * r_tt / F_tt
        }

        # the site's own NNGP conditional mean given its neighbours
        parent_mean <- mu[s]
        if (length(neighbors[[s]]) > 0) {
            nb <- neighbors[[s]]
            parent_mean <- parent_mean + as.numeric(B_s[[s]] %*% (effect[nb] - mu[nb]))
        }

        V_s <- 1 / (lik_prec[s] + 1 / F_s[s] + sum_B_F_inv_B)
        m_s <- lik_lin[s] +
            (1 / F_s[s]) * parent_mean +
            sum_B_F_inv_a +
            sum_B_F_inv_B * mu[s]

        effect[s] <- stats::rnorm(1, V_s * m_s, sqrt(V_s))
    }

    effect
}

#' Sequential Gibbs update of psi_k, Polya-Gamma likelihood with NNGP prior
#' @noRd
mcmc_draw_psi_nngp <- function(psi_k, mu_k, kappa_k, omega_k,
                               neighbors, neighbors_inverse, pos_in_neighbors,
                               B_s, F_s) {
    mcmc_draw_nngp_effect(effect = psi_k,
                          mu = mu_k,
                          lik_prec = omega_k,
                          lik_lin = kappa_k,
                          neighbors = neighbors,
                          neighbors_inverse = neighbors_inverse,
                          pos_in_neighbors = pos_in_neighbors,
                          B_s = B_s,
                          F_s = F_s)
}

#' Sequential Gibbs update of the shared spatial intercept delta
#'
#' @param delta Current delta vector in NNGP order
#' @param b_s Per-location precision sum, sum_t 1 / v_{s,t}
#' @param e_s Per-location weighted residual sum, sum_t r_{s,t} / v_{s,t}
#' @noRd
mcmc_draw_delta_nngp <- function(delta, b_s, e_s,
                                 neighbors, neighbors_inverse, pos_in_neighbors,
                                 B_s, F_s) {
    mcmc_draw_nngp_effect(effect = delta,
                          mu = rep(0, length(delta)),
                          lik_prec = b_s,
                          lik_lin = e_s,
                          neighbors = neighbors,
                          neighbors_inverse = neighbors_inverse,
                          pos_in_neighbors = pos_in_neighbors,
                          B_s = B_s,
                          F_s = F_s)
}


# ============================================================================
# 5. NNGP KRIGING (PREDICTION)
# ============================================================================

#' Krige one or more NNGP fields to new locations
#'
#' Vectorised over fields (columns) so the neighbour solve is done once per
#' prediction location per MCMC draw, rather than once per stick-breaking
#' component.
#'
#' @param field_ord S_train x J matrix of field values, in NNGP order
#' @param mu_ord S_train x J matrix of field means at training locations, in
#'   NNGP order
#' @param mu_pred S_pred x J matrix of field means at prediction locations
#' @param tau Variance scale (tau2)
#' @param theta Range parameter
#' @param neighbors_pred List from get_neighbors_ref()
#' @param dist_matrices_pred List from get_dist_matrices_ref()
#' @param cov_kern Kernel function of (distance, theta)
#' @return List with mean (S_pred x J) and var (length S_pred; the conditional
#'   variance does not depend on the field)
#' @noRd
nngp_krige <- function(field_ord, mu_ord, mu_pred, tau, theta,
                       neighbors_pred, dist_matrices_pred, cov_kern) {

    field_ord <- as.matrix(field_ord)
    mu_ord <- as.matrix(mu_ord)
    mu_pred <- as.matrix(mu_pred)

    n_pred <- length(neighbors_pred)
    J <- ncol(field_ord)

    resid_ord <- field_ord - mu_ord

    mean_out <- matrix(0, nrow = n_pred, ncol = J)
    var_out <- numeric(n_pred)

    for (i in seq_len(n_pred)) {
        nb <- neighbors_pred[[i]]
        S_i <- tau * cov_kern(distance = dist_matrices_pred[[i]], theta = theta)
        cross <- S_i[1, -1, drop = FALSE]
        w_i <- cross %*% solve(S_i[-1, -1, drop = FALSE])

        mean_out[i, ] <- mu_pred[i, ] + as.numeric(w_i %*% resid_ord[nb, , drop = FALSE])
        var_out[i] <- max(S_i[1, 1] - as.numeric(w_i %*% t(cross)), 0)
    }

    list(mean = mean_out, var = var_out)
}
