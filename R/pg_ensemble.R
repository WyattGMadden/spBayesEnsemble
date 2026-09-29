#' Fit spatial Bayesian ensemble
#'
#' This function fits a spatial Bayesian ensemble using model predictive
#' means and standard deviations, and (optional) ensemble weight covariates.
#' Component densities are computed internally as N(y | model_est, model_sd^2).
#'
#' @param y Vector of observations (N)
#' @param model_est Matrix or data frame of model predictive means (N x M).
#'   Column names, if provided, are used as model names; otherwise names
#'   "model_1", ..., "model_M" are assigned.
#' @param model_sd Matrix or data frame of model predictive standard
#'   deviations (N x M). Columns are paired with model_est positionally
#'   (column j of model_sd belongs to column j of model_est); its column
#'   names, if any, are ignored.
#' @param X Matrix of covariates for the ensemble weight GP mean (S x P), one
#'   row per unique spatial location, ordered by space_id
#' @param space_id Spatial location ID vector (N). Must be integers
#'   1, ..., S with no gaps.
#' @param coords Matrix of x y coordinates, with colnames(coords) == c("x", "y"), (N, 2)
#' @param intercept Include shared (across member models) spatial intercept to means
#' @param nngp Use a nearest neighbor Gaussian process (NNGP) in place of a
#'   full Gaussian process for the stick-breaking logits and, when requested,
#'   the shared spatial intercept
#' @param number_neighbors Number of nearest neighbors to use in the NNGP
#' @param covariance Covariance function, one of "exponential", "matern",
#'   "custom"
#' @param covariance_kernel Custom covariance function, used when
#'   covariance = "custom". Must be a function with "distance" and "theta"
#'   arguments.
#' @param matern_nu Nu parameter for the Matern covariance (0.5, 1.5, or 2.5)
#' @param n_iter Number of iterations used in MCMC
#' @param beta_prior_var Variance of normal prior placed on betas
#' @param model_names Optional character vector of model names (length M),
#'   applied positionally to the columns of model_est/model_sd and used to
#'   name weight output columns (e.g. "cmaqgrm" yields "cmaqgrm_weight" in
#'   pg_pred() output). Overrides colnames(model_est) if both are present.
#'   If NULL, colnames(model_est) are used, or "model_1", ..., "model_M"
#'   if model_est is unnamed.
#' @param verbose Print MCMC progress
#' @param verbose_iter Print progress every verbose_iter iterations
#'
#' @return A list containing MCMC output
#'
#' @examples
#' # pg_ensemble()
#'
#'
#' @export
pg_ensemble <- function(
    y, model_est, model_sd, X, space_id, coords, 
    intercept = FALSE,
    nngp = FALSE,
    number_neighbors = 10,
    covariance = "exponential",
    covariance_kernel = NULL,
    matern_nu = 1.5,
    n_iter = 1000, 
    beta_prior_var = 100, 
    model_names = NULL,
    verbose = TRUE,
    verbose_iter = 100) {

    if (!is.null(model_names)) {
        model_est <- as.matrix(model_est)
        if (length(model_names) != ncol(model_est)) {
            stop("model_names must have one name per column of model_est.")
        }
        colnames(model_est) <- model_names
    }
    aligned <- align_est_sd(model_est, model_sd)
    model_est <- aligned$model_est
    model_sd <- aligned$model_sd
    model_names <- aligned$model_names

    if (nrow(model_est) != length(y) || length(space_id) != length(y)) {
        stop("y, space_id, and the rows of model_est/model_sd must all have length N.")
    }

    cov_kern <- get_cov_kern(covariance = covariance,
                             matern_nu = matern_nu,
                             covariance_kernel = covariance_kernel)

    dens <- stats::dnorm(y, mean = model_est, sd = model_sd) |>
        matrix(nrow = length(y), ncol = ncol(model_est)) |>
        t()
    locs <- cbind(space_id = space_id, coords) |>
        unique() |>
        as.data.frame()
    locs <- locs[order(locs$space_id), ]
    X <- as.matrix(X)
    P <- ncol(X)
    S <- nrow(locs)
    M <- nrow(dens)

    if (!all(sort(unique(space_id)) == seq_len(S))) {
        stop("space_id must be integers 1, ..., S with no gaps.")
    }
    if (nrow(X) != S) {
        stop("X must have one row per unique spatial location (S = ", S,
             "), ordered by space_id.")
    }

    #######################################
    ######## spatial structure ############
    #######################################

    distmat <- NULL
    nngp_info <- NULL

    if (nngp) {
        if (verbose) {
            cat("Setting up NNGP structures (m =", number_neighbors, "neighbors)\n")
        }
        nngp_info <- nngp_setup(coords = as.matrix(locs[, c("x", "y")]),
                                space_id = locs$space_id,
                                m = number_neighbors)
        neighbors <- nngp_info$neighbors
        neighbors_inverse <- nngp_info$neighbors_inverse
        pos_in_neighbors <- nngp_info$pos_in_neighbors
        dist_matrices <- nngp_info$dist_matrices
        ord <- nngp_info$coord_ordering
        rord <- nngp_info$coord_reverse_ordering
        X_ord <- X[ord, , drop = FALSE]
    } else {
        distmat <- stats::dist(locs[, c("x", "y")]) |>
            as.matrix()
    }

    weights_all <- array(1, dim = c(M, S, n_iter))

    #weights_all[, , 1] <- counts_by_spat[, c("weight1", "weight2", "weight3")] |>
    #    as.matrix()
    weights_all[, , 1] <- matrix(1/M, nrow = M, ncol = S)

    #psi1 <- matrix(0, nrow = S, ncol = n_iter)
    #psi2 <- matrix(0, nrow = S, ncol = n_iter)
    psi <- array(0, dim = c(M - 1, S, n_iter))

    betas <- array(0, dim = c(M - 1, P, n_iter))
    betas[, , 1] <- 0

    # init gp params
    tau2 <- rep(0, n_iter)
    tau2[1] <- 1
    tau2_a <- 0.001
    tau2_b <- 0.001

    rho <- rep(0, n_iter)
    rho[1] <- 20
    rho_step <- 2
    rho_mu <- 3
    rho_sd <- 1

    # init ensemble intercept
    delta <- matrix(0, nrow = S, ncol = n_iter)
    tau2_delta <- rep(0, n_iter)
    tau2_delta[1] <- 1
    tau2_delta_a <- 0.001
    tau2_delta_b <- 0.001
    rho_delta <- rep(20, n_iter)
    rho_delta_step <- 2
    rho_delta_mu <- 3
    rho_delta_sd <- 1

    for (i in 2:n_iter) {

        #recompute densities with current spatial intercept
        if (intercept) {
            delta_obs <- delta[space_id, i - 1]
            dens <- stats::dnorm(y, mean = model_est + delta_obs, sd = model_sd) |>
                matrix(nrow = length(y), ncol = ncol(model_est)) |>
                t()
        }

        #calculate Z
        obs_weights <- weights_all[, , i - 1][, space_id]
        weight_dens <- dens * obs_weights
        zprobs <- sweep(weight_dens, 2, colSums(weight_dens), FUN = "/")

        




        #######################################
        ########gen latent variable z##########
        #######################################

        #vectorized rmultim, over prob matrix
        z_draw <- apply(zprobs, 
                        2, 
                        function(x) stats::rmultinom(n = 1, size = 1, prob = x))

        #counts of z wrt space_id
        z_counts <- apply(z_draw, 1, function(x) tapply(x, space_id, sum))
        
        #cumsums of z counts
        z_counts_rev_cum <- z_counts
        for (j in (ncol(z_counts) - 1):1) {
            z_counts_rev_cum[, j] <- z_counts_rev_cum[, j + 1] + z_counts_rev_cum[, j]
        }
        z_counts_rev_cum <- z_counts_rev_cum[, 1:(ncol(z_counts_rev_cum) - 1), drop = FALSE]

         
        #######################################
        ######## spatial intercept (delta) ####
        #######################################
        if (intercept) {
            #assigned-model residual and variance per observation
            est_assigned <- colSums(z_draw * t(model_est))
            var_assigned <- colSums(z_draw * t(model_sd^2))
            resid <- y - est_assigned

            b_s <- as.numeric(tapply(1 / var_assigned, space_id, sum))
            e_s <- as.numeric(tapply(resid / var_assigned, space_id, sum))

            if (!nngp) {

                K_delta <- cov_kern(distance = distmat, theta = rho_delta[i - 1])
                K_delta_inv <- solve(K_delta)

                V_delta <- solve(diag(b_s) + K_delta_inv / tau2_delta[i - 1])
                m_delta <- V_delta %*% e_s
                delta[, i] <- t(mvtnorm::rmvnorm(n = 1, mean = m_delta, sigma = V_delta))

                #quadratic form under the unit kernel
                SSS_delta <- t(delta[, i]) %*% K_delta_inv %*% delta[, i]

            } else {

                #unit-kernel B and F for the current range
                BF_d <- get_nngp_B_and_F(dist_matrices = dist_matrices,
                                         theta = rho_delta[i - 1],
                                         cov_kern = cov_kern,
                                         neighbors = neighbors,
                                         tau = 1)

                delta_ord <- mcmc_draw_delta_nngp(
                    delta = delta[ord, i - 1],
                    b_s = b_s[ord],
                    e_s = e_s[ord],
                    neighbors = neighbors,
                    neighbors_inverse = neighbors_inverse,
                    pos_in_neighbors = pos_in_neighbors,
                    B_s = BF_d$B,
                    F_s = tau2_delta[i - 1] * BF_d$F
                    )
                delta[, i] <- delta_ord[rord]

                #quadratic form under the unit kernel
                SSS_delta <- nngp_quadform(delta_ord, BF_d$B, BF_d$F, neighbors)
            }

            #update tau2_delta
            SSS_delta <- SSS_delta / 2
            tau2_delta[i] <- 1 / stats::rgamma(1, S / 2 + tau2_delta_a, SSS_delta + tau2_delta_b)

            #update rho_delta (MH)
            rho_delta_prop <- stats::rlnorm(1, log(rho_delta[i - 1]), rho_delta_step)
            if (!nngp) {

                C_delta_curr <- tau2_delta[i] * K_delta
                C_delta_prop <- tau2_delta[i] * cov_kern(distance = distmat,
                                                         theta = rho_delta_prop)
                lik_curr_d <- mvtnorm::dmvnorm(delta[, i], rep(0, S), C_delta_curr, log = TRUE)
                lik_prop_d <- mvtnorm::dmvnorm(delta[, i], rep(0, S), C_delta_prop, log = TRUE)

            } else {

                BF_d_prop <- get_nngp_B_and_F(dist_matrices = dist_matrices,
                                              theta = rho_delta_prop,
                                              cov_kern = cov_kern,
                                              neighbors = neighbors,
                                              tau = 1)
                lik_curr_d <- nngp_loglik(delta_ord, BF_d$B, BF_d$F,
                                          neighbors, tau = tau2_delta[i])
                lik_prop_d <- nngp_loglik(delta_ord, BF_d_prop$B, BF_d_prop$F,
                                          neighbors, tau = tau2_delta[i])
            }

            ratio_d <- lik_prop_d +
                stats::dlnorm(rho_delta_prop, rho_delta_mu, rho_delta_sd, log = TRUE) +
                log(rho_delta_prop) -
                lik_curr_d -
                stats::dlnorm(rho_delta[i - 1], rho_delta_mu, rho_delta_sd, log = TRUE) -
                log(rho_delta[i - 1])

            if (log(stats::runif(1)) < ratio_d) {
                rho_delta[i] <- rho_delta_prop
            } else {
                rho_delta[i] <- rho_delta[i - 1]
            }
        }


        #######################################
        ########weights sample##########
        #######################################

        #calculate gp cov (unit kernel; tau2 applied where needed)
        if (!nngp) {
            K_unit <- cov_kern(distance = distmat, theta = rho[i - 1])
            K_unit_inv <- solve(K_unit)
            covar_inv <- K_unit_inv / tau2[i - 1]
        } else {
            BF <- get_nngp_B_and_F(dist_matrices = dist_matrices,
                                   theta = rho[i - 1],
                                   cov_kern = cov_kern,
                                   neighbors = neighbors,
                                   tau = 1)
            F_scaled <- tau2[i - 1] * BF$F
        }

        #draw pg variables
        pg_vars <- matrix(0, nrow = M - 1, ncol = S)
        for (j in 1:(M - 1)) {
            pg_vars[j, ] <- rpg_safe(h = z_counts_rev_cum[, j],
                                     z = psi[j, , i - 1])
        }
        
        k_z <- matrix(0, nrow = M - 1, ncol = S)
        for (j in 1:(M - 1)) {
            k_z[j, ] <- z_counts[, j] - (z_counts_rev_cum[, j]) / 2
        }


        MMM <- matrix(0, nrow = M - 1, ncol = S)
        for (j in 1:(M - 1)) {
            mu_j <- as.numeric(X %*% betas[j, , i - 1])
            if (!nngp) {
                omega <- diag(pg_vars[j, ])
                Sigma <- solve(omega + covar_inv)
                MMM[j, ] <- Sigma %*% (k_z[j, ] + covar_inv %*% mu_j)
                psi[j, , i] <- t(mvtnorm::rmvnorm(n = 1, mean = MMM[j, ], sigma = Sigma))
            } else {
                psi_ord <- mcmc_draw_psi_nngp(
                    psi_k = psi[j, ord, i - 1],
                    mu_k = mu_j[ord],
                    kappa_k = k_z[j, ord],
                    omega_k = pg_vars[j, ord],
                    neighbors = neighbors,
                    neighbors_inverse = neighbors_inverse,
                    pos_in_neighbors = pos_in_neighbors,
                    B_s = BF$B,
                    F_s = F_scaled
                    )
                psi[j, , i] <- psi_ord[rord]
            }
        }


        #update betas
        #X' C^{-1} X does not depend on j, so form it once
        if (!nngp) {
            XtCinvX <- t(X) %*% covar_inv %*% X
        } else {
            Cinv_X <- apply(X_ord, 2, function(v) {
                nngp_precision_multiply(v, BF$B, F_scaled, neighbors)
            })
            Cinv_X <- matrix(Cinv_X, nrow = S)
            XtCinvX <- crossprod(X_ord, Cinv_X)
        }
        var_beta <- solve(XtCinvX + diag(1 / beta_prior_var, P))

        for (j in 1:(M - 1)) {
            if (!nngp) {
                XtCinv_psi <- t(X) %*% (covar_inv %*% psi[j, , i])
            } else {
                Cinv_psi <- nngp_precision_multiply(psi[j, ord, i], BF$B,
                                                    F_scaled, neighbors)
                XtCinv_psi <- crossprod(X_ord, Cinv_psi)
            }
            mean_beta <- var_beta %*% XtCinv_psi
            betas[j, , i] <- t(mvtnorm::rmvnorm(n = 1, mean = mean_beta, sigma = var_beta))

        }

        
        #update tau
        SSS <- 0
        for (j in 1:(M - 1)) {
            resid_j <- psi[j, , i] - as.numeric(X %*% betas[j, , i])
            if (!nngp) {
                SSS_j <- t(resid_j) %*% K_unit_inv %*% resid_j
            } else {
                SSS_j <- nngp_quadform(resid_j[ord], BF$B, BF$F, neighbors)
            }
            SSS <- SSS + SSS_j
        }
        SSS <- SSS / 2


        tau2[i] <- 1 / stats::rgamma(1, (S * (M - 1)) / 2 + tau2_a, SSS + tau2_b)

        #update rho
        rho_prop <- stats::rlnorm(1, log(rho[i - 1]), rho_step)

        if (!nngp) {
            SSS_curr <- tau2[i] * K_unit
            SSS_prop <- tau2[i] * cov_kern(distance = distmat, theta = rho_prop)
        } else {
            BF_prop <- get_nngp_B_and_F(dist_matrices = dist_matrices,
                                        theta = rho_prop,
                                        cov_kern = cov_kern,
                                        neighbors = neighbors,
                                        tau = 1)
        }

        lik_curr <- 0
        lik_prop <- 0
        for (j in 1:(M - 1)) {
            resid_j <- psi[j, , i] - as.numeric(X %*% betas[j, , i])
            if (!nngp) {
                lik_curr <- lik_curr + mvtnorm::dmvnorm(resid_j, rep(0, S), SSS_curr, log = TRUE)
                lik_prop <- lik_prop + mvtnorm::dmvnorm(resid_j, rep(0, S), SSS_prop, log = TRUE)
            } else {
                lik_curr <- lik_curr + nngp_loglik(resid_j[ord], BF$B, BF$F,
                                                   neighbors, tau = tau2[i])
                lik_prop <- lik_prop + nngp_loglik(resid_j[ord], BF_prop$B, BF_prop$F,
                                                   neighbors, tau = tau2[i])
            }
        }


        ratio <- lik_prop + 
            stats::dlnorm(rho_prop, rho_mu, rho_sd, log = TRUE) + 
            log(rho_prop) -
            lik_curr - 
            stats::dlnorm(rho[i - 1], rho_mu, rho_sd, log = TRUE) - 
            log(rho[i - 1])

        if (log(stats::runif(1)) < ratio) {
            rho[i] <- rho_prop
        } else {
            rho[i] <- rho[i - 1]
        }

        #update weights

        #inverse logit function
        logit_psi <- matrix(ilogit(psi[, , i]), nrow = M - 1)
        for (j in 1:M) {
            if (j > 1) {
                for (k in (1:(j - 1))) {
                    weights_all[j, , i] <- weights_all[j, , i] * (1 - logit_psi[k, ])
                }
            }
            if (j < M) {
                weights_all[j, , i] <- weights_all[j, , i] * logit_psi[j, ]
            }
        }


        if (verbose && (i %% verbose_iter == 0)) {
            cat(paste("     Iteration", i, "of", n_iter, "\n"))
        }

    }

    return(
        list(
            weights_all = weights_all, 
            psi = psi, 
            betas = betas, 
            tau2 = tau2, 
            rho = rho,
            intercept = intercept,
            delta = delta,
            tau2_delta = tau2_delta,
            rho_delta = rho_delta,
            locs = locs,
            distmat = distmat,
            X = X,
            model_names = model_names,
            nngp = nngp,
            nngp_info = nngp_info,
            cov_kern = cov_kern
        )
    )
}


#' Validate and align model estimate and sd matrices (internal)
#'
#' Coerces model_est and model_sd to matrices, checks dimension agreement,
#' and resolves model names. Columns of model_sd are paired with model_est
#' positionally; column names on model_sd are ignored. If model_names is
#' supplied (prediction time) and colnames(model_est) contain all of
#' model_names, both matrices are permuted together to model_names order;
#' if the names share nothing with model_names (e.g. raw data column
#' names), columns are taken positionally; a partial overlap is an error.
#'
#' @param model_est Matrix or data frame of model predictive means (N x M)
#' @param model_sd Matrix or data frame of model predictive sds (N x M)
#' @param model_names Optional character vector of model names to align
#'   columns to (used at prediction time). If NULL, names are taken from
#'   colnames(model_est), or defaults "model_1", ..., "model_M" are assigned.
#'
#' @return A list with elements model_est, model_sd, and model_names, with
#'   columns of both matrices in model_names order.
#'
#' @noRd
align_est_sd <- function(model_est, model_sd, model_names = NULL) {
    model_est <- as.matrix(model_est)
    model_sd <- as.matrix(model_sd)

    if (!identical(dim(model_est), dim(model_sd))) {
        stop("model_est and model_sd must have identical dimensions.")
    }

    est_named <- !is.null(colnames(model_est))

    if (est_named && anyDuplicated(colnames(model_est))) {
        stop("model_est column names must be unique.")
    }

    if (is.null(model_names)) {
        model_names <- if (est_named) {
            colnames(model_est)
        } else {
            paste0("model_", seq_len(ncol(model_est)))
        }
    } else {
        if (ncol(model_est) != length(model_names)) {
            stop("model_est and model_sd must each have one column per model (",
                 length(model_names), " models in pg_fit).")
        }
        if (est_named) {
            n_matched <- sum(model_names %in% colnames(model_est))
            if (n_matched == length(model_names)) {
                # permute est and sd together to preserve positional pairing
                perm <- match(model_names, colnames(model_est))
                model_est <- model_est[, perm, drop = FALSE]
                model_sd <- model_sd[, perm, drop = FALSE]
            } else if (n_matched > 0) {
                stop("colnames(model_est) partially match pg_fit$model_names. ",
                     "Provide all matching names (any order) or non-matching ",
                     "names/no names (columns then taken in pg_fit$model_names order).")
            }
            # n_matched == 0: raw column names (e.g. "*_model_estimate");
            # columns are assumed to be in pg_fit$model_names order
        }
    }

    colnames(model_est) <- model_names
    colnames(model_sd) <- model_names

    list(model_est = model_est, model_sd = model_sd, model_names = model_names)
}
