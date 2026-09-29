#' Predict spatially-varying ensemble weights at new locations
#'
#' Kriges the latent stick-breaking GPs (psi) from a fitted ensemble object
#' to new locations and returns posterior draws of the ensemble weights, one
#' draw per stored MCMC iteration. For intercept fits, the shared spatial
#' intercept (delta) is kriged from the same iteration, so weight and delta
#' draws are aligned by iteration.
#'
#' @param pg_fit Ensemble fit object created with pg_ensemble()
#' @param X Matrix of covariates (S_pred x P), one row per unique prediction
#'   location, ordered by space_id
#' @param space_id Spatial location ID vector (N_pred)
#' @param coords Matrix of prediction x y coordinates, with
#'   colnames(coords) == c("x", "y"), (N_pred, 2)
#' @param debug If TRUE, also return draws of the kriged stick-breaking
#'   logits. Default FALSE.
#'
#' @return A list with elements
#'   \describe{
#'     \item{weight_preds}{Array (n_iter x M x S_pred) of weight draws.}
#'     \item{locs}{Data frame of unique prediction locations (space_id, x, y),
#'       sorted by space_id, in the order of the third dimension of
#'       weight_preds.}
#'     \item{delta_pred}{Matrix (n_iter x S_pred) of spatial intercept draws.
#'       Present only when the fit used intercept = TRUE.}
#'     \item{logit_preds}{Array (n_iter x (M-1) x S_pred) of stick-breaking
#'       logit draws. Present only when debug = TRUE.}
#'   }
#'
#' @examples
#' # pg_weight_pred()
#'
#' @export
pg_weight_pred <- function(pg_fit, X, space_id, coords, debug = FALSE) {
    X_pred <- as.matrix(X)
    X_obs <- pg_fit$X
    psi <- pg_fit$psi
    betas <- pg_fit$betas
    tau2 <- pg_fit$tau2
    rho <- pg_fit$rho

    nngp <- isTRUE(pg_fit$nngp)
    cov_kern <- if (is.null(pg_fit$cov_kern)) exponential_kernel else pg_fit$cov_kern

    # unique observed locations (order must match pg_fit$distmat / psi columns)
    locs <- pg_fit$locs

    pred_locs <- unique(cbind(space_id = space_id, coords)) |>
        as.data.frame()

    pred_locs <- pred_locs[order(pred_locs$space_id), ]

    M <- dim(psi)[1] + 1
    n_iter <- dim(psi)[3]
    n_pred <- nrow(pred_locs)

    if (nrow(X_pred) != n_pred) {
        stop("X must have one row per unique prediction location (", n_pred,
             "), ordered by space_id.")
    }

    coords_pred <- as.matrix(pred_locs[, c("x", "y")])
    coords_full <- as.matrix(locs[, c("x", "y")])

    if (!nngp) {
        # distances: pred -> obs
        distmat_obs <- pg_fit$distmat
        dx <- outer(coords_pred[, 1], coords_full[, 1], "-")
        dy <- outer(coords_pred[, 2], coords_full[, 2], "-")
        distmat_obs_preds <- sqrt(dx^2 + dy^2)
    } else {
        nngp_info <- pg_fit$nngp_info
        ord <- nngp_info$coord_ordering
        neighbors_pred <- get_neighbors_ref(
            ordered_coords = nngp_info$ordered_coords,
            pred_coords = coords_pred,
            m = nngp_info$number_neighbors
            )
        dist_matrices_pred <- get_dist_matrices_ref(
            ordered_coords = nngp_info$ordered_coords,
            coords_pred = coords_pred,
            neighbors = neighbors_pred
            )
        X_obs_ord <- X_obs[ord, , drop = FALSE]
        S_train <- nrow(X_obs_ord)
    }

    intercept <- isTRUE(pg_fit$intercept)
    if (intercept) {
        delta <- pg_fit$delta
        tau2_delta <- pg_fit$tau2_delta
        rho_delta <- pg_fit$rho_delta
        delta_pred_mat <- matrix(0, nrow = n_iter, ncol = n_pred)
    }

    # calculate weights
    weight_preds <- array(0, dim = c(n_iter, M, n_pred))
    if (debug) {
        logit_preds <- array(0, dim = c(n_iter, M - 1, n_pred))
    }

    for (t in 1:n_iter) {

        mu_pred_mat <- matrix(0, nrow = n_pred, ncol = M - 1)

        if (!nngp) {

            sigma12 <- tau2[t] * cov_kern(distance = distmat_obs_preds, theta = rho[t])
            sigma22 <- tau2[t] * cov_kern(distance = distmat_obs, theta = rho[t])
            A <- sigma12 %*% solve(sigma22)
            kriged_var_pred <- pmax(tau2[t] - rowSums(A * sigma12), 0)

            for (j in 1:(M - 1)) {
                kriged_mu_pred <- as.vector(X_pred %*% betas[j, , t]) +
                    A %*% (psi[j, , t] - as.vector(X_obs %*% betas[j, , t]))

                mu_pred_mat[, j] <- stats::rnorm(n_pred, kriged_mu_pred, sqrt(kriged_var_pred))
            }

            if (intercept) {
                s12_d <- tau2_delta[t] * cov_kern(distance = distmat_obs_preds,
                                                  theta = rho_delta[t])
                s22_d <- tau2_delta[t] * cov_kern(distance = distmat_obs,
                                                  theta = rho_delta[t])
                A_delta <- s12_d %*% solve(s22_d)
                # conditional mean (delta is zero-mean, no X beta add-back)
                mu_delta <- as.vector(A_delta %*% delta[, t])
                # conditional variance: diag(Sigma11 - A Sigma21)
                var_delta <- tau2_delta[t] - rowSums(A_delta * s12_d)
                var_delta <- pmax(var_delta, 0)
                delta_pred_mat[t, ] <- mu_delta + stats::rnorm(n_pred, 0, sqrt(var_delta))
            }

        } else {

            psi_t <- psi_slice(psi, t, M - 1)
            beta_t <- beta_slice(betas, t, M - 1)

            krige_psi <- nngp_krige(
                field_ord = psi_t[ord, , drop = FALSE],
                mu_ord = X_obs_ord %*% beta_t,
                mu_pred = X_pred %*% beta_t,
                tau = tau2[t],
                theta = rho[t],
                neighbors_pred = neighbors_pred,
                dist_matrices_pred = dist_matrices_pred,
                cov_kern = cov_kern
                )
            kriged_mu_pred <- krige_psi$mean
            kriged_var_pred <- krige_psi$var
            mu_pred_mat <- stats::rnorm(
                n_pred * (M - 1), 
                kriged_mu_pred, 
                sqrt(kriged_var_pred)
                ) |>
                matrix(nrow = n_pred, ncol = M - 1)

            if (intercept) {
                krige_delta <- nngp_krige(
                    field_ord = matrix(delta[ord, t], ncol = 1),
                    mu_ord = matrix(0, nrow = S_train, ncol = 1),
                    mu_pred = matrix(0, nrow = n_pred, ncol = 1),
                    tau = tau2_delta[t],
                    theta = rho_delta[t],
                    neighbors_pred = neighbors_pred,
                    dist_matrices_pred = dist_matrices_pred,
                    cov_kern = cov_kern
                    )
                delta_pred_mat[t, ] <- krige_delta$mean[, 1] +
                    stats::rnorm(n_pred, 0, sqrt(krige_delta$var))
            }
        }

        v_list <- vector("list", M - 1)
        w_list <- vector("list", M)

        for (j in 1:(M - 1)) {
            v_list[[j]] <- ilogit(mu_pred_mat[, j])
            if (debug) {
                logit_preds[t, j, ] <- as.vector(mu_pred_mat[, j])
            }
        }

        w_list[[1]] <- v_list[[1]]
        weight_preds[t, 1, ] <- w_list[[1]]
        for (j in 2:M) {
            w_list[[j]] <- Reduce("*", lapply(1:(j - 1), function(k) 1 - v_list[[k]]))
            if (j < M) {
                w_list[[j]] <- w_list[[j]] * v_list[[j]]
            }
            weight_preds[t, j, ] <- w_list[[j]]
        }

        if (t %% 100 == 0) {
            print(paste("Prediction iteration", t, "of", n_iter))
        }
    }



    dimnames(weight_preds) <- list(NULL, pg_fit$model_names, NULL)



    to_return <- list(
        weight_preds = weight_preds,
        locs = pred_locs
    )
    if (intercept) {
        to_return$delta_pred <- delta_pred_mat
    }
    if (debug) {
        to_return$logit_preds <- logit_preds
    }


    return(to_return)
}


#' Ensemble predictions at new locations
#'
#' Computes ensemble predictions by combining individual model
#' predictive distributions with spatially-varying weights kriged from a
#' fitted ensemble object.
#'
#'
#' @param pg_fit Ensemble fit object created with pg_ensemble()
#' @param X Matrix of covariates (S_pred x P), one row per unique prediction
#'   location, ordered by space_id
#' @param space_id Spatial location ID vector (N_pred)
#' @param coords Matrix of prediction x y coordinates, with
#'   colnames(coords) == c("x", "y"), (N_pred, 2)
#' @param model_est Matrix or data frame of model predictive means
#'   (N_pred x M). If column names contain all of pg_fit$model_names,
#'   columns are matched by name (any order); otherwise (e.g. raw data
#'   column names, or no names) columns are assumed to be in the same
#'   order as pg_fit$model_names.
#' @param model_sd Matrix or data frame of model predictive standard
#'   deviations (N_pred x M). Columns are paired with model_est
#'   positionally; its column names, if any, are ignored. Any reordering
#'   applied to model_est is applied to model_sd as well.
#' @param debug If TRUE, also return the intermediary kriged surfaces used
#'   internally: the stick-breaking logits (logit1, ..., logit(M-1)) and,
#'   for an intercept fit, the spatial intercept (delta, delta_sd).
#'   Default FALSE.
#'
#' @return A data frame with one row per prediction row, containing
#'   space_id, x, y, the ensemble prediction (pred), the ensemble
#'   predictive sd (pred_sd), and the posterior mean kriged weight for
#'   each model in columns named "\{model_name\}_weight". When
#'   debug = TRUE, additional columns logit1..logit(M-1), and delta,
#'   delta_sd for an intercept fit.
#'
#' @examples
#' # pg_pred()
#'
#' @export
pg_pred <- function(pg_fit, X, space_id, coords, model_est, model_sd, debug = FALSE) {
    model_names <- pg_fit$model_names

    aligned <- align_est_sd(model_est, model_sd, model_names = model_names)
    model_est <- aligned$model_est
    model_sd <- aligned$model_sd

    if (nrow(model_est) != length(space_id)) {
        stop("model_est and model_sd must have one row per element of space_id.")
    }

    weights_pred <- pg_weight_pred(
        pg_fit = pg_fit,
        X = X,
        space_id = space_id,
        coords = coords,
        debug = debug
    )

    idx <- match(space_id, weights_pred$locs$space_id)

    M <- length(model_names)
    n_iter <- dim(weights_pred$weight_preds)[1]
    sd2 <- model_sd^2
    Ey <- numeric(length(space_id))
    Ey2 <- numeric(length(space_id))

    for (t in seq_len(n_iter)) {
        W <- t(matrix(weights_pred$weight_preds[t, , idx], nrow = M))  # N_pred x M
        mu <- if (pg_fit$intercept) {
            model_est + weights_pred$delta_pred[t, idx]
        } else {
            model_est
        }
        Ey <- Ey + rowSums(W * mu)
        Ey2 <- Ey2 + rowSums(W * (sd2 + mu^2))
    }

    pred_mean <- Ey / n_iter
    pred_sd <- sqrt(pmax(Ey2 / n_iter - pred_mean^2, 0))
    out <- data.frame(
        space_id = space_id,
        x = coords[, "x"],
        y = coords[, "y"],
        pred = pred_mean,
        pred_sd = pred_sd
    )


    w_mean <- apply(weights_pred$weight_preds, c(3, 2), mean)   # S_pred x 
    w_mean <- matrix(w_mean, ncol = M)[idx, , drop = FALSE]
    colnames(w_mean) <- paste0(model_names, "_weight")
    out <- cbind(out, as.data.frame(w_mean))

    if (pg_fit$intercept) {
        out$delta <- colMeans(weights_pred$delta_pred)[idx]
     }
    if (debug) {
        l_mean <- apply(weights_pred$logit_preds, c(3, 2), mean)    # S_pred x (M-1)
        l_mean <- matrix(l_mean, ncol = M - 1)[idx, , drop = FALSE]
        colnames(l_mean) <- paste0("logit", seq_len(M - 1))
        out <- cbind(out, as.data.frame(l_mean))
    }

    return(out)
}


logit <- function(x) log(x / (1 - x))
ilogit <- function(x) 1 / (1 + exp(-x))
