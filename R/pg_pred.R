#' Predict spatially-varying ensemble weights at new locations
#'
#' Kriges the latent stick-breaking GPs (psi) from a fitted ensemble object
#' to new locations and returns posterior mean ensemble weights.
#'
#' @param pg_fit Ensemble fit object created with pg_ensemble()
#' @param X Matrix of covariates (N_pred x P)
#' @param space_id Spatial location ID vector (N_pred)
#' @param coords Matrix of prediction x y coordinates, with
#'   colnames(coords) == c("x", "y"), (N_pred, 2)
#'
#' @return A data frame of posterior mean weights (one column per model,
#'   named by pg_fit$model_names), with x, y, and space_id columns
#'   appended. One row per unique prediction location.
#'
#' @examples
#' # pg_weight_pred()
#'
#' @export
pg_weight_pred <- function(pg_fit, X, space_id, coords) {
    X_pred <- X
    X_obs <- pg_fit$X
    psi <- pg_fit$psi
    betas <- pg_fit$betas
    tau2 <- pg_fit$tau2
    rho <- pg_fit$rho

    # unique observed locations (order must match pg_fit$distmat / psi columns)
    locs <- pg_fit$coords |>
        unique() |>
        as.data.frame()
    pred_locs <- unique(cbind(space_id = space_id, coords)) |>
        as.data.frame()

    # distances: pred -> obs
    distmat_obs <- pg_fit$distmat
    coords_pred <- as.matrix(pred_locs[, c("x", "y")])
    coords_full <- as.matrix(locs[, c("x", "y")])
    dx <- outer(coords_pred[, 1], coords_full[, 1], "-")
    dy <- outer(coords_pred[, 2], coords_full[, 2], "-")
    distmat_obs_preds <- sqrt(dx^2 + dy^2)

    M <- dim(psi)[1] + 1
    n_iter <- dim(psi)[3]
    n_pred <- nrow(pred_locs)

    intercept <- isTRUE(pg_fit$intercept)
    if (intercept) {
        delta <- pg_fit$delta
        tau2_delta <- pg_fit$tau2_delta
        rho_delta <- pg_fit$rho_delta
        delta_pred_mat <- matrix(0, nrow = n_pred, ncol = n_iter)
    }

    # accumulate posterior predictive mean of weights
    wsum <- matrix(0, nrow = M, ncol = n_pred)

    for (t in 1:n_iter) {
        sigma12 <- tau2[t] * exp(-distmat_obs_preds / rho[t])
        sigma22 <- tau2[t] * exp(-distmat_obs / rho[t])
        A <- sigma12 %*% solve(sigma22)

        if (intercept) {
            s12_d <- tau2_delta[t] * exp(-distmat_obs_preds / rho_delta[t])
            s22_d <- tau2_delta[t] * exp(-distmat_obs / rho_delta[t])
            A_delta <- s12_d %*% solve(s22_d)
            # conditional mean (delta is zero-mean, no X beta add-back)
            mu_delta <- as.vector(A_delta %*% delta[, t])
            # conditional variance: diag(Sigma11 - A Sigma21)
            var_delta <- tau2_delta[t] - rowSums(A_delta * s12_d)
            var_delta <- pmax(var_delta, 0)
            delta_pred_mat[, t] <- mu_delta + stats::rnorm(n_pred, 0, sqrt(var_delta))
        }

        v_list <- vector("list", M - 1)
        w_list <- vector("list", M)

        for (j in 1:(M - 1)) {
            mu_j <- as.vector(X_pred %*% betas[j, , t]) +
                A %*% (psi[j, , t] - as.vector(X_obs %*% betas[j, , t]))
            v_list[[j]] <- ilogit(mu_j)
        }

        w_list[[1]] <- v_list[[1]]
        wsum[1, ] <- wsum[1, ] + w_list[[1]]
        for (j in 2:M) {
            w_list[[j]] <- Reduce("*", lapply(1:(j - 1), function(k) 1 - v_list[[k]]))
            if (j < M) {
                w_list[[j]] <- w_list[[j]] * v_list[[j]]
            }
            wsum[j, ] <- wsum[j, ] + w_list[[j]]
        }

        if (t %% 100 == 0) {
            print(paste("Prediction iteration", t, "of", n_iter))
        }
    }

    wmean <- t(wsum / n_iter)
    colnames(wmean) <- pg_fit$model_names

    weights_pred <- as.data.frame(wmean)
    weights_pred$x <- pred_locs$x
    weights_pred$y <- pred_locs$y
    weights_pred$space_id <- pred_locs$space_id

    if (intercept) {
        weights_pred$delta <- rowMeans(delta_pred_mat)
        weights_pred$delta_sd <- apply(delta_pred_mat, 1, stats::sd)
    }

    return(weights_pred)
}


#' Ensemble predictions at new locations
#'
#' Computes ensemble predictions by combining individual model
#' predictive distributions with spatially-varying weights kriged from a
#' fitted ensemble object.
#'
#'
#' @param pg_fit Ensemble fit object created with pg_ensemble()
#' @param X Matrix of covariates (N_pred x P)
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
#'
#' @return A data frame with one row per prediction row, containing
#'   space_id, x, y, the ensemble prediction (pred), the ensemble
#'   predictive sd (pred_sd), and the kriged weight for each model in
#'   columns named "\{model_name\}_weight".
#'
#' @examples
#' # pg_pred()
#'
#' @export
pg_pred <- function(pg_fit, X, space_id, coords, model_est, model_sd) {
    model_names <- pg_fit$model_names

    aligned <- align_est_sd(model_est, model_sd, model_names = model_names)
    model_est <- aligned$model_est
    model_sd <- aligned$model_sd

    if (nrow(model_est) != length(space_id)) {
        stop("model_est and model_sd must have one row per element of space_id.")
    }

    # kriged posterior-mean weights at unique prediction locations
    weights_pred <- pg_weight_pred(
        pg_fit = pg_fit,
        X = X,
        space_id = space_id,
        coords = coords
    )

    intercept <- isTRUE(pg_fit$intercept)
    if (intercept) {
        idx <- match(space_id, weights_pred$space_id)
        delta_row <- weights_pred$delta[idx]
        delta_sd_row <- weights_pred$delta_sd[idx]
    } else {
        delta_row <- 0
        delta_sd_row <- 0
    }

    # map location-level weights to prediction rows
    w <- as.matrix(
        weights_pred[match(space_id, weights_pred$space_id),
                     model_names,
                     drop = FALSE]
    )

    # mixture mean and variance:
    # E[y] = sum_k w_k mu_k
    # Var[y] = sum_k w_k (sd_k^2 + mu_k^2) - E[y]^2
    shifted_est <- model_est + delta_row
    pred_mean <- rowSums(w * shifted_est)
    pred_var <- rowSums(w * (model_sd^2 + shifted_est^2)) - pred_mean^2 + delta_sd_row^2
    pred_sd <- sqrt(pmax(pred_var, 0))

    out <- data.frame(
        space_id = space_id,
        x = coords[, "x"],
        y = coords[, "y"],
        pred = pred_mean,
        pred_sd = pred_sd
    )
    colnames(w) <- paste0(model_names, "_weight")
    out <- cbind(out, as.data.frame(w))

    return(out)
}


logit <- function(x) log(x / (1 - x))
ilogit <- function(x) 1 / (1 + exp(-x))
