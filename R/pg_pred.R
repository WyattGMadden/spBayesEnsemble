
#' Fit the Geostatistical Regression Model (GRM)
#'
#' This function fits a spatial Bayesian ensemble using model densities and (optional) ensemble weight covariates
#'
#' @param pg_fit Ensemble fit object created with pg_ensemble()
#' @param X Matrix of covariates (N_pred x P)
#' @param space_id spatial location ID vector (N_pred) 
#' @param coords Matrix of prediction x y coordinates, with colnames(coords) == c("x", "y"), (N_pred, 2)
#'
#' @return A list containing MCMC output 
#'
#' @examples
#' # pg_pred()
#' 
#' 
#' @export
pg_pred <- function(pg_fit, X, space_id, coords) {
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

    pred_locs <- unique(cbind(space_id = space_id, coords))

    # distances: pred -> obs
    distmat_obs <- pg_fit$distmat
    coords_pred <- as.matrix(pred_locs[, c("x", "y")])
    coords_full <- as.matrix(locs[, c("x", "y")])

    dx <- outer(coords_pred[, 1], coords_full[, 1], "-")
    dy <- outer(coords_pred[, 2], coords_full[, 2], "-")
    distmat_obs_preds <- sqrt(dx^2 + dy^2)

    M <- dim(psi)[1] + 1
    S_obs <- dim(psi)[2]
    n_iter <- dim(psi)[3]
    n_pred <- nrow(pred_locs)

    # accumulate posterior predictive mean of weights
    wsum <- matrix(
        0,
        nrow = M,
        ncol = n_pred
    )

    for (t in 1:n_iter) {

        sigma12 <- tau2[t] * exp(-distmat_obs_preds / rho[t])
        sigma22 <- tau2[t] * exp(-distmat_obs / rho[t])

        A <- sigma12 %*% solve(sigma22)

        mu_list <- vector("list", M - 1)
        eta_list <- vector("list", M - 1)
        v_list <- vector("list", M - 1)
        w_list <- vector("list", M)

        for (j in 1:(M - 1)) {
            mu_list[[j]] <- as.vector(X_pred %*% betas[j, , t]) +
                A %*% (psi[j, , t] - as.vector(X_obs %*% betas[j, , t]))

            eta_list[[j]] <- mu_list[[j]]
            v_list[[j]] <- ilogit(eta_list[[j]])

        }

        #w1 <- v1
        #w2 <- v2 * (1 - v1)
        #w3 <- (1 - v1) * (1 - v2)
        w_list[[1]] <- v_list[[1]]
        wsum[1, ] <- wsum[1, ] + w_list[[1]]
        for (j in 2:M) {
            w_list[[j]] <- Reduce("*", lapply(1:(j-1), function(k) 1 - v_list[[k]]))
            if (j < M) {
                w_list[[j]] <- w_list[[j]] * v_list[[j]]
            }
            wsum[j, ] <- wsum[j, ] + w_list[[j]]
        }

        if (t %% 100 == 0) {
            print(paste("Prediction iteration", t, "of", n_iter))
        }
    }

    wmean <- wsum / n_iter
    wmean <- t(wmean)
    colnames(wmean) <- pg_fit$model_names
    weights_pred <- as.data.frame(wmean)

    weights_pred$x <- pred_locs$x
    weights_pred$y <- pred_locs$y
    weights_pred$space_id <- pred_locs$space_id

    return(weights_pred)
}

