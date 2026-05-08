#' Fit spatial Bayesian ensemble
#'
#' This function fits a spatial Bayesian ensemble using model densities and (optional) ensemble weight covariates
#'
#' @param dens Matrix of model densities (N x M)
#' @param X Matrix of covariates (N x P)
#' @param space_id Spatial location ID vector (N) 
#' @param coords Matrix of x y coordinates, with colnames(coords) == c("x", "y"), (N, 2)
#' @param n_iter Number of iterations used in MCMC
#'
#' @return A list containing MCMC output 
#'
#' @examples
#' # pg_ensemble()
#' 
#' 
#' @export
pg_ensemble <- function(dens, X, space_id, coords, n_iter) {

    model_names <- colnames(dens)
    dens <- t(dens)
    locs <- cbind(space_id = space_id, coords) |>
        unique() |>
        as.data.frame()
    distmat <- stats::dist(locs[, c("x", "y")]) |>
        as.matrix()
    P <- ncol(X)
    S <- nrow(locs)
    M <- nrow(dens)


    weights_all <- array(1, dim = c(M, S, n_iter))

    #weights_all[, , 1] <- counts_by_spat[, c("weight1", "weight2", "weight3")] |>
    #    as.matrix()
    weights_all[, , 1] <- matrix(1/M, nrow = M, ncol = S)

    #psi1 <- matrix(0, nrow = S, ncol = n_iter)
    #psi2 <- matrix(0, nrow = S, ncol = n_iter)
    psi <- array(0, dim = c(M - 1, S, n_iter))

    betas <- array(0, dim = c(M - 1, P, n_iter))
    betas[, , 1] <- 0
    betas_mu <- 0
    betas_sd <- 10

    #init gp params
    tau2 <- rep(0, n_iter)
    tau2[1] <- 1
    tau2_a <- 0.001
    tau2_b <- 0.001

    rho <- rep(0, n_iter)
    rho[1] <- 20
    rho_step <- 2
    rho_mu <- 3
    rho_sd <- 1

    for (i in 2:n_iter) {


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
        z_counts_rev_cum <- z_counts_rev_cum[, 1:(ncol(z_counts_rev_cum) - 1)]


        #######################################
        ########weights sample##########
        #######################################

        #calculate gp cov
        covar <- tau2[i - 1] * exp(-1 / rho[i - 1] * distmat)

        #draw pg variables
        pg_vars <- matrix(0, nrow = M - 1, ncol = S)
        for (j in 1:(M - 1)) {
            pg_vars[j, ] <- BayesLogit::rpg.devroye(num = S,
                                                    h = z_counts_rev_cum[, j],
                                                    z = psi[j, , i - 1])
        }
        
        k_z <- matrix(0, nrow = M - 1, ncol = S)
        for (j in 1:(M - 1)) {
            k_z[j, ] <- z_counts[, j] - (z_counts_rev_cum[, j]) / 2
        }


        MMM <- matrix(0, nrow = M - 1, ncol = S)
        for (j in 1:(M - 1)) {
            omega <- diag(pg_vars[j, ])
            Sigma <- solve(omega + solve(covar))
            MMM[j, ] <- Sigma %*% (k_z[j, ] + solve(covar) %*% (X %*% betas[j, , i - 1]))
            psi[j, , i] <- t(mvtnorm::rmvnorm(n = 1, mean = MMM[j, ], sigma = Sigma))
        }


        #update betas
        for (j in 1:(M - 1)) {
            var_beta <- solve(t(X) %*% solve(covar) %*% X + diag(1 / betas_sd^2, P))
            mean_beta <- var_beta %*% (t(X) %*% solve(covar) %*% psi[j, , i])
            betas[j, , i] <- t(mvtnorm::rmvnorm(n = 1, mean = mean_beta, sigma = var_beta))

        }
        #update tau
        SSS <- 0
        for (j in 1:(M - 1)) {
            SSS_j <- t(psi[j, , i] - X %*% betas[j, , i]) %*% solve(exp(-distmat / rho[i - 1])) %*% (psi[j, , i] - X %*% betas[j, , i])
            SSS <- SSS + SSS_j
        }
        SSS <- SSS / 2


        tau2[i] <- 1 / stats::rgamma(1, (S * (M - 1)) / 2 + tau2_a, SSS + tau2_b)
        covar <- tau2[i] * exp(-distmat / rho[i - 1])
        #tau2_all[i] <- 1

        #update rho
        rho_prop <- stats::rlnorm(1, log(rho[i - 1]), rho_step)
        SSS_curr <- covar
        SSS_prop <- tau2[i] * exp(-distmat / rho_prop)

        lik_curr <- 0
        lik_prop <- 0
        for (j in 1:(M - 1)) {
            lik_curr <- lik_curr + mvtnorm::dmvnorm(t(psi[j, , i]), X %*% betas[j, , i], SSS_curr, log = TRUE)
            lik_prop <- lik_prop + mvtnorm::dmvnorm(t(psi[j, , i]), X %*% betas[j, , i], SSS_prop, log = TRUE)
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
        logit_psi <- ilogit(psi[, , i])
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


        if (i %% 100 == 0) {
            print(i)
        }

    }

    return(
        list(
            weights_all = weights_all, 
            psi = psi, 
            betas = betas, 
            tau2 = tau2, 
            rho = rho,
            coords = coords,
            distmat = distmat,
            X = X,
            model_names = model_names
        )
    )
}

