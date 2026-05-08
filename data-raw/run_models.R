library(tidyverse)
bart_cv <- function(dat,
                   Y_name,
                   X_names,
                   cv.object,
                   space.id, 
                   time.id, 
                   spacetime.id,
                   just.fit.i = NULL) {

    cv.id <- cv.object$cv.id

    Y.cv <- data.frame(time_id = time.id, 
                      space_id = space.id, 
                      spacetime_id = spacetime.id,
                      obs = dat[, Y_name][[1]],
                      estimate = NA, 
                      sd = NA)
  
  
    for (cv.i in 1:cv.object$num.folds) {
        
        #just.fit.i override
        if (!is.null(just.fit.i)) cv.i <- just.fit.i
    
        print(paste0("Performing CV Experiment ---- Fold ", cv.i))

        if (cv.object$type == "spatial_buffered") {

            train.id <- cv.id != cv.i & cv.id != 0 & (!cv.object$drop.matrix[, cv.i])

        } else {

            train.id <- cv.id != cv.i & cv.id != 0

        }

        test.id.temp <- cv.id == cv.i
        #remove any test observations that are not within the training observations time range
        #these will be NA's in final cv predictions
        test.id.remove <- (min(time.id[train.id]) > time.id | max(time.id[train.id]) < time.id) & test.id.temp
        test.id <- test.id.temp & !test.id.remove

        time.id.train <- time.id[train.id]
        time.id.test <- time.id[test.id]

    
        space.id.train <- space.id[train.id]
        space.id.test <- space.id[test.id]

        #grm requires space.id to be from 1:max(space_id)
        #spatial cross validation breaks this assumption (missing space_id values)
        #here we create temporary space.id values that are from 1:length(unique(space.id))
        space.id.train.key <- sort(unique(space.id.train))
        space.id.test.key <- sort(unique(space.id.test))
        space.id.train.temp <- sapply(space.id.train, 
                                     function(x) which(space.id.train.key == x))
        space.id.test.temp <- sapply(space.id.test,
                                    function(x) which(space.id.test.key == x))
        spacetime.id.train <- spacetime.id[train.id]
        spacetime.id.test <- spacetime.id[test.id]
   
        fit.cv <- dbarts::bart2(formula = as.formula(paste0(Y_name, " ~ . ")),
                        n.samples = 500,
                        data = dat[train.id, c(Y_name, X_names)],
                        n.cuts = 100,
                        keepTrees = T)
    

        cv.test <- predict(fit.cv, 
                           newdata = dat[test.id, X_names])



        sigma2_bart <- mean(fit.cv$sigma^2)


        Y.cv$estimate[test.id] <- colMeans(cv.test)
        Y.cv$sd[test.id] <- sqrt(apply(cv.test, 2, var) + sigma2_bart)

        if (!is.null(just.fit.i)) break

    }
 
    Y.cv$upper_95 <- Y.cv$estimate + 1.96 * Y.cv$sd
    Y.cv$lower_95 <- Y.cv$estimate - 1.96 * Y.cv$sd
  
    return(Y.cv)
}

################################
###prepare model fitting data###
################################

load("~/ensembleDownscaleR/data/modis_aqs_matched.rda")
la_aqs <- modis_aqs_matched |>
    mutate(time_id = as.integer(as.factor(time_id)),
           space_id = as.integer(as.factor(space_id)),
           spacetime_id = as.integer(as.factor(spacetime_id)))


##########################
###fit component models###
##########################
cv_info <- ensembleDownscaleR::create_cv(
    space.id = la_aqs$space_id,
    time.id = la_aqs$time_id, 
    spacetime.id = la_aqs$spacetime_id,
    type = "spatial"
)


cmaq_cv <- ensembleDownscaleR::grm_cv(
    Y = la_aqs$pm25,
    X = la_aqs$cmaq,
    L = la_aqs[, c("elevation", "forestcover",
                   "hwy_length", "lim_hwy_length", 
                   "local_rd_length", "point_emi_any")],
    M = la_aqs[, c("tmp", "wind")],
    coords = la_aqs[, c("x", "y")],
    space.id = la_aqs$space_id,
    time.id = la_aqs$time_id,
    spacetime.id = la_aqs$spacetime_id,
    cv.object = cv_info,
    n.iter = 5000,
    burn = 1000,
    thin = 4
)

la_aqs$cmaqgrm_model_estimate <- cmaq_cv$estimate
la_aqs$cmaqgrm_model_sd <- cmaq_cv$sd
la_aqs$cmaqgrm_model_density <- dnorm(
    la_aqs$pm25, 
    mean = la_aqs$cmaqgrm_model_estimate, 
    sd = la_aqs$cmaqgrm_model_sd
)



modis_cv <- ensembleDownscaleR::grm_cv(
    Y = la_aqs$pm25,
    X = la_aqs$aod,
    L = la_aqs[, c("elevation", "forestcover",
                   "hwy_length", "lim_hwy_length", 
                   "local_rd_length", "point_emi_any")],
    M = la_aqs[, c("tmp", "wind", "cmaq", "tempaod", 
                   "windaod", "elevationaod")],
    coords = la_aqs[, c("x", "y")],
    space.id = la_aqs$space_id,
    time.id = la_aqs$time_id,
    spacetime.id = la_aqs$spacetime_id,
    cv.object = cv_info,
    n.iter = 5000,
    burn = 1000,
    thin = 4
)


la_aqs$modisgrm_model_estimate <- modis_cv$estimate
la_aqs$modisgrm_model_sd <- modis_cv$sd
la_aqs$modisgrm_model_density <- dnorm(
    la_aqs$pm25, 
    mean = la_aqs$modisgrm_model_estimate, 
    sd = la_aqs$modisgrm_model_sd
)



Y_name <- "pm25"
X_names <- c("cmaq", "aod", "elevation", "forestcover",
             "hwy_length", "lim_hwy_length", "local_rd_length", 
             "point_emi_any", "tmp", "wind", "tempaod", 
             "windaod", "elevationaod",
             "x", "y")
bart_fit_cv <- bart_cv(
    dat = la_aqs[, c(Y_name, X_names)],
    Y_name = Y_name,
    X_names = X_names,
    cv.object = cv_info,
    space.id = la_aqs$space_id,
    time.id = la_aqs$time_id,
    spacetime.id = la_aqs$spacetime_id
)
la_aqs$bart_model_estimate <- bart_fit_cv$estimate
la_aqs$bart_model_sd <- bart_fit_cv$sd
la_aqs$bart_model_density <- dnorm(
    la_aqs$pm25, 
    mean = la_aqs$bart_model_estimate, 
    sd = la_aqs$bart_model_sd
)

la_aqs <- la_aqs |>
    filter(!is.na(modisgrm_model_estimate),
           !is.na(cmaqgrm_model_estimate),
           !is.na(bart_model_estimate))

la_aqs$time_id <- as.integer(as.factor(la_aqs$time_id))
la_aqs$space_id <- as.integer(as.factor(la_aqs$space_id))
la_aqs$spacetime_id <- as.integer(as.factor(la_aqs$spacetime_id))

#rmse check
sqrt(mean((la_aqs$cmaqgrm_model_estimate - la_aqs$pm25)^2))
sqrt(mean((la_aqs$modisgrm_model_estimate - la_aqs$pm25)^2))
sqrt(mean((la_aqs$bart_model_estimate - la_aqs$pm25)^2))

#coverage check
mean(la_aqs$pm25 > la_aqs$cmaqgrm_model_estimate - 1.96 * la_aqs$cmaqgrm_model_sd &
         la_aqs$pm25 < la_aqs$cmaqgrm_model_estimate + 1.96 * la_aqs$cmaqgrm_model_sd)
mean(la_aqs$pm25 > la_aqs$modisgrm_model_estimate - 1.96 * la_aqs$modisgrm_model_sd &
            la_aqs$pm25 < la_aqs$modisgrm_model_estimate + 1.96 * la_aqs$modisgrm_model_sd)
mean(la_aqs$pm25 > la_aqs$bart_model_estimate - 1.96 * la_aqs$bart_model_sd &
            la_aqs$pm25 < la_aqs$bart_model_estimate + 1.96 * la_aqs$bart_model_sd)

usethis::use_data(la_aqs, overwrite = TRUE)

#############################
###prepare prediction data###
#############################
load("../data/la_aqs.rda")
cmaq_fit <- ensembleDownscaleR::grm(
    Y = la_aqs$pm25,
    X = la_aqs$cmaq,
    L = la_aqs[, c("elevation", "forestcover",
                   "hwy_length", "lim_hwy_length", 
                   "local_rd_length", "point_emi_any")],
    M = la_aqs[, c("tmp", "wind")],
    coords = la_aqs[, c("x", "y")],
    space.id = la_aqs$space_id,
    time.id = la_aqs$time_id,
    spacetime.id = la_aqs$spacetime_id,
    n.iter = 5000,
    burn = 1000,
    thin = 4
)

modis_fit <- ensembleDownscaleR::grm(
    Y = la_aqs$pm25,
    X = la_aqs$aod,
    L = la_aqs[, c("elevation", "forestcover",
                   "hwy_length", "lim_hwy_length", 
                   "local_rd_length", "point_emi_any")],
    M = la_aqs[, c("tmp", "wind", "cmaq", "tempaod", 
                   "windaod", "elevationaod")],
    coords = la_aqs[, c("x", "y")],
    space.id = la_aqs$space_id,
    time.id = la_aqs$time_id,
    spacetime.id = la_aqs$spacetime_id,
    n.iter = 5000,
    burn = 1000,
    thin = 4
)

bart_fit <- dbarts::bart2(formula = as.formula(paste0(Y_name, " ~ . ")),
                        data = la_aqs[, c(Y_name, X_names)],
                        keepTrees = T)



load("~/ensembleDownscaleR/data/modis_full.rda")
la_grid <- modis_full

##preds
la_grid <- la_grid |>
    mutate(time_id = as.integer(as.factor(time_id)),
           space_id = as.integer(as.factor(space_id)),
           spacetime_id = as.integer(as.factor(spacetime_id)))


cmaq_pred <- ensembleDownscaleR::grm_pred(
    grm.fit = cmaq_fit,
    X = la_grid$cmaq,
    L = la_grid[, c("elevation", "forestcover",
                   "hwy_length", "lim_hwy_length", 
                   "local_rd_length", "point_emi_any")],
    M = la_grid[, c("tmp", "wind")],
    coords = la_grid[, c("x", "y")],
    space.id = la_grid$space_id,
    time.id = la_grid$time_id,
    spacetime.id = la_grid$spacetime_id,
    verbose = T
)

modis_pred <- ensembleDownscaleR::grm_pred(
    grm.fit = modis_fit,
    X = la_grid$aod,
    L = la_grid[, c("elevation", "forestcover",
                   "hwy_length", "lim_hwy_length", 
                   "local_rd_length", "point_emi_any")],
    M = la_grid[, c("tmp", "wind", "cmaq", "tempaod", 
                   "windaod", "elevationaod")],
    coords = la_grid[, c("x", "y")],
    space.id = la_grid$space_id,
    time.id = la_grid$time_id,
    spacetime.id = la_grid$spacetime_id,
    verbose = T
)

bart_pred <- predict(
    bart_fit,
    newdata = la_grid[, X_names]
)

la_grid$cmaqgrm_model_estimate <- cmaq_pred$estimate
la_grid$cmaqgrm_model_sd <- cmaq_pred$sd

la_grid$modisgrm_model_estimate <- modis_pred$estimate
la_grid$modisgrm_model_sd <- modis_pred$sd

sigma2_bart <- mean(bart_fit$sigma^2)
la_grid$bart_model_estimate <- colMeans(bart_pred)
la_grid$bart_model_sd <- sqrt(apply(bart_pred, 2, var) + sigma2_bart)


usethis::use_data(
    la_grid, 
    overwrite = TRUE
)
