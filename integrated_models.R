#' ---
#' title: "Integrated model demo"
#' author: "Saras Windecker & Nick Golding"
#' date: "21 June, 2022"
#' output: html_document
#' ---
#'

# run integrated density model
library(raster)

# mask the stack
stack_all <- readRDS('spatial-data/processed-layers/static_stack.rds')
# we want the tree as an offset rather than in the linear XB because
# we care about the effect of zero trees. If it was only in the XB it would
# move the intercept effect 
stack_all$nveg_alltrees_sq <- stack_all$nveg_alltrees^2
abund_covariates <- c('nveg_nativeshrub', 'annprecip', 
                      'ndvi', 'nveg_alltrees_sq') 
ncov <- length(abund_covariates)

# extract center and scale for full region
scaled_covariates <- scale(as.data.frame(stack_all)[, abund_covariates])
scale <- attr(scaled_covariates, 'scaled:scale')
center <- attr(scaled_covariates, 'scaled:center')

# add dist to coast? soil type? select a few that are relevant. urban cover?
# posterior predictive checks. or dharma. randomised quantile residuals. 
# simulate new counts from model. get for each 
# data point, and plot against each covariate. also plot in space. 
# can slo plot residuals to other covariates 
# not included to check them. 

# presence-only data
po <- targets::tar_read(density_data)
spatial_po  <- SpatialPoints(po)
extract_po <- as.data.frame(extract(stack_all, spatial_po))
new_po <- po[which(complete.cases(as.data.frame((extract_po)))),]

ppm <- ppmify::ppmify(new_po,
                      area = stack_all[[1]],
                      covariates = stack_all,
                      density = 1 / 10,
                      method = 'grid')

X1 <- scale(ppm[, abund_covariates], 
            center = center,
            scale = scale)
offset1 <- ppm[, 'nveg_alltrees']
# log(0) = -Inf
offset1[offset1 == 0] <- .0001

bias_covariates <- c('city_accessibility') 
ncov_bias <- length(bias_covariates)
R <- scale(ppm[, bias_covariates])

# spotlighting data
counts <- targets::tar_read(abundance_data)
# plot(stack_all[[1]])
# points(counts[,c('lon', 'lat')],
#        col = new_counts$cluster, pch = 16)
spatial_counts <- SpatialPoints(counts[, c('lon', 'lat')])

X2 <- scale(
  extract(stack_all[[abund_covariates]], 
          spatial_counts), 
  center = center,
  scale = scale)

offset2 <- extract(
  stack_all[['nveg_alltrees']], 
  spatial_counts)
offset2[offset2 == 0] <- .0001

# predictions across whole region
pred_extract <- as.data.frame(stack_all)
complete <- complete.cases(pred_extract)
X_pred_withNA <- pred_extract[complete,]

X_pred <- scale(X_pred_withNA[, abund_covariates],
                center = center,
                scale = scale)

offset_pred <- as.data.frame(X_pred_withNA[, 'nveg_alltrees'])
offset_pred <- as.vector(offset_pred[,1])
offset_pred[offset_pred == 0] <- .0001

# time to detection scat data
occ <- targets::tar_read(scat_data)
spatial_occ  <- SpatialPoints(occ[, c('lon', 'lat')])
extract_occ <- as.data.frame(extract(
  stack_all[[abund_covariates]], spatial_occ))
new_occ <- occ[which(complete.cases(as.data.frame((extract_occ)))),]
new_spatial_occ <- SpatialPoints(new_occ[, c('lon', 'lat')])

X3 <- scale(extract(stack_all[[abund_covariates]], new_spatial_occ), 
            center = center,
            scale = scale)
offset3 <- extract(stack_all[['nveg_alltrees']], new_spatial_occ)
offset3[offset3 == 0] <- .0001



n_lambda <- 100

win.data <- list(
  
  ncov = ncov,
  
  # presence-only data
  n_y1 = nrow(ppm),
  y1 = ppm$points,
  X1 = X1, 
  offset1 = offset1,
  R = R,
  ncov_bias = ncov_bias,
  area = ppm$weights,
  
  # spotlight count data
  n_y2 = nrow(new_counts),
  y2 = new_counts$count,
  X2 = X2,
  offset2 = offset2,
  
  # time to detection scat data
  n_y3 = nrow(new_occ),
  ttd = new_occ$ttd,
  d = new_occ$censoring_indicator,
  X3 = X3,
  offset3 = offset3,
  nminutes = 5,
  Tmax = 5.1,
  
  # prediction 
  n_pred = nrow(X_pred),
  X_pred = X_pred,
  offset_pred = offset_pred,
  n_lambda = n_lambda,
  lambda = seq(0.001, 20, length.out = n_lambda)
)




sink("explore_scripts/IM_IPPpo_PoissonCount_ttdScat_28May.txt")
cat("


# Define model and write into R working directory

model {

  ## Priors ##

  # Priors for abundance of possums
  beta0 ~ dnorm(0, 0.1)

  # Prior for slopes in abundance predictors
  for(b in 1:ncov) {
    beta[b] ~ dnorm(0, 0.1) #stdev 3.16
  }

  # Priors for possum spatial sampling bias
  beta_bias0 ~ dnorm(0, 0.1)

  # Prior for slopes in bias predictors
  for(b in 1:ncov_bias) {
    beta_bias[b] ~ dnorm(0, 0.1)
  }

  # Prior for scat detection
  beta_det0 ~ dnorm(0, 0.1)

  # Prior for abundance of scat
  # truncated positive because we expect scat abundance to increase 
  # with possum abundance 
  beta_scatabund0 ~ dnorm(0, 0.1)T(0,)
  beta_scatabund1 ~ dnorm(0, 0.1)T(0,)
  # truncated negative 
  beta_scatabund2 ~ dnorm(0, 1)T(,0)

  ## Likelihood ##

  # presence-only possum data is y1
  for (i in 1:n_y1){

   # note no variance parameter
   y1[i] ~ dpois(lambda1[i]*area[i]*bias[i])

   # where offset is alltrees because we want to vary possums by trees
   log(lambda1[i]) <- beta0 + log(offset1[i]) + inprod(beta[], X1[i, ])

   # our estimate is conditional on sampling effort (which is area of
   # the 50 m radius circle)
   # R is the cov matrix for bias
   log(bias[i]) <- beta_bias0 + inprod(beta_bias[], R[i,])

  }

# spotlight count data of possums is y2
for (i in 1:n_y2){

   y2[i] ~ dpois(lambda2[i])
   log(lambda2[i]) <- beta0 + log(offset2[i]) + inprod(beta[], X2[i, ])

}

# scat time to detection data is y3
for (i in 1:n_y3) {               # Loop over sites

  # True state model for the partially observed true state
  # True occupancy z at site i
  z[i] ~ dbern(psi[i])
  cloglog(psi[i]) <- log_abund_scat[i]

  ttd[i] ~ dcat(q[i,])

  # we added intercept here to allow the prob of scats for single possum to vary
  # we removed the polynomial because it was not monotonically increasing
  # and we would expect scats to only ever increase with more possums
  log(lambda3[i]) <- beta0 + log(offset3[i]) + inprod(beta[], X3[i,])
  log_abund_scat[i] <- beta_scatabund0 + beta_scatabund1*log(lambda3[i])
                        + beta_scatabund2*log(lambda3[i])^2
  log(encounter[i]) <- beta_det0 + log_abund_scat[i]

  # Observation model for the actual observations
  # Discretised exponential time to detection
  # we use dcat instead of dexp because of the need to discretise the
  # exponential distribution for time presented in minute intervals

  # integration of the exponential distribution of q probabilities
  # for each interval of nminutes
  for (k in 1:nminutes) {
    q[i,k] <- -exp(-encounter[i]*k) + exp(-encounter[i]*(k-1))
  }

  # Accomodation of z = 0 and censoring
  d[i] ~ dbern(theta[i])       # Model for censoring indicator
  theta[i] <- z[i] * step(ttd[i] - (Tmax)) + (1 - z[i])


  ## this is an alternative scat data model that ignores ttd. we could run
  ## this as an alternative and compare the improvement in predictive performance
  ## and/or the increase in explained variation.
  # # presence/absence model for observed scat (detection probability is wrapped
  # # up in observed scat - possum abundance model)
  # d[i] ~ dbern(1 - prob_observed_scat[i])
  # cloglog(prob_observed_scat[i]) <- log_abund_observed_scat[i]
  # log_abund_observed_scat[i] <- beta_scatabund1*log(lambda3[i])
  # log(lambda3[i]) <- beta0 + log(offset3[i]) + inprod(beta[], X3[i,])
  }

  ## Derived quantities
  for (j in 1:n_pred){
    log(pred[j]) <- beta0 + log(offset_pred[j]) + inprod(beta[], X_pred[j, ])
  }

  # predict relationship of scat abundance to possum abundance, and
  # scat enouncter rate to possum abundance. to plot relationships with seq. lambdas.
  for (i in 1:n_lambda){
    log(pred_abund_scat[i]) <- beta_scatabund0 + beta_scatabund1*log(lambda[i])
                                + beta_scatabund2*log(lambda[i])^2
    log(pred_encounter[i]) <- beta_det0 + log(pred_abund_scat[i])
  }
}


", fill = TRUE)
sink()

# # Initialize with z = 1 throughout and 
# # all missings due to censoring, rather than non-occurrence
# zst <- rep(1, times = nrow(new_occ))
# 
# ttdst <- as.numeric(new_occ$ttd)
# I <- which(is.na(new_occ$ttd))
# ttdst[I] <- 5 + 1
# ttdst[-I] <- NA
# 
# qst <- matrix(rep(0.5, times = nrow(new_occ)*5),
#               ncol = 5)
# 
# inits <- function(){ list(beta0 = rnorm(1), 
#                           beta = rnorm(ncov),
#                           beta_bias0 = rnorm(1),
#                           beta_bias = rnorm(ncov_bias),
#                           z = zst,
#                           ttd = ttdst#,
#                           beta_det0 = rnorm(1),
#                           beta_scatabund1 = rnorm(1),
#                           beta_scatabund2 = rnorm(1)
#                           ) }


# inits <- function(){ list(beta0 = rnorm(1), 
#                           beta = rnorm(ncov),
#                           beta_bias0 = rnorm(1),
#                           beta_bias = rnorm(ncov_bias),
#                           beta_det0 = rnorm(1),
#                           beta_scatabund0 = rnorm(1),
#                           beta_scatabund1 = exp(rnorm(1)),
#                           # beta_scatabund2 = rnorm(1)
# ) }


# parameters 
params <- c("beta", "beta0", "pred", "q", "pred_abund_scat",
            "pred_encounter")

# MCMC settings

# ni <- 200   ;   nt <- 1   ;   nb <- 100   ;   nc <- 3
ni <- 10000   ;   nt <- 50   ;   nb <- 5000   ;   nc <- 3

# Run JAGS, check convergence and summarize posteriors:
out <- jagsUI::jags(win.data, inits = NULL, params,
                    "explore_scripts/IM_IPPpo_PoissonCount_ttdScat_28May.txt",
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    codaOnly = TRUE, parallel = TRUE)
saveRDS(out, "explore_scripts/out_quadcorrect_28May.rds")
# out <- readRDS("explore_scripts/out_quadcorrect_preliminar_28May.rds")

mask <- stack_all[[1]]
mask[!is.na(mask)] <- 0
mask[complete] <- out$mean$pred
mask[mask>25] <- NA # Removing outliers  
plot(mask)
writeRaster(mask, "explore_scripts/pred_abund_quadcorrect_28May.tif")

# mask[complete] <- out$sd$pred
# plot(mask)


# Plot of predicted abundance of scats for the 100 simulated betas
plot(out$mean$pred_abund_scat ~ win.data$lambda)


# crop_pred <- crop(pred, c(144.5, 145.3, -38.6, -37.9))


# Plot of predicted encounters of scats for the 100 simulated betas
plot(out$mean$pred_encounter ~ win.data$lambda)


# Plot of the mean predicted probability of encountering a scat 
# if present, during the first minute of sampling

plot(out$mean$q[,1])

# Histogram of the mean predicted probability of encountering a scat 
# if present, during the complete five minutes of sampling
hist(apply(out$mean$q, 1, sum))

# Statistics of the mean predicted probability of encountering a scat 
# if present, during the complete five minutes of sampling
summary(apply(out$mean$q, 1, sum))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3132  0.9650  0.9891  0.9699  0.9980  0.9998
