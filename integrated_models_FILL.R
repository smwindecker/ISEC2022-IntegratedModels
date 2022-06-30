#' ---
#' title: "Integrated model - for filling in"
#' author: "Saras Windecker & Nick Golding"
#' date: "20 June, 2022"
#' output: html_document
#' ---
#' 

#' packages
library(raster)

#' covariate for distribution of possums
tree_cover <- readRDS('data/nveg_ras.rds')

# extract center and scale for full region
tree_cover_scaled <- scale(getValues(tree_cover))
scale <- attr(tree_cover_scaled, 'scaled:scale')
center <- attr(tree_cover_scaled, 'scaled:center')

#' PRESENCE-ONLY DATA
po_raw <- read.csv('data/po.csv')
po_ras <- rasterize(po_raw, tree_cover, 
                    fun = 'count', background = 0)
po <- getValues(po_ras)
table(po)

# vis data
plot(tree_cover)
points(SpatialPoints(po_raw), pch = 20)

# area for po data is sqkm
area_ras <- area(tree_cover)
area_po <- getValues(area_ras)

# covariates for po data
tree_cover_po <- as.vector(tree_cover_scaled)

access_ras <- readRDS('data/access_ras.rds')
access <- as.vector(scale(getValues(access_ras)))

#' COUNT DATA
counts <- read.csv('data/counts.csv')
counts_sp <- SpatialPoints(counts[, c('lon', 'lat')])

# area of count data is search radius 50m. convert to sqkm here
area_counts <- (pi * 50 ^ 2) / 1e6

tree_cover_counts <- as.vector(scale(
  extract(tree_cover, 
          counts_sp), 
  center = center,
  scale = scale))

#' DETECTION DATA
scat_det <- read.csv('data/scat_det.csv')
scat_det_sp  <- SpatialPoints(scat_det[, c('lon', 'lat')])

# scat search radius 10m
area_scat_det <- (pi * 10 ^ 2) / 1e6

tree_cover_scat_det <- as.vector(scale(
  extract(tree_cover, 
          scat_det_sp), 
  center = center,
  scale = scale))

#' predictions across whole region
tree_cover_pred <- as.vector(tree_cover_scaled)
access_pred <- access

#' create jags data
win.data <- list(
  
  # presence-only data
  n_po = length(po), 
  po = po,
  tree_cover_po = tree_cover_po, 
  area_po = area_po,
  access = access,
  
  # spotlight count data
  n_counts = length(counts),
  counts = counts$count,
  tree_cover_counts = tree_cover_counts,
  area_counts = area_counts,
  
  # time to detection scat data 
  n_scat_det = length(scat_det),
  tree_cover_scat_det = tree_cover_scat_det,
  area_scat_det = area_scat_det,
  
  # prediction
  n_pred = length(tree_cover_pred),
  tree_cover_pred = tree_cover_pred,
  access_pred = access_pred
)

sink("jags_mod.txt")
cat("

# Define model and write into R working directory

model {

  ## Priors ##

  # Priors for abundance of possums
  alpha ~ dnorm(0, 0.1)

  # Prior for slope of abundance predictor
  beta ~ dnorm(0, 0.1) 

  # Prior for possum spatial sampling bias
  alpha_bias ~ dnorm(0, 0.1)

  # Prior for slope in bias predictors
  beta_bias ~ dnorm(0, 0.1)

  # Prior for scat detection
  alpha_scat_det ~ dnorm(0, 0.1)

  ## Likelihood ##

  # count data of possums
  for (i in 1:n_counts){

  }

  # presence-only possum data
  for (j in 1:n_po){
   
  }

  # scat detection data
  for (k in 1:n_scat_det) { 

  }

  ## Derived quantities ##
  
  for (i in 1:n_pred){

    log(lambda_pred[i]) <- alpha + beta*tree_cover_pred[i]
    log(bias_pred[i]) <- alpha_bias + beta_bias*access_pred[i]

  }

}

", fill = TRUE)
sink()

#' MCMC
# Run JAGS, check convergence and summarize posteriors:
out <- jagsUI::jags(data = win.data, inits = NULL,
                    parameters.to.save = c("beta", "alpha",
                                           "alpha_scat_det",
                                           "alpha_bias", "beta_bias",
                                           "lambda_pred"),
                    model.file = "jags_mod.txt",
                    n.chains = 3, n.thin = 1, n.iter = 5000, n.burnin = 2500,
                    codaOnly = TRUE, parallel = FALSE)
# saveRDS(out, "jags_mod.rds")
# out <- readRDS('jags_mod.rds')

#' prediction
pred_map <- tree_cover
values(pred_map) <- as.vector(out$mean$lambda_pred)
pred_map <- mask(pred_map, tree_cover)
plot(pred_map)
