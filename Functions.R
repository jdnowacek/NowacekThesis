## Load packages

library(tidyverse)
library(dsims)
library(Distance)
library(dssd)
library(dsm)
library(sf)
library(mgcv)

# Functions

## Generates a study region based on upper and lower x and y bounds 
## and allows for units to be specified

# make_region <- function(lower_x_bound,
#                         upper_x_bound,
#                         lower_y_bound,
#                         upper_y_bound,
#                         units) {
#   
#   # create survey sf object
#   simsquare <- sf::st_bbox(c(xmin = lower_x_bound, xmax = upper_x_bound, 
#                              ymin = lower_y_bound, ymax = upper_y_bound)) |> 
#     sf::st_as_sfc() |> 
#     sf::st_as_sf()
#   
#   # make region with dssd
#   region <- make.region(region.name = "region",
#                               shape = simsquare,
#                               units = units)
#   
#   return(region)
# }

make_region <- function(region_sf = NULL, 
                        lower_x_bound = 0, upper_x_bound = 0, 
                        lower_y_bound = 0, upper_y_bound = 0, 
                        units = "m") { 
  
    # if not custom region, use bounds to create rectangle
    if (is.null(region_sf)) {
    region_sf <- sf::st_bbox(c(xmin = lower_x_bound, xmax = upper_x_bound, 
                               ymin = lower_y_bound, ymax = upper_y_bound)) |> 
      sf::st_as_sfc() |> 
      sf::st_as_sf()
  }
  
  region <- make.region(region.name = "region",
                        shape = region_sf,
                        units = units)
  return(region)
}

## Make theoretical density of animals
## using iterative hotspots added to the densities

make_hotspot_density <- function(region, 
                                 hotspots, 
                                 x_space = density_grid_spacing,
                                 plot_density = FALSE) {
  
  # Create basic uniform density
  density_true <- dsims::make.density(region = region,
                                      x.space = x_space,
                                      constant = 1)
  
  for (hotspot in hotspots) {
    density_true <- dsims::add.hotspot(object = density_true,
                                       centre = hotspot$centre,
                                       sigma = hotspot$sigma,
                                       amplitude = hotspot$amplitude)
    
    # Added visual check
    if (plot_density) {
      # Uses dsims native plot method
      plot(density_true, region) 
    }
    
    return(density_true)
  }
  
  ## Hotspots should come in a list like this
  # my_hotspots <- list(
  # list(centre = c(2000, 8000), sigma = 800, amplitude = 1.2),
  # list(centre = c(5000, 2000), sigma = 6000, amplitude = 0.5),
  # list(centre = c(8000, 9000), sigma = 1000, amplitude = 2)
  # )
  
  return(density_true)
}

## Generates;
## true density surface using hotspot density above,
## population description from that density and region,
## detectability from the scale parameter and truncation distance
## population iteration from the population description, detecability, and region

generate_simulated_pop <- function(region,
                                     N = true_N,
                                     x_space = density_grid_spacing,
                                     hotspots = my_hotspots,
                                     scale_param = scale_parameter,
                                     truncation = trunc_dist) {
  
  density_true <- make_hotspot_density(region = region,
                               hotspots = hotspots,
                               x_space = x_space)
  
  density_surface_true <- density_true@density.surface[[1]] |>
    mutate(
      area = as.numeric(sf::st_area(geometry)),
      density = density * N / sum(density * area)
    ) |>
    select(strata, density, x, y, geometry)
  
  density_true <- make.density(
    region = region,
    x.space = x_space,
    density.surface = density_surface_true
  )
  
  pop.desc_true <- make.population.description(
    region = region,
    density = density_true,
    N = N,
    fixed.N = TRUE
  )
  
  detect_true <- make.detectability(
    key.function = "hn",
    scale.param = scale_param,
    truncation = truncation
  )
  
  pop_true <- generate.population(
    object = pop.desc_true,
    detectability = detect_true,
    region = region
  )
  
  list(
    density = density_true,
    population_description = pop.desc_true,
    population = pop_true,
    detectability = detect_true
  )
}

## Generates survey data from the population specified in generate simulated pop
## Allows for different survey types (points or lines)
## returns survey design object, 
## transects from that design,
## 'dist data' which contains the results of the survey specified

generate_survey_data <- function(region,
                                 realized_population, 
                                 # output from (generate_simulated_pop_object)$population
                                 angle = design_angle,
                                 transect_type = points_or_lines,
                                 spacing = design_spacing,
                                 truncation = trunc_dist) {
  
  design <- make.design(
    region        = region,
    transect.type = transect_type,  # Evaluates the literal variable passed ("point" or "line")
    design        = "systematic",
    spacing       = spacing,
    edge.protocol = "minus",
    design.angle  = angle,
    truncation    = truncation
  )
  
  transects <- generate.transects(design)
  
  # Conditional routing for dsims strict Survey classes 
  if (transect_type == "point") {
    survey <- new(
      Class = "Survey.PT", 
      population = realized_population,
      transect = transects,
      rad.truncation = truncation 
    )
  } else {
    survey <- new(
      Class = "Survey.LT", 
      population = realized_population,
      transect = transects,
      perp.truncation = truncation
    )
  }
  
  observed_survey <- run.survey(survey, region = region)
  
  list(
    design = design,
    transects = transects,
    dist_data = observed_survey@dist.data
  )
}

## Uses distribution data from generate survey data, and design to:
## generate estimates and SEs from distance sampling model
## fit a density surface
## predict from that density surface
## make a density and population iteration from that surface
## Calculate analytical variances from the survey

get_fit <- function(dist_data,
                    transects,
                    region,
                    transect_type = points_or_lines,
                    spacing = design_spacing,
                    truncation = trunc_dist,
                    x_space = density_grid_spacing,
                    y_space = density_grid_spacing) {
  
  ## ds() model ----- 
  
  # Fit the distance sampling model on the 'true' data from the survey above
  m1 <- ds(
    data = dist_data,
    transect = transect_type,
    key = "hn",
    adjustment = NULL,
    truncation = truncation,
    quiet = TRUE
  )
  
  # Store estimates of total count from the ds() model
  N_hat <- as.numeric(m1$dht$individuals$N$Estimate)
  N_hat <- if (length(N_hat) > 0) as.integer(round(N_hat)) else NA
  
  # Safely extract standard errors (fall back to NA if numeric(0))
  se_ds <- as.numeric(m1$dht$individuals$N$se)
  se_ds <- if (length(se_ds) > 0) se_ds else NA
  
  # Store value of the scale parameter estimate from the ds() model
  sigma_hat <- as.numeric(exp(coef(m1$ddf)$scale["(Intercept)", "estimate"]))
  
  ## data wrangling for dsm() ----- 
  
  # re-organizes dist_data as obsdata for use in the density surface models
  obsdata <- dist_data |>
    filter(!is.na(distance)) |>
    mutate(object = as.character(object),
           Sample.Label = as.character(Sample.Label),
           size = 1) |>
    select(object, Sample.Label, size, distance)
  
  # pulls the samplers from the transects object for sf data
  samplers <- transects@samplers
  
  # pulls the geometry of the transects for lines or points
  # calculates effort using line length
  if (transect_type == "line") {
    # Suppress the centroid warning for constant attributes
    # sf::st_agr(samplers) <- "constant" 
    sampler_centroids <- sf::st_centroid(samplers)
    sampler_xy <- sf::st_coordinates(sampler_centroids)
    effort <- as.numeric(sf::st_length(samplers))
  } else {
    sampler_xy <- sf::st_coordinates(samplers)
    effort <- 1
  }
  
  # produces segdata, which is geometry compatible for dsm models
  # arranges lines so that they work for overlapping var. est. methods
  segdata <- samplers |>
    sf::st_drop_geometry() |>
    mutate(
      Sample.Label = as.character(transect),
      Effort = effort,
      x = sampler_xy[, "X"],
      y = sampler_xy[, "Y"]
    ) |>
    select(Sample.Label, Effort, x, y) |> 
    arrange(x) # CRITICAL: Orders segments spatially for overlapping variances
  
  # Check if there are valid counts for the DSM
  if (nrow(obsdata) == 0) {
    return(list(
      detection_model = m1, dsm = NULL, N_hat = N_hat,
      se_m1 = se_ds, sigma_hat = sigma_hat,
      obsdata = obsdata, segdata = segdata, density_surface = NULL,
      density = NULL, population_description = NULL,
      analytical_variances = NULL
    ))
  }
  
  ## dsm() models ----- 
  
  # fits the density surface models
  # different for lines than for points bc we have already grouped the lines
  # to the point that there is no information for the y part of the smooth
  if (transect_type == "line") {
    dsm1 <- dsm(
      count ~ s(x), # 1D smooth for lines
      ddf.obj = m1,
      segment.data = segdata,
      observation.data = obsdata,
      family = quasipoisson(link = "log"),
      method = "REML",
      convert.units = 1
    )
  } else {
    dsm1 <- dsm(
      count ~ s(x, y),     # 2D smooth for points
      ddf.obj = m1,
      segment.data = segdata,
      observation.data = obsdata,
      family = quasipoisson(link = "log"),
      method = "REML",
      convert.units = 1
    )
  }
  
  ## prediction from dsm() ----- 
  
  prediction_grid <- make.density(
    region = region,
    x.space = x_space,
    y.space = y_space,
    constant = 1
  )
  
  pred_grid <- prediction_grid@density.surface[[1]] |>
    mutate(area = as.numeric(sf::st_area(geometry)))
  
  pred_data <- sf::st_drop_geometry(pred_grid)
  
  predictions <- predict(dsm1, newdata = pred_data, off.set = pred_grid$area, type = "response")
  
  pred_grid <- pred_grid |> 
    mutate(N_hat_pred = as.numeric(predictions))
  
  ## dsm surface from preds ----- 
  
  if (transect_type == "line") {
    
    dsm_surface <- pred_grid |> 
      group_by(strata, x) |> 
      summarize(
        N_hat_pred = sum(N_hat_pred, na.rm = TRUE),
        area = sum(area, na.rm = TRUE),
        geometry = sf::st_union(geometry),
        .groups = "drop"
      ) |> 
      mutate(density = pmax(N_hat_pred / area, .Machine$double.eps), y = 0) |> 
      select(strata, density, x, y, N_hat_pred, area, geometry)
    
  } else {
    dsm_surface <- pred_grid |> 
      mutate(density = pmax(N_hat_pred / area, .Machine$double.eps)) |> 
      select(strata, density, x, y, N_hat_pred, area, geometry)
  }
  
  ## density, pop from surface -----
  
  est.density <- make.density(
    region = region, 
    x.space = x_space, 
    y.space = y_space, 
    density.surface = dsm_surface)
  
  pop.desc <- make.population.description(
    region = region, 
    density = est.density, 
    N = N_hat, 
    fixed.N = TRUE)
  
  ## variance estimators ----
  
  # empty object for variance estimate storage
  analytical_variances <- NULL
  
  # sets up a section where variances from striplet estimator,
  # as well as 2009 Fewster variance estimators can be calculated
  # uses fewster code from her email
  if (transect_type == "line") {
    
    obs_counts <- obsdata |> 
      group_by(Sample.Label) |> 
      summarize(count = sum(size), 
                .groups = "drop")
    
    line_data <- segdata |> 
      left_join(obs_counts, by = "Sample.Label") |> 
      mutate(count = replace_na(count, 0))
    
    nspotted <- line_data$count
    lvec <- line_data$Effort
    L <- sum(lvec)
    k <- length(lvec)
    ntot <- sum(nspotted)
    
    ### Empirical Estimators (R2, R3, S1, S2, O1, O2) -----
    
    ## R2, R3
    var.R2 <- (k * sum(lvec^2 * (nspotted/lvec - ntot/L)^2)) / (L^2 * (k - 1))
    var.R3 <- 1 / (L * (k - 1)) * sum(lvec * (nspotted/lvec - ntot/L)^2)
    
    ## Stratified (S1, S2)
    H_strat <- floor(k/2)
    k.h <- rep(2, H_strat)
    if(k %% 2 > 0) k.h[H_strat] <- 3
    end.strat <- cumsum(k.h)
    begin.strat <- cumsum(k.h) - k.h + 1
    
    sum.S1 <- 0; sum.S2 <- 0
    for(h in 1:H_strat) {
      n.strat <- nspotted[begin.strat[h]:end.strat[h]]
      l.strat <- lvec[begin.strat[h]:end.strat[h]]
      nbar.strat <- mean(n.strat)
      lbar.strat <- mean(l.strat)
      
      inner.S1 <- sum((n.strat - nbar.strat - (ntot/L) * (l.strat - lbar.strat))^2)
      sum.S1 <- sum.S1 + k.h[h] / (k.h[h] - 1) * inner.S1
      
      L.strat <- sum(l.strat)
      var.strat.S2 <- k.h[h] / (L.strat^2 * (k.h[h] - 1)) * sum(l.strat^2 * (n.strat/l.strat - nbar.strat/lbar.strat)^2)
      sum.S2 <- sum.S2 + L.strat^2 * var.strat.S2
    }
    var.S1 <- sum.S1 / L^2
    var.S2 <- sum.S2 / L^2
    
    ## Overlapping (O1, O2)
    lvec.1 <- lvec[-k]; lvec.2 <- lvec[-1]
    nvec.1 <- nspotted[-k]; nvec.2 <- nspotted[-1]
    ervec.1 <- nvec.1/lvec.1; ervec.2 <- nvec.2/lvec.2
    
    var.O1 <- k / (2 * L^2 * (k - 1)) * sum((nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2)
    var.O2 <- (2 * k) / (L^2 * (k - 1)) * sum(((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1 - ervec.2)^2)
    
    
    ## Striplet Variance -----
    muvec <- dsm_surface$N_hat_pred 
    midvec <- dsm_surface$x
    musum <- sum(muvec, na.rm = TRUE)
    
    bbox <- sf::st_bbox(region@region)
    y_length <- as.numeric(bbox["ymax"] - bbox["ymin"])
    x_min <- as.numeric(bbox["xmin"])
    x_max <- as.numeric(bbox["xmax"])
    
    B <- 50 
    bvec <- seq(0, spacing, length.out = B)
    Lbvec <- rep(0, B)
    Abvec <- rep(0, B)
    
    for (b_idx in 1:B) {
      b_val <- bvec[b_idx]
      lines_grid <- seq(x_min + b_val, x_max, by = spacing)
      
      
      
      # old version
      # Lbvec[b_idx] <- length(lines_grid) * y_length
      
      # new version
      # Build theoretical vertical sf lines extending across the bounding box
      line_strings <- lapply(lines_grid, function(x) {
        sf::st_linestring(matrix(c(x, bbox["ymin"], x, bbox["ymax"]), ncol = 2, byrow = TRUE))
      })
      lines_sf <- sf::st_sfc(line_strings, crs = sf::st_crs(region@region))
      
      # Clip the theoretical lines to the exact boundary of the irregular region
      clipped_lines <- suppressWarnings(sf::st_intersection(lines_sf, region@region))
      
      # Calculate exact line length inside the polygon bounds
      Lbvec[b_idx] <- as.numeric(sum(sf::st_length(clipped_lines)))
      
      
      
      min_dist <- sapply(midvec, function(m) min(abs(m - lines_grid)))
      
      g_b <- exp(- (min_dist^2) / (2 * sigma_hat^2))
      g_b[min_dist > truncation] <- 0 
      
      Abvec[b_idx] <- sum(muvec * g_b, na.rm = TRUE)
    }
    
    var_ER_striplet <- mean((Abvec + (Abvec^2) * (1 - 1/musum)) / (Lbvec^2)) - (mean(Abvec / Lbvec))^2
    
    ## The Delta Method: Combining var(ER) with var(Detection)
    
    # Safely extract the detection probability variance from the 'ds' model
    cv_N_m1_sq <- (se_ds / N_hat)^2
    cv_ER_m1_sq <- (as.numeric(m1$dht$individuals$summary$cv.ER[1]))^2
    cv_Pa_sq <- max(0, cv_N_m1_sq - cv_ER_m1_sq) # Isolate detection variance
    
    # Calculate global Density
    Region_Area <- sum(pred_grid$area, na.rm = TRUE)
    D_hat <- N_hat / Region_Area
    erhat <- ntot / L
    
    # Delta Method Combiner Function
    apply_delta <- function(var_ER) {
      cv_ER_sq <- var_ER / (erhat^2)
      var_N <- (N_hat^2) * (cv_ER_sq + cv_Pa_sq)
      var_D <- (D_hat^2) * (cv_ER_sq + cv_Pa_sq)
      
      # FIX: Convert Variances back to Standard Errors for comparison
      list(
        var_N = var_N, 
        var_D = var_D,
        se_N = sqrt(var_N), 
        se_D = sqrt(var_D)
      )
    }
    
    # Apply to all estimators
    delta_R2 <- apply_delta(var.R2)
    delta_R3 <- apply_delta(var.R3)
    delta_S1 <- apply_delta(var.S1)
    delta_S2 <- apply_delta(var.S2)
    delta_O1 <- apply_delta(var.O1)
    delta_O2 <- apply_delta(var.O2)
    delta_striplet <- apply_delta(var_ER_striplet)
    
    # Bundle all variances into a comprehensive flat list
    analytical_variances <- list(
      # Encounter Rate Variances (Raw)
      var_ER_R2 = var.R2,
      var_ER_R3 = var.R3,
      var_ER_S1 = var.S1,
      var_ER_S2 = var.S2,
      var_ER_O1 = var.O1,
      var_ER_O2 = var.O2,
      var_ER_striplet = var_ER_striplet,
      
      # Final Abundance Standard Errors se(N) -- USE THESE FOR COMPARISON
      se_N_R2 = delta_R2$se_N,
      se_N_R3 = delta_R3$se_N,
      se_N_S1 = delta_S1$se_N,
      se_N_S2 = delta_S2$se_N,
      se_N_O1 = delta_O1$se_N,
      se_N_O2 = delta_O2$se_N,
      se_N_striplet = delta_striplet$se_N,
      
      # Final Density Standard Errors se(D)
      se_D_R2 = delta_R2$se_D,
      se_D_R3 = delta_R3$se_D,
      se_D_S1 = delta_S1$se_D,
      se_D_S2 = delta_S2$se_D,
      se_D_O1 = delta_O1$se_D,
      se_D_O2 = delta_O2$se_D,
      se_D_striplet = delta_striplet$se_D
    )
  }
  
  # Return final integrated list
  list(
    detection_model = m1, dsm = dsm1, N_hat = N_hat, se_m1 = se_ds, 
    sigma_hat = sigma_hat, obsdata = obsdata, segdata = segdata, 
    density_surface = dsm_surface, density = est.density, 
    population_description = pop.desc, 
    analytical_variances = analytical_variances
  )
}


## uses the population from get_fit, sigma hat from the ds() model,
## and the survey design to generate an estimated N_hat 
## Run this multiple times to estimate a SE by bootstraps

get_bootstrap <- function(region,
                          population_description,
                          sigma_hat,
                          transect_type = points_or_lines,
                          reps = bootstrap_reps,
                          angle = design_angle,
                          spacing = design_spacing,
                          truncation = trunc_dist) {
  
  
  detect <- make.detectability(
    key.function = "hn",
    scale.param = sigma_hat,
    truncation = truncation
  )
  
  
  analyses <- make.ds.analysis(
    dfmodel = ~ 1,
    key = "hn",
    truncation = truncation,
    criteria = "AIC"
  )
  
  design <- make.design(
    region        = region,
    transect.type = transect_type,  # Evaluates the literal variable passed ("point" or "line")
    design        = "systematic",
    spacing       = spacing,
    edge.protocol = "minus",
    design.angle  = angle,
    truncation    = truncation
  )
  
  sim_obj <- make.simulation(
    reps = reps,
    design = design, 
    population.description = population_description,
    detectability = detect,
    ds.analysis = analyses
  )
  
  sim_obj <- run.simulation(sim_obj)
  
  N_hat <- as.numeric(
    sim_obj@results[["individuals"]][["N"]][1, "Estimate", seq_len(reps)]
  )
  
  data.frame(
    replicate = seq_len(reps),
    N_hat = N_hat
  )
}

## EOF