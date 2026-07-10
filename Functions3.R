# Functions3.R

## Load packages

library(tidyverse)
library(dsims)
library(Distance)
library(dssd)
library(dsm)
library(sf)
library(mgcv)
library(qrng) 

# Generates a study region based on either a custom region sf object or 
# upper and lower x and y bounds, allows for units to be specified

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
## using iterative hotspots added to a uniform density

make_hotspot_density <- function(region, 
                                 hotspots, 
                                 x_space = density_grid_spacing) {
  
  # Create basic uniform density
  density_true <- make.density(region = region,
                               x.space = x_space,
                               constant = 1)
  
  for (hotspot in hotspots) {
    density_true <- add.hotspot(object = density_true,
                                centre = hotspot$centre,
                                sigma = hotspot$sigma,
                                amplitude = hotspot$amplitude)
  }
  
  return(density_true)
}

## Hotspots should come in a list like this
# my_hotspots <- list(
# list(centre = c(2000, 8000), sigma = 800, amplitude = 1.2),
# list(centre = c(5000, 2000), sigma = 6000, amplitude = 0.5),
# list(centre = c(8000, 9000), sigma = 1000, amplitude = 2)
# )

## Generates;
## true density surface using hotspot density above,
## population description from that density and region,
## detectability from the scale parameter and truncation distance
## population iteration from the population description, detectability, and region

generate_simulated_pop <- function(region,
                                   density_true,
                                   N = true_N,
                                   scale_param = scale_parameter,
                                   truncation = trunc_dist) {
  
  # Extract the spacing directly from the provided density object
  x_space <- density_true@x.space
  
  density_surface_true <- density_true@density.surface[[1]] |>
    mutate(
      area = as.numeric(sf::st_area(geometry)),
      density = density * N / sum(density * area) # might not need this anymore
      # nevermind, it ensures that the number of animals matches true_N 
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

survey_data <- function(region,
                        realized_population, 
                        # output from (generate_simulated_pop_object)$population
                        angle = design_angle,
                        transect_type = points_or_lines,
                        spacing = design_spacing,
                        truncation = trunc_dist,
                        visits = 1) {
  
  design <- make.design(
    region        = region,
    transect.type = transect_type,  # points or lines
    design        = "systematic",
    spacing       = spacing,
    edge.protocol = "minus",
    design.angle  = angle,
    truncation    = truncation
  )
  
  transects <- generate.transects(design)
  
  # first example of ifelse for points or lines 
  if (transect_type == "line") {
    survey <- new(
      Class = "Survey.LT", 
      population = realized_population,
      transect = transects,
      perp.truncation = truncation
    )
  } else {
    survey <- new(
      Class = "Survey.PT", 
      population = realized_population,
      transect = transects,
      rad.truncation = truncation
    )
  }
  
  observed_survey <- run.survey(survey, region = region)
  
  dist_data <- observed_survey@dist.data
  
  # Re-organize dist_data as obsdata for use in the density surface models
  obsdata <- dist_data |> 
    filter(!is.na(distance)) |> 
    mutate(object = as.character(object),
           Sample.Label = as.character(Sample.Label),
           size = 1) |> 
    select(object, Sample.Label, size, distance)
  
  samplers <- transects@samplers
  
  # calculate effort and extract geometry
  if (transect_type == "line") {
    # Suppress the centroid warning for constant attributes safely
    suppressWarnings({
      sampler_centroids <- sf::st_centroid(samplers)
    })
    sampler_xy <- sf::st_coordinates(sampler_centroids)
    effort <- as.numeric(sf::st_length(samplers))
  } else {
    sampler_xy <- sf::st_coordinates(samplers)
    
    # Ensure visits vector matches the number of samplers if a vector is provided
    if (length(visits) > 1 && length(visits) != nrow(samplers)) {
      stop("Length of 'visits' must be 1 or match the number of point samplers.")
    }
    effort <- visits
  }
  
  segdata <- data.frame(
    x = sampler_xy[, 1],
    y = sampler_xy[, 2],
    Effort = effort,
    Sample.Label = as.character(samplers$transect) 
  )
  
  # Ensure orig_transect exists for the grouping logic later
  if (transect_type == "line") {
    segdata <- segdata |> 
      mutate(orig_transect = as.numeric(Sample.Label))
  }
  
  list(
    design = design,
    transects = transects,
    dist_data = observed_survey@dist.data,
    obsdata = obsdata,
    segdata = segdata
  )
  
}


## fits a ds() model to the above survey data and stores its results
## represents the initial model that would be fit during a real analysis

## calculate with R2, R3, etc, store estimates
## fit with dht() to check agreement with manual fewster 2009 with delta method

fit_ds <- function(region,
                   dist_data, 
                   obsdata, 
                   segdata, 
                   transect_type = points_or_lines, 
                   truncation = trunc_dist) {
  
  # Fit the base model
  m1 <- ds(
    data = dist_data,
    transect = transect_type,
    key = "hn",
    adjustment = NULL,
    truncation = truncation,
    quiet = TRUE
  )
  
  # Extract N_hat, Standard Error, and Sigma
  N_hat <- as.numeric(m1$dht$individuals$N$Estimate)
  N_hat <- if (length(N_hat) > 0) as.integer(round(N_hat)) else NA
  
  # se_ds <- as.numeric(m1$dht$individuals$N$se)
  # se_ds <- if (length(se_ds) > 0) se_ds else NA
  
  sigma_hat <- as.numeric(exp(coef(m1$ddf)$scale["(Intercept)", "estimate"]))
  
  # Extract cv_Pa_sq safely
  cv_p_val <- m1$dht$individuals$summary$cv.p[1]
  if (!is.null(cv_p_val) && length(cv_p_val) > 0) {
    cv_Pa_sq <- (as.numeric(cv_p_val))^2
  } else {
    ddf_sum <- summary(m1$ddf)
    cv_Pa_sq <- (ddf_sum$average.p.se / ddf_sum$average.p)^2
  }
  
  # Calculate global Density
  Region_Area <- as.numeric(sf::st_area(region@region))
  D_hat <- N_hat / Region_Area
  
  # Delta Method Helper
  apply_delta <- function(var_ER, erhat) {
    cv_ER_sq <- var_ER / (erhat^2)
    var_N <- (N_hat^2) * (cv_ER_sq + cv_Pa_sq)
    var_D <- (D_hat^2) * (cv_ER_sq + cv_Pa_sq)
    
    list(
      var_N = var_N, 
      var_D = var_D,
      se_N = sqrt(max(0, var_N)), 
      se_D = sqrt(max(0, var_D))
    )
  }
  
  # calculate observed counts, effort data
  obs_counts <- obsdata |> 
    group_by(Sample.Label) |> 
    summarize(count = sum(size), .groups = "drop")
  
  if (transect_type == "line") {
    effort_data <- segdata |> 
      left_join(obs_counts, by = "Sample.Label") |> 
      mutate(count = replace_na(count, 0)) |> 
      group_by(orig_transect) |> 
      summarize(
        count = sum(count),
        Effort = sum(Effort),
        x = mean(x), 
        .groups = "drop") |> 
      arrange(x) 
  } else {
    effort_data <- segdata |> 
      left_join(obs_counts, by = "Sample.Label") |> 
      mutate(count = replace_na(count, 0))
  }
  
  nspotted <- effort_data$count
  lvec <- effort_data$Effort
  L <- sum(lvec)
  k <- length(lvec) 
  ntot <- sum(nspotted)
  
  erhat_obs <- ntot / L
  
  analytical_variances <- NULL
  
  if (transect_type == "line") {
    
    # Empirical Estimators (R2, R3, S1, S2, O1, O2) 
    
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
    
    # apply delta method to each variance
    delta_R2 <- apply_delta(var.R2, erhat_obs)
    delta_R3 <- apply_delta(var.R3, erhat_obs)
    delta_S1 <- apply_delta(var.S1, erhat_obs)
    delta_S2 <- apply_delta(var.S2, erhat_obs)
    delta_O1 <- apply_delta(var.O1, erhat_obs)
    delta_O2 <- apply_delta(var.O2, erhat_obs)
    
    # Bundle all variances into a comprehensive flat list
    analytical_variances <- list(
      var_ER_R2 = var.R2,
      var_ER_R3 = var.R3,
      var_ER_S1 = var.S1,
      var_ER_S2 = var.S2,
      var_ER_O1 = var.O1,
      var_ER_O2 = var.O2,
      
      se_N_R2 = delta_R2$se_N,
      se_N_R3 = delta_R3$se_N,
      se_N_S1 = delta_S1$se_N,
      se_N_S2 = delta_S2$se_N,
      se_N_O1 = delta_O1$se_N,
      se_N_O2 = delta_O2$se_N,
      
      se_D_R2 = delta_R2$se_D,
      se_D_R3 = delta_R3$se_D,
      se_D_S1 = delta_S1$se_D,
      se_D_S2 = delta_S2$se_D,
      se_D_O1 = delta_O1$se_D,
      se_D_O2 = delta_O2$se_D
    )
  } else {
    
    # p2 <- ds(
    #    data = dist_data,
    #    transect = "point",
    #    key = "hn",
    #    adjustment = NULL,
    #    er_var = "P2",
    #    truncation = truncation,
    #    quiet = TRUE
    # )
    # 
    # se_P2 <- p2$dht$individuals$N$se
    # 
    # p3 <- ds(
    #    data = dist_data,
    #    transect = "point",
    #    key = "hn",
    #    adjustment = NULL,
    #    er_var = "P3",
    #    truncation = truncation,
    #    quiet = TRUE
    # )
    # 
    # se_P3 <- p3$dht$individuals$N$se
    
    nvec <- effort_data$count
    tvec <- effort_data$Effort
    
    k <- length(tvec)
    T <- sum(tvec)
    
    erhat_P2 <- mean(nvec / tvec)
    erhat_P3 <- sum(nvec) / T
    
    var.P2 <- (1 / (k * (k - 1))) *
      sum((nvec / tvec - erhat_P2)^2)
    
    var.P3 <- (1 / (T * (k - 1))) *
      sum(tvec * (nvec / tvec - erhat_P3)^2)
    
    # apply delta method to each variance
    delta_P2 <- apply_delta(var.P2, erhat_obs)
    delta_P3 <- apply_delta(var.P3, erhat_obs)
    
    # Bundle all variances into a comprehensive flat list
    analytical_variances <- list(
      var_ER_P2 = var.P2,
      var_ER_P3 = var.P3,
      
      se_N_P2 = delta_P2$se_N,
      se_N_P3 = delta_P3$se_N,
      
      se_D_P2 = delta_P2$se_D,
      se_D_P3 = delta_P3$se_D
    )
  }
  
  list(
    ds_model = m1,
    N_hat = N_hat,
    D_hat = D_hat,
    # se_ds = se_ds,
    sigma_hat = sigma_hat,
    cv_Pa_sq = cv_Pa_sq,
    analytical_variances = analytical_variances
  )
}


## fitting dsm() and generating estimated surface 
fit_dsm_and_surface <- function(obsdata, 
                                segdata, 
                                ds_model, 
                                region, 
                                N_hat,
                                transect_type = points_or_lines, 
                                x_space = density_grid_spacing, 
                                y_space = density_grid_spacing) {
  
  if (nrow(obsdata) == 0) return(NULL)
  
  # 1. Fit DSM based on geometry
  if (transect_type == "line") {
    dsm1 <- dsm(count ~ s(x), 
                ddf.obj = ds_model, 
                segment.data = segdata,
                observation.data = obsdata, 
                family = quasipoisson(link = "log"),
                method = "REML", 
                convert.units = 1)
  } else {
    dsm1 <- dsm(count ~ s(x, y), 
                ddf.obj = ds_model, 
                segment.data = segdata,
                observation.data = obsdata, 
                family = quasipoisson(link = "log"),
                method = "REML", 
                convert.units = 1)
  }
  
  # 2. Create prediction grid
  prediction_grid <- make.density(region = region, x.space = x_space, y.space = y_space, constant = 1)
  pred_grid <- prediction_grid@density.surface[[1]] |> 
    mutate(area = as.numeric(sf::st_area(geometry)))
  
  pred_data <- sf::st_drop_geometry(pred_grid)
  
  # 3. Predict the surface and apply SAFE calculations
  predictions <- predict(dsm1, newdata = pred_data, off.set = pred_grid$area, type = "response")
  
  pred_grid <- pred_grid |> 
    mutate(
      N_hat_pred = as.numeric(predictions),
      N_hat_pred = replace_na(N_hat_pred, 0),
      # prevents N_hat_pred / 0 from creating NaNs
      density = if_else(area > 0, N_hat_pred / area, 0), 
      density = pmax(density, .Machine$double.eps),
      # final barrier to catch stray NAs so they don't break the bootstrap prob vectors
      density = if_else(is.na(density) | is.nan(density), .Machine$double.eps, density) 
    )
  
  # 4. Generate dsims S4 Objects
  est.density <- make.density(
    region = region, 
    x.space = x_space, 
    y.space = y_space, 
    density.surface = pred_grid
  )
  
  pop.desc <- make.population.description(
    region = region, 
    density = est.density, 
    N = N_hat,          # Uses the N_hat passed into the function
    fixed.N = TRUE
  )
  
  # 5. Build the final dsm_surface object (Handling 1D vs 2D shapes)  
  if (transect_type == "line") {
    dsm_surface <- pred_grid |> 
      group_by(strata, x) |> 
      summarize(
        N_hat_pred = sum(N_hat_pred, na.rm = TRUE), # sum over transect segments vertically
        area = sum(area, na.rm = TRUE),
        geometry = sf::st_union(geometry),
        .groups = "drop"
      ) |>
      mutate(
        # Same safe division applied to lines after merging geometries
        density = if_else(area > 0, N_hat_pred / area, 0),
        density = pmax(density, .Machine$double.eps),
        density = if_else(is.na(density) | is.nan(density), .Machine$double.eps, density),
        y = 0
      ) |> 
      select(strata, density, x, y, N_hat_pred, area, geometry)
    
  } else {
    dsm_surface <- pred_grid |> 
      select(strata, density, x, y, N_hat_pred, area, geometry)
  }
  
  # 6. Return comprehensive list needed for downstream simulation
  list(
    dsm = dsm1,
    obsdata = obsdata, 
    segdata = segdata,
    density_surface = dsm_surface,
    density = est.density,
    population_description = pop.desc,
    pred_grid = pred_grid
  )
}


striplet_boxlet <- function(region, 
                            segdata, 
                            dsm_model, 
                            dsm_surface, 
                            ds_results, 
                            transect_type = points_or_lines, 
                            spacing = design_spacing, 
                            truncation = trunc_dist) {
  
  if (is.null(dsm_model)) return(NULL)
  
  # Unpack needed values from the DS fit results
  m1 <- ds_results$ds_model
  N_hat <- ds_results$N_hat
  D_hat <- ds_results$D_hat
  cv_Pa_sq <- ds_results$cv_Pa_sq
  sigma_hat <- ds_results$sigma_hat
  
  # Local Delta Method
  apply_delta <- function(var_ER, erhat) {
    cv_ER_sq <- var_ER / (erhat^2)
    var_N <- (N_hat^2) * (cv_ER_sq + cv_Pa_sq)
    var_D <- (D_hat^2) * (cv_ER_sq + cv_Pa_sq)
    
    list(se_N = sqrt(max(0, var_N)), se_D = sqrt(max(0, var_D)))
  }
  
  spatial_variances <- list()
  
  if (transect_type == "line") {
    
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
      
      line_strings <- lapply(lines_grid, function(x) {
        sf::st_linestring(matrix(c(x, bbox["ymin"], x, bbox["ymax"]), ncol = 2, byrow = TRUE))
      })
      lines_sf <- sf::st_sfc(line_strings, crs = sf::st_crs(region@region))
      clipped_lines <- suppressWarnings(sf::st_intersection(lines_sf, region@region))
      
      Lbvec[b_idx] <- as.numeric(sum(sf::st_length(clipped_lines)))
      min_dist <- sapply(midvec, function(m) min(abs(m - lines_grid)))
      
      g_b <- exp(- (min_dist^2) / (2 * sigma_hat^2))
      g_b[min_dist > truncation] <- 0 
      
      Abvec[b_idx] <- sum(muvec * g_b, na.rm = TRUE)
    }
    
    # var_ER_striplet <- mean((Abvec + (Abvec^2) * (1 - 1/musum)) / (Lbvec^2)) - (mean(Abvec / Lbvec))^2
    # erhat_striplet <- mean(Abvec / Lbvec, na.rm = TRUE)
    
    var_ER_striplet <- mean((Abvec + (Abvec^2) * (1 - 1/musum)) / (Lbvec^2), na.rm = TRUE) - (mean(Abvec / Lbvec, na.rm = TRUE))^2
    erhat_striplet <- mean(Abvec / Lbvec, na.rm = TRUE)
    
    delta_striplet <- apply_delta(var_ER_striplet, erhat_striplet)
    spatial_variances$se_N_striplet <- delta_striplet$se_N
    spatial_variances$se_D_striplet <- delta_striplet$se_D
    
  } else if (transect_type == "point") {
    
    grid_cells <- 50 
    region_sf <- region@region
    region_bbox <- sf::st_bbox(region_sf)
    
    x_range <- as.numeric(region_bbox["xmax"] - region_bbox["xmin"])
    y_range <- as.numeric(region_bbox["ymax"] - region_bbox["ymin"])
    max_range <- max(x_range, y_range)
    square_dim <- max_range / grid_cells
    
    boxlets <- sf::st_make_grid(region_sf, cellsize = c(square_dim, square_dim), square = TRUE)
    boxlets_sf <- sf::st_sf(geometry = boxlets) |> dplyr::mutate(box_id = dplyr::row_number())
    
    boxlet_cents <- sf::st_centroid(boxlets_sf)
    inside_region <- sf::st_intersects(boxlet_cents, region_sf, sparse = FALSE)[, 1]
    boxlets_sf <- boxlets_sf[inside_region, ]
    
    boxlets_sf <- boxlets_sf |> dplyr::mutate(area = as.numeric(sf::st_area(geometry)))
    boxlet_coords <- sf::st_coordinates(sf::st_centroid(boxlets_sf))
    boxlets_sf$x <- boxlet_coords[, "X"]
    boxlets_sf$y <- boxlet_coords[, "Y"]
    
    pred_data <- sf::st_drop_geometry(boxlets_sf)
    pred_N <- predict(dsm_model, newdata = pred_data, off.set = pred_data$area, type = "response")
    p_j <- pred_N / sum(pred_N, na.rm = TRUE)
    
    p_j[is.na(p_j)] <- 0 
    p_j[p_j < 0] <- 0 
    
    fitted_probs <- m1$ddf$fitted
    if (is.null(fitted_probs) || length(fitted_probs) == 0) {
      P_a <- summary(m1$ddf)$average.p
    } else {
      P_a <- length(fitted_probs) / sum(1 / fitted_probs)
    }
    
    shift_res <- spacing / 5 
    b_x <- seq(0, spacing - shift_res, by = shift_res)
    b_y <- seq(0, spacing - shift_res, by = shift_res)
    shifts <- expand.grid(x = b_x, y = b_y)
    B <- nrow(shifts) 
    
    A_b <- numeric(B)
    
    base_x <- seq(region_bbox["xmin"] - spacing, region_bbox["xmax"] + spacing, by = spacing)
    base_y <- seq(region_bbox["ymin"] - spacing, region_bbox["ymax"] + spacing, by = spacing)
    base_grid <- expand.grid(X = base_x, Y = base_y)
    
    for (i in 1:B) {
      shifted_grid <- base_grid
      shifted_grid$X <- shifted_grid$X + shifts$x[i]
      shifted_grid$Y <- shifted_grid$Y + shifts$y[i]
      
      shifted_points <- sf::st_as_sf(shifted_grid, coords = c("X", "Y"))
      shifted_buffers <- sf::st_buffer(shifted_points, dist = truncation)
      shifted_survey_area <- sf::st_union(shifted_buffers)
      
      active_boxlets <- sf::st_intersects(boxlet_cents[inside_region], shifted_survey_area, sparse = FALSE)[, 1]
      
      Q_b <- sum(p_j[active_boxlets], na.rm = TRUE) * P_a
      A_b[i] <- N_hat * Q_b
    }
    
    mean_A <- mean(A_b)
    var_n <- (1/B) * sum(A_b + (A_b^2) * (1 - 1/N_hat)) - mean_A^2
    var_n_boxlet <- max(0, var_n) 
    
    k_points <- nrow(segdata)
    erhat_boxlet <- mean_A / k_points
    var_ER_boxlet <- var_n_boxlet / (k_points^2)
    
    delta_boxlet <- apply_delta(var_ER_boxlet, erhat_boxlet)
    spatial_variances$se_N_boxlet <- delta_boxlet$se_N
    spatial_variances$se_D_boxlet <- delta_boxlet$se_D
  }
  
  return(spatial_variances)
}




## uses the population from fit_ds, sigma hat from the ds() model,
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


## uses the population from get_fit, sigma hat from the ds() model,
## and the survey design to generate an estimated N_hat 
## Run this multiple times to estimate a SE by bootstraps

get_bootstrap_invcdf <- function(region,
                                 population_description,
                                 sigma_hat,
                                 transect_type = points_or_lines,
                                 reps = bootstrap_reps,
                                 angle = design_angle,
                                 spacing = design_spacing,
                                 truncation = trunc_dist) {
  
  # Pull the density surface from the population description
  density_surface <- population_description@density@density.surface[[1]]
  N_total <- as.integer(round(population_description@N))
  
  # Calculate expected N for each cell based on density * area
  density_surface <- density_surface |>
    mutate(N_expected = density * area)
  
  # Sum over Y to get the PDF, then "integrate" to get the CDF over X
  marg_x <- density_surface |>
    group_by(x) |>
    summarize(N_x = sum(N_expected, na.rm = TRUE), .groups = "drop") |>
    mutate(prob_x = N_x / sum(N_x), cdf_x = cumsum(prob_x))
  
  # Split and calculate Conditional CDFs for Y given X
  cond_y_list <- density_surface |>
    group_by(x) |>
    mutate(prob_y = N_expected / sum(N_expected, na.rm = TRUE), cdf_y = cumsum(prob_y)) |>
    split(~x)
  
  # Create a baseline uniform density to act as a placeholder for a valid S4 initialization
  dummy_density <- dsims::make.density(region = region, x.space = spacing, constant = 1)
  dummy_pop_desc <- dsims::make.population.description(region = region, density = dummy_density, N = N_total, fixed.N = TRUE)
  detect <- make.detectability(key.function = "hn", scale.param = sigma_hat, truncation = truncation)
  
  # Create the survey design 
  design <- make.design(
    region        = region,
    transect.type = transect_type, 
    design        = "systematic",
    spacing       = spacing,
    edge.protocol = "minus",
    design.angle  = angle,
    truncation    = truncation
  )
  
  N_hat_results <- rep(NA, reps)
  
  # Bootstrap loop
  for (b in 1:reps) {
    
    # place animals
    # Use "Owen" scrambling to randomly permute the QMC sequence on every loop iteration
    U <- sobol(n = N_total, d = 2, randomize = "Owen", seed = b)
    
    animal_x <- numeric(N_total)
    animal_y <- numeric(N_total)
    
    # Perform Inverse CDF Sampling
    for (i in 1:N_total) {
      x_idx <- which(marg_x$cdf_x >= U[i, 1])[1]
      exact_x <- marg_x$x[x_idx]
      animal_x[i] <- exact_x
      
      y_dist <- cond_y_list[[as.character(exact_x)]]
      y_idx <- which(y_dist$cdf_y >= U[i, 2])[1]
      animal_y[i] <- y_dist$y[y_idx]
    }
    
    # build s4 object for the generated population
    realized_population <- dsims::generate.population(object = dummy_pop_desc, 
                                                      detectability = detect, 
                                                      region = region)
    
    # Overwrite
    native_pop <- realized_population@population
    native_pop$x <- animal_x
    native_pop$y <- animal_y
    realized_population@population <- as.data.frame(native_pop)
    
    # run survey
    transects <- generate.transects(design)
    
    if (transect_type == "point") {
      survey <- new("Survey.PT", population = realized_population, transect = transects, rad.truncation = truncation)
    } else {
      survey <- new("Survey.LT", population = realized_population, transect = transects, perp.truncation = truncation)
    }
    
    survey_run <- suppressWarnings(run.survey(survey, region = region))
    obs_data <- survey_run@dist.data
    
    # fit model
    if (nrow(obs_data) > 0) {
      tryCatch({
        capture.output(
          m2 <- ds(data = obs_data, 
                   transect = transect_type, key = "hn", 
                   adjustment = NULL, 
                   truncation = truncation, 
                   quiet = TRUE)
        )
        N_hat_results[b] <- as.numeric(m2$dht$individuals$N$Estimate)
        
      }, error = function(e) {
        N_hat_results[b] <- NA
      })
    }
  }
  
  data.frame(replicate = seq_len(reps), N_hat = N_hat_results)
}


## uses the population from get_fit, sigma hat from the ds() model,
## and the survey design to generate an estimated N_hat 
## Run this multiple times to estimate a SE by bootstraps

get_bootstrap_disc_density <- function(region,
                                       population_description,
                                       sigma_hat,
                                       transect_type = points_or_lines,
                                       reps = bootstrap_reps,
                                       angle = design_angle,
                                       spacing = design_spacing,
                                       truncation = trunc_dist) {
  
  # Pull the density surface from the population description
  density_surface <- population_description@density@density.surface[[1]]
  N_total <- as.integer(round(population_description@N))
  
  # Calculate expected N for each cell based on density * area
  density_surface <- density_surface |>
    mutate(N_expected = density * area)
  
  # Calculate the probability of an animal falling into each specific cell
  prob_vec <- density_surface$N_expected / sum(density_surface$N_expected, na.rm = TRUE)
  cell_dimensions <- sqrt(density_surface$area)
  
  # Create a baseline uniform density to act as a placeholder for a valid S4 initialization
  dummy_density <- dsims::make.density(region = region, x.space = spacing, constant = 1)
  dummy_pop_desc <- dsims::make.population.description(
    region = region,
    density = dummy_density,
    N = N_total,
    fixed.N = TRUE
  )
  
  detect <- make.detectability(key.function = "hn", scale.param = sigma_hat, truncation = truncation)
  
  # Create the survey design (This dictates how the grid shifts randomly)
  design <- make.design(
    region        = region,
    transect.type = transect_type, 
    design        = "systematic",
    spacing       = spacing,
    edge.protocol = "minus",
    design.angle  = angle,
    truncation    = truncation
  )
  
  # Empty vector to store variance results
  N_hat_results <- rep(NA, reps)
  
  # Bootstrap loop: Evaluate BOTH spatial layout variance and process variance
  for (b in 1:reps) {
    
    # generate animal placement in the boxlets
    cell_counts <- as.vector(rmultinom(n = 1, size = N_total, prob = prob_vec))
    
    base_x <- rep(density_surface$x, cell_counts)
    base_y <- rep(density_surface$y, cell_counts)
    base_sizes <- rep(cell_dimensions, cell_counts)
    
    # Scatter the animals uniformly within the boundaries of their assigned boxlets
    animal_x <- base_x + runif(N_total, min = -base_sizes/2, max = base_sizes/2)
    animal_y <- base_y + runif(N_total, min = -base_sizes/2, max = base_sizes/2)
    
    # build population based on those placements
    realized_population <- dsims::generate.population(
      object = dummy_pop_desc,
      detectability = detect,
      region = region
    )
    
    # Overwrite the uniform random X and Y coordinates with your new randomized points
    native_pop <- realized_population@population
    native_pop$x <- animal_x
    native_pop$y <- animal_y
    realized_population@population <- as.data.frame(native_pop)
    
    # run surveys
    transects <- generate.transects(design)
    
    if (transect_type == "point") {
      survey <- new("Survey.PT", 
                    population = realized_population, 
                    transect = transects, 
                    rad.truncation = truncation)
    } else {
      survey <- new("Survey.LT", 
                    population = realized_population, 
                    transect = transects, 
                    perp.truncation = truncation)
    }
    
    survey_run <- suppressWarnings(run.survey(survey, region = region))
    obs_data <- survey_run@dist.data
    
    # fit ds model to extract n_hat value
    if (nrow(obs_data) > 0) {
      tryCatch({
        
        capture.output(
          m2 <- ds(data = obs_data, 
                   transect = transect_type, 
                   key = "hn", 
                   adjustment = NULL, 
                   truncation = truncation, 
                   quiet = TRUE)
        )
        
        N_hat_results[b] <- as.numeric(m2$dht$individuals$N$Estimate)
        
      }, error = function(e) {
        
        N_hat_results[b] <- NA
        
      })
    }
  }
  
  data.frame(replicate = seq_len(reps), N_hat = N_hat_results)
}

## EOF