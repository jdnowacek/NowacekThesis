# Functions.R

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
  density_true <- dsims::make.density(region = region,
                                      x.space = x_space,
                                      constant = 1)
  
  for (hotspot in hotspots) {
    density_true <- dsims::add.hotspot(object = density_true,
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

## fits a ds() model to the above survey data and stores its results
## represents the initial model that would be fit during a real analysis
fit_ds_model <- function(dist_data,
                         transect_type = points_or_lines,
                         truncation = trunc_dist) {
  
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
  
  list(
    detection_model = m1,
    N_hat = N_hat,
    se_m1 = se_ds,
    sigma_hat = sigma_hat
  )
}

## data wrangling for, fitting, and generating estimated surface from dsm()
fit_dsm_model <- function(dist_data,
                          transects,
                          region,
                          detection_model,
                          N_hat,
                          transect_type = points_or_lines,
                          truncation = trunc_dist,
                          x_space = density_grid_spacing,
                          y_space = density_grid_spacing) {
  
  m1 <- detection_model
  
  # data wrangling for dsm()
  
  samplers <- transects@samplers
  
  if (transect_type == "line") {
    ## Line Transects:
    # split transects into segments (using grid spacing as a starting value)
    # segment_length <- y_space/3
    segment_length <- trunc_dist # maybe twice the truncation distance?
    
    segdata_list <- list()
    
    obsdata <- dist_data |> 
      filter(!is.na(distance))
    
    new_labels <- character(nrow(obsdata))
    
    # Loop through each transect to segment it
    for (i in seq_len(nrow(samplers))) {
      geom <- samplers$geometry[i]
      tr_label <- as.character(samplers$transect[i])
      
      # Calculate total line length (dropping units for dsm)
      L <- as.numeric(sf::st_length(geom))
      
      # Determine number of segments for this transect
      num_segs <- max(1, round(L / segment_length))
      eff <- L / num_segs
      
      # Find segment centroids using proportional spacing
      # (e.g., for 2 segments, sample at 0.25 and 0.75 along the line)
      fractions <- (seq_len(num_segs) - 0.5) / num_segs
      cents <- sf::st_line_sample(geom, sample = fractions)
      cents_xy <- sf::st_coordinates(cents)
      
      # Create segdata entries for these new segments
      temp_seg <- data.frame(
        Sample.Label = paste0(tr_label, "_", seq_len(num_segs)),
        Effort = eff,
        x = cents_xy[, "X"],
        y = cents_xy[, "Y"],
        orig_transect = tr_label,
        stringsAsFactors = FALSE
      )
      segdata_list[[i]] <- temp_seg
      
      # Map observations on this transect to their closest segment centroid
      obs_idx <- which(as.character(obsdata$Sample.Label) == tr_label)
      
      if (length(obs_idx) > 0) {
        for (j in obs_idx) {
          obs_x <- obsdata$x[j]
          obs_y <- obsdata$y[j]
          
          # Euclidean distance from the animal to all segments on this line
          dists <- sqrt((temp_seg$x - obs_x)^2 + (temp_seg$y - obs_y)^2)
          best_seg <- which.min(dists)
          
          new_labels[j] <- temp_seg$Sample.Label[best_seg]
        }
      }
    }
    
    # Combine segmented data and arrange for variance estimators
    segdata <- bind_rows(segdata_list) |>
      arrange(x, y) |> 
      select(Sample.Label, Effort, x, y, orig_transect)
    
    # Update obsdata with the new granular Segment Labels
    obsdata$Sample.Label <- new_labels
    obsdata <- obsdata |>
      mutate(object = as.character(object), size = 1) |>
      select(object, Sample.Label, size, distance)
    
  } else {
    ## Point transects:
    # Points are already discrete, no geometry splitting needed
    sampler_xy <- sf::st_coordinates(samplers)
    
    segdata <- samplers |>
      sf::st_drop_geometry() |>
      mutate(
        Sample.Label = as.character(transect),
        Effort = 1, # Effort is always 1 for point transects in dsm
        x = sampler_xy[, "X"],
        y = sampler_xy[, "Y"]
      ) |>
      arrange(x, y) |>
      select(Sample.Label, Effort, x, y)
    
    obsdata <- dist_data |>
      filter(!is.na(distance)) |>
      mutate(
        object = as.character(object),
        Sample.Label = as.character(Sample.Label),
        size = 1
      ) |>
      select(object, Sample.Label, size, distance)
  }
  
  # Check if there are valid counts to prevent dsm() crash
  if (nrow(obsdata) == 0) {
    return(list(
      dsm = NULL,
      obsdata = obsdata, 
      segdata = segdata,
      density_surface = NULL,
      density = NULL,
      population_description = NULL,
      pred_grid = NULL
    ))
  }
  
  
  ## Fitting the dsm() models 
  
  # Because the lines have now been spatially segmented, 
  # BOTH points and lines can universally use the 2D s(x, y) spatial smooth!
  dsm1 <- dsm(
    count ~ s(x, y), 
    ddf.obj = m1,
    segment.data = segdata,
    observation.data = obsdata,
    family = quasipoisson(link = "log"),
    method = "REML",
    convert.units = 1
  )
  
  
  ## prediction from dsm() 
  
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
    mutate(N_hat_pred = as.numeric(predictions),
           density = pmax(N_hat_pred / area, .Machine$double.eps))
  
  
  
  ## density, pop from predicted surface 
  
  est.density <- make.density(
    region = region, 
    x.space = x_space, 
    y.space = y_space, 
    density.surface = pred_grid
  )
  
  pop.desc <- make.population.description(
    region = region, 
    density = est.density, 
    N = N_hat, 
    fixed.N = TRUE
  )
  
  
  
  ## building a dsm surface object from preds  
  
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
      select(strata, density, x, y, N_hat_pred, area, geometry)
  }
  
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



## Calculate the variance estimators

## Calculate the variance estimators
calculate_variance_estimators <- function(obsdata,
                                          segdata,
                                          region,
                                          detection_model,
                                          dsm_model,
                                          dsm_surface,
                                          pred_grid,
                                          N_hat,
                                          sigma_hat,
                                          transect_type = points_or_lines,
                                          spacing = design_spacing,
                                          truncation = trunc_dist) {
  
  # If the simulated survey yielded no observations, the models won't fit
  # Return NULL safely instead of breaking.
  if (is.null(dsm_model)) {
    return(NULL)
  }
  
  m1 <- detection_model
  dsm1 <- dsm_model
  
  ##
  ## Delta Method Prep (Shared by both Lines and Points)
  ##
  
  # Safely extract the squared CV of the detection probability 
  cv_p_val <- m1$dht$individuals$summary$cv.p[1]
  if (!is.null(cv_p_val)) {
    cv_Pa_sq <- (as.numeric(cv_p_val))^2
  } else {
    ddf_sum <- summary(m1$ddf)
    cv_Pa_sq <- (ddf_sum$average.p.se / ddf_sum$average.p)^2
  }
  
  # Calculate global Density
  Region_Area <- sum(pred_grid$area, na.rm = TRUE)
  D_hat <- N_hat / Region_Area
  
  # Delta Method Combiner Function
  apply_delta <- function(var_ER, erhat) {
    cv_ER_sq <- var_ER / (erhat^2)
    var_N <- (N_hat^2) * (cv_ER_sq + cv_Pa_sq)
    var_D <- (D_hat^2) * (cv_ER_sq + cv_Pa_sq)
    
    # Convert Variances back to Standard Errors for comparison
    list(
      var_N = var_N, 
      var_D = var_D,
      se_N = sqrt(var_N), 
      se_D = sqrt(var_D)
    )
  }
  
  ##
  ## Stage 1: Construct spatial support (muvec)
  ##
  
  # Both striplets and boxlets use the predicted abundance from the DSM surface.
  # For lines, the surface was grouped by x (forming vertical strips).
  # For points, the surface retains distinct x, y coordinates (forming boxlets).
  muvec <- dsm_surface$N_hat_pred 
  midvec_x <- dsm_surface$x
  midvec_y <- dsm_surface$y
  musum <- sum(muvec, na.rm = TRUE)
  
  bbox <- sf::st_bbox(region@region)
  x_min <- as.numeric(bbox["xmin"])
  x_max <- as.numeric(bbox["xmax"])
  y_min <- as.numeric(bbox["ymin"])
  y_max <- as.numeric(bbox["ymax"])
  
  ##
  ## Stage 2: Construct the Detection Matrix (gmat)
  ##
  
  if (transect_type == "line") {
    
    B_shifts <- 50 
    bvec <- seq(0, spacing, length.out = B_shifts)
    B <- length(bvec)
    
    gmat <- matrix(0, nrow = length(muvec), ncol = B)
    effort_b <- numeric(B)
    
    for (b_idx in 1:B) {
      b_val <- bvec[b_idx]
      lines_grid <- seq(x_min + b_val, x_max, by = spacing)
      
      # Exact line length calculation for effort inside region
      line_strings <- lapply(lines_grid, function(x) {
        sf::st_linestring(matrix(c(x, y_min, x, y_max), ncol = 2, byrow = TRUE))
      })
      lines_sf <- sf::st_sfc(line_strings, crs = sf::st_crs(region@region))
      clipped_lines <- suppressWarnings(sf::st_intersection(lines_sf, region@region))
      effort_b[b_idx] <- as.numeric(sum(sf::st_length(clipped_lines)))
      
      # Distance to nearest line for each striplet
      min_dist <- sapply(midvec_x, function(m) min(abs(m - lines_grid)))
      
      g_b <- exp(- (min_dist^2) / (2 * sigma_hat^2))
      g_b[min_dist > truncation] <- 0 
      
      gmat[, b_idx] <- g_b
    }
    
  } else if (transect_type == "point") {
    
    shift_res <- spacing / 5
    b_x <- seq(0, spacing - shift_res, by = shift_res)
    b_y <- seq(0, spacing - shift_res, by = shift_res)
    shifts <- expand.grid(x = b_x, y = b_y)
    B <- nrow(shifts)
    
    gmat <- matrix(0, nrow = length(muvec), ncol = B)
    effort_b <- numeric(B)
    
    # Establish baseline infinite point grid
    base_x <- seq(x_min - spacing, x_max + spacing, by = spacing)
    base_y <- seq(y_min - spacing, y_max + spacing, by = spacing)
    base_grid <- expand.grid(X = base_x, Y = base_y)
    
    for (i in 1:B) {
      shifted_x <- base_grid$X + shifts$x[i]
      shifted_y <- base_grid$Y + shifts$y[i]
      
      # Count valid effort (total points falling perfectly inside the region boundaries)
      shifted_points <- sf::st_as_sf(data.frame(X = shifted_x, Y = shifted_y), 
                                     coords = c("X", "Y"), crs = sf::st_crs(region@region))
      intersects <- sf::st_intersects(shifted_points, region@region, sparse = FALSE)
      effort_b[i] <- sum(intersects)
      
      # Distance to nearest point for each boxlet 
      # Optimized: fast mathematical Euclidean matching over discrete grid logic
      nearest_x <- round((midvec_x - shifts$x[i]) / spacing) * spacing + shifts$x[i]
      nearest_y <- round((midvec_y - shifts$y[i]) / spacing) * spacing + shifts$y[i]
      min_dist <- sqrt((midvec_x - nearest_x)^2 + (midvec_y - nearest_y)^2)
      
      # Determine continuous probability inclusion 
      g_b <- exp(- (min_dist^2) / (2 * sigma_hat^2))
      g_b[min_dist > truncation] <- 0 
      
      gmat[, i] <- g_b
    }
  }
  
  ##
  ## Stage 3: Fewster Variance Engine (Identical for both Boxlets & Striplets)
  ##
  
  # Protect against theoretical shift permutations completely outside the region
  valid_b <- effort_b > 0
  effort_b <- effort_b[valid_b]
  gmat <- gmat[, valid_b, drop = FALSE]
  
  # A_b is the expected count for each grid alignment. Vectorized for extreme speed.
  Abvec <- as.numeric(crossprod(gmat, muvec))
  
  # Master analytical variance equation 
  var_ER_analytical <- mean((Abvec + (Abvec^2) * (1 - 1/musum)) / (effort_b^2)) - (mean(Abvec / effort_b))^2
  var_ER_analytical <- max(0, var_ER_analytical)
  
  erhat_analytical <- mean(Abvec / effort_b, na.rm = TRUE)
  delta_analytical <- apply_delta(var_ER_analytical, erhat_analytical)
  
  
  ##
  ## Stage 4: Empirical Estimators & Output Structuring
  ##
  
  if (transect_type == "line") {
    
    obs_counts <- obsdata |> 
      group_by(Sample.Label) |> 
      summarize(count = sum(size), .groups = "drop")
    
    line_data <- segdata |> 
      left_join(obs_counts, by = "Sample.Label") |> 
      mutate(count = replace_na(count, 0)) |> 
      group_by(orig_transect) |> 
      summarize(
        count = sum(count),
        Effort = sum(Effort),
        x = mean(x), 
        .groups = "drop") |> 
      arrange(x) 
    
    nspotted <- line_data$count
    lvec <- line_data$Effort
    L <- sum(lvec)
    k <- length(lvec) 
    ntot <- sum(nspotted)
    
    ## Empirical Estimators (R2, R3, S1, S2, O1, O2)
    var.R2 <- (k * sum(lvec^2 * (nspotted/lvec - ntot/L)^2)) / (L^2 * (k - 1))
    var.R3 <- 1 / (L * (k - 1)) * sum(lvec * (nspotted/lvec - ntot/L)^2)
    
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
    
    lvec.1 <- lvec[-k]; lvec.2 <- lvec[-1]
    nvec.1 <- nspotted[-k]; nvec.2 <- nspotted[-1]
    ervec.1 <- nvec.1/lvec.1; ervec.2 <- nvec.2/lvec.2
    
    var.O1 <- k / (2 * L^2 * (k - 1)) * sum((nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2)
    var.O2 <- (2 * k) / (L^2 * (k - 1)) * sum(((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1 - ervec.2)^2)
    
    erhat_obs <- ntot / L
    
    delta_R2 <- apply_delta(var.R2, erhat_obs)
    delta_R3 <- apply_delta(var.R3, erhat_obs)
    delta_S1 <- apply_delta(var.S1, erhat_obs)
    delta_S2 <- apply_delta(var.S2, erhat_obs)
    delta_O1 <- apply_delta(var.O1, erhat_obs)
    delta_O2 <- apply_delta(var.O2, erhat_obs)
    
    analytical_variances <- list(
      var_ER_R2 = var.R2,
      var_ER_R3 = var.R3,
      var_ER_S1 = var.S1,
      var_ER_S2 = var.S2,
      var_ER_O1 = var.O1,
      var_ER_O2 = var.O2,
      var_ER_striplet = var_ER_analytical,
      
      se_N_R2 = delta_R2$se_N,
      se_N_R3 = delta_R3$se_N,
      se_N_S1 = delta_S1$se_N,
      se_N_S2 = delta_S2$se_N,
      se_N_O1 = delta_O1$se_N,
      se_N_O2 = delta_O2$se_N,
      se_N_striplet = delta_analytical$se_N,
      
      se_D_R2 = delta_R2$se_D,
      se_D_R3 = delta_R3$se_D,
      se_D_S1 = delta_S1$se_D,
      se_D_S2 = delta_S2$se_D,
      se_D_O1 = delta_O1$se_D,
      se_D_O2 = delta_O2$se_D,
      se_D_striplet = delta_analytical$se_D
    )
    
  } else if (transect_type == "point") {
    
    # Points lack the traditional empirical estimators, returning purely the boxlet 
    analytical_variances <- list(
      var_ER_boxlet = var_ER_analytical,
      se_N_boxlet = delta_analytical$se_N,
      se_D_boxlet = delta_analytical$se_D
    )
  }
  
  return(analytical_variances)
}

# calculate_variance_estimators <- function(obsdata,
#                                           segdata,
#                                           region,
#                                           detection_model,
#                                           dsm_model,
#                                           dsm_surface,
#                                           pred_grid,
#                                           N_hat,
#                                           sigma_hat,
#                                           transect_type = points_or_lines,
#                                           spacing = design_spacing,
#                                           truncation = trunc_dist) {
#   
#   # If the simulated survey yielded no observations, the models won't fit
#   # Return NULL safely instead of breaking.
#   if (is.null(dsm_model)) {
#     return(NULL)
#   }
#   
#   m1 <- detection_model
#   dsm1 <- dsm_model
#   
#   # empty object for variance estimate storage
#   analytical_variances <- NULL
#   
#   
#   ## Delta Method Prep (Shared by both Lines and Points)
#   
#   # Safely extract the squared CV of the detection probability 
#   # Checks the dht summary first, falls back to the ddf summary if NULL
#   cv_p_val <- m1$dht$individuals$summary$cv.p[1]
#   if (!is.null(cv_p_val)) {
#     cv_Pa_sq <- (as.numeric(cv_p_val))^2
#   } else {
#     ddf_sum <- summary(m1$ddf)
#     cv_Pa_sq <- (ddf_sum$average.p.se / ddf_sum$average.p)^2
#   }
#   
#   # Calculate global Density
#   Region_Area <- sum(pred_grid$area, na.rm = TRUE)
#   D_hat <- N_hat / Region_Area
#   
#   # Delta Method Combiner Function
#   # NOTE: Group-size variance is intentionally omitted from this equation 
#   # because the simulation forces `size = 1` (individuals) in the obsdata.
#   apply_delta <- function(var_ER, erhat) {
#     cv_ER_sq <- var_ER / (erhat^2)
#     var_N <- (N_hat^2) * (cv_ER_sq + cv_Pa_sq)
#     var_D <- (D_hat^2) * (cv_ER_sq + cv_Pa_sq)
#     
#     # Convert Variances back to Standard Errors for comparison
#     list(
#       var_N = var_N, 
#       var_D = var_D,
#       se_N = sqrt(var_N), 
#       se_D = sqrt(var_D)
#     )
#   }
#   
#   
#   ##
#   ## Line Transect Estimators
#   ##
#   if (transect_type == "line") {
#     
#     obs_counts <- obsdata |> 
#       group_by(Sample.Label) |> 
#       summarize(count = sum(size), 
#                 .groups = "drop")
#     
#     line_data <- segdata |> 
#       left_join(obs_counts, by = "Sample.Label") |> 
#       mutate(count = replace_na(count, 0)) |> 
#       group_by(orig_transect) |> 
#       summarize(
#         count = sum(count),
#         Effort = sum(Effort),
#         x = mean(x), # Get the average X-coordinate of the whole line
#         .groups = "drop") |> 
#       arrange(x) # Re-sorts the whole lines spatially left-to-right
#     
#     nspotted <- line_data$count
#     lvec <- line_data$Effort
#     L <- sum(lvec)
#     k <- length(lvec) # 'k' is back to being the true number of lines!
#     ntot <- sum(nspotted)
#     
#     ## Empirical Estimators (R2, R3, S1, S2, O1, O2)
#     
#     ## R2, R3
#     var.R2 <- (k * sum(lvec^2 * (nspotted/lvec - ntot/L)^2)) / (L^2 * (k - 1))
#     var.R3 <- 1 / (L * (k - 1)) * sum(lvec * (nspotted/lvec - ntot/L)^2)
#     
#     ## Stratified (S1, S2)
#     H_strat <- floor(k/2)
#     k.h <- rep(2, H_strat)
#     if(k %% 2 > 0) k.h[H_strat] <- 3
#     end.strat <- cumsum(k.h)
#     begin.strat <- cumsum(k.h) - k.h + 1
#     
#     sum.S1 <- 0; sum.S2 <- 0
#     for(h in 1:H_strat) {
#       n.strat <- nspotted[begin.strat[h]:end.strat[h]]
#       l.strat <- lvec[begin.strat[h]:end.strat[h]]
#       nbar.strat <- mean(n.strat)
#       lbar.strat <- mean(l.strat)
#       
#       inner.S1 <- sum((n.strat - nbar.strat - (ntot/L) * (l.strat - lbar.strat))^2)
#       sum.S1 <- sum.S1 + k.h[h] / (k.h[h] - 1) * inner.S1
#       
#       L.strat <- sum(l.strat)
#       var.strat.S2 <- k.h[h] / (L.strat^2 * (k.h[h] - 1)) * sum(l.strat^2 * (n.strat/l.strat - nbar.strat/lbar.strat)^2)
#       sum.S2 <- sum.S2 + L.strat^2 * var.strat.S2
#     }
#     var.S1 <- sum.S1 / L^2
#     var.S2 <- sum.S2 / L^2
#     
#     ## Overlapping (O1, O2)
#     lvec.1 <- lvec[-k]; lvec.2 <- lvec[-1]
#     nvec.1 <- nspotted[-k]; nvec.2 <- nspotted[-1]
#     ervec.1 <- nvec.1/lvec.1; ervec.2 <- nvec.2/lvec.2
#     
#     var.O1 <- k / (2 * L^2 * (k - 1)) * sum((nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2)
#     var.O2 <- (2 * k) / (L^2 * (k - 1)) * sum(((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1 - ervec.2)^2)
#     
#     
#     ## Striplet Variance
#     muvec <- dsm_surface$N_hat_pred 
#     midvec <- dsm_surface$x
#     musum <- sum(muvec, na.rm = TRUE)
#     
#     bbox <- sf::st_bbox(region@region)
#     y_length <- as.numeric(bbox["ymax"] - bbox["ymin"])
#     x_min <- as.numeric(bbox["xmin"])
#     x_max <- as.numeric(bbox["xmax"])
#     
#     B <- 50 
#     bvec <- seq(0, spacing, length.out = B)
#     Lbvec <- rep(0, B)
#     Abvec <- rep(0, B)
#     
#     for (b_idx in 1:B) {
#       b_val <- bvec[b_idx]
#       lines_grid <- seq(x_min + b_val, x_max, by = spacing)
#       
#       # Build theoretical vertical sf lines extending across the bounding box
#       line_strings <- lapply(lines_grid, function(x) {
#         sf::st_linestring(matrix(c(x, bbox["ymin"], x, bbox["ymax"]), ncol = 2, byrow = TRUE))
#       })
#       lines_sf <- sf::st_sfc(line_strings, crs = sf::st_crs(region@region))
#       
#       # Clip the theoretical lines to the exact boundary of the irregular region
#       clipped_lines <- suppressWarnings(sf::st_intersection(lines_sf, region@region))
#       
#       # Calculate exact line length inside the polygon bounds
#       Lbvec[b_idx] <- as.numeric(sum(sf::st_length(clipped_lines)))
#       
#       min_dist <- sapply(midvec, function(m) min(abs(m - lines_grid)))
#       
#       g_b <- exp(- (min_dist^2) / (2 * sigma_hat^2))
#       g_b[min_dist > truncation] <- 0 
#       
#       Abvec[b_idx] <- sum(muvec * g_b, na.rm = TRUE)
#     }
#     
#     var_ER_striplet <- mean((Abvec + (Abvec^2) * (1 - 1/musum)) / (Lbvec^2)) - (mean(Abvec / Lbvec))^2
#     
#     # Define Encounter Rates (Observed vs. Model)
#     # Traditional estimators use the raw observed encounter rate:
#     erhat_obs <- ntot / L
#     
#     # The Striplet estimator should theoretically use the mean corresponding
#     # directly to the striplet model rather than the raw counts:
#     erhat_striplet <- mean(Abvec / Lbvec, na.rm = TRUE)
#     
#     # Apply to all empirical estimators using the observed encounter rate
#     delta_R2 <- apply_delta(var.R2, erhat_obs)
#     delta_R3 <- apply_delta(var.R3, erhat_obs)
#     delta_S1 <- apply_delta(var.S1, erhat_obs)
#     delta_S2 <- apply_delta(var.S2, erhat_obs)
#     delta_O1 <- apply_delta(var.O1, erhat_obs)
#     delta_O2 <- apply_delta(var.O2, erhat_obs)
#     
#     # Apply to the striplet estimator using the strictly aligned striplet mean
#     delta_striplet <- apply_delta(var_ER_striplet, erhat_striplet)
#     
#     # Bundle all variances into a comprehensive flat list
#     analytical_variances <- list(
#       # Encounter Rate Variances (Raw)
#       var_ER_R2 = var.R2,
#       var_ER_R3 = var.R3,
#       var_ER_S1 = var.S1,
#       var_ER_S2 = var.S2,
#       var_ER_O1 = var.O1,
#       var_ER_O2 = var.O2,
#       var_ER_striplet = var_ER_striplet,
#       
#       # Final Abundance Standard Errors se(N) -- USE THESE FOR COMPARISON
#       se_N_R2 = delta_R2$se_N,
#       se_N_R3 = delta_R3$se_N,
#       se_N_S1 = delta_S1$se_N,
#       se_N_S2 = delta_S2$se_N,
#       se_N_O1 = delta_O1$se_N,
#       se_N_O2 = delta_O2$se_N,
#       se_N_striplet = delta_striplet$se_N,
#       
#       # Final Density Standard Errors se(D)
#       se_D_R2 = delta_R2$se_D,
#       se_D_R3 = delta_R3$se_D,
#       se_D_S1 = delta_S1$se_D,
#       se_D_S2 = delta_S2$se_D,
#       se_D_O1 = delta_O1$se_D,
#       se_D_O2 = delta_O2$se_D,
#       se_D_striplet = delta_striplet$se_D
#     )
#     
#     
#     ##
#     ## Point Transect Estimators 
#     ##
#   } else if (transect_type == "point") {
#     
#     ### Boxlet Variance ----
#     grid_res <- 100
#     
#     # Tessellate the Region (The Boxlets)
#     region_sf <- region@region
#     region_bbox <- sf::st_bbox(region_sf)
#     
#     # Create a fine grid of 2D 'boxlets' across the region
#     boxlets <- sf::st_make_grid(region_sf, cellsize = c(grid_res, grid_res), square = TRUE)
#     boxlets_sf <- sf::st_sf(geometry = boxlets) |> 
#       dplyr::mutate(
#         box_id = dplyr::row_number(),
#         area = as.numeric(sf::st_area(geometry))
#       )
#     
#     boxlet_centroids <- sf::st_centroid(boxlets_sf)
#     boxlet_coords <- sf::st_coordinates(boxlet_centroids)
#     boxlets_sf$x <- boxlet_coords[, "X"]
#     boxlets_sf$y <- boxlet_coords[, "Y"]
#     
#     ## Estimate Boxlet Probabilities (p_j)
#     
#     # Predict the spatial trend using the 2D GAM s(x,y)
#     pred_data <- sf::st_drop_geometry(boxlets_sf)
#     pred_N <- predict(dsm1, newdata = pred_data, off.set = pred_data$area, type = "response")
#     
#     # Normalize predictions into multinomial probabilities
#     p_j <- pred_N / sum(pred_N, na.rm = TRUE)
#     
#     # Safely compute Average Detection Probability (P_a)
#     # The standard math is: P_a = n / N_covered = length(fitted) / sum(1 / fitted)
#     fitted_probs <- m1$ddf$fitted
#     if (is.null(fitted_probs) || length(fitted_probs) == 0) {
#       P_a <- summary(m1$ddf)$average.p
#     } else {
#       P_a <- length(fitted_probs) / sum(1 / fitted_probs)
#     }
#     
#     ## Simulate the Grid Shifts (b)
#     
#     # Define the uniform sampling frame for the start point 'b'
#     shift_res <- spacing / 5 # The number of shifts per axis (adjust for speed vs. precision)
#     b_x <- seq(0, spacing - shift_res, by = shift_res)
#     b_y <- seq(0, spacing - shift_res, by = shift_res)
#     shifts <- expand.grid(x = b_x, y = b_y)
#     B <- nrow(shifts) # Total number of possible grid alignments
#     
#     # Pre-allocate A(b) vector (expected counts per alignment)
#     A_b <- numeric(B)
#     
#     # Create a theoretical base grid of points spanning far outside the region 
#     # to ensure full coverage when shifted
#     base_x <- seq(region_bbox["xmin"] - spacing, region_bbox["xmax"] + spacing, by = spacing)
#     base_y <- seq(region_bbox["ymin"] - spacing, region_bbox["ymax"] + spacing, by = spacing)
#     base_grid <- expand.grid(X = base_x, Y = base_y)
#     
#     ## Calculate Q(b) for Each Shift
#     
#     for (i in 1:B) {
#       # Shift the theoretical point grid by b_x and b_y
#       shifted_grid <- base_grid
#       shifted_grid$X <- shifted_grid$X + shifts$x[i]
#       shifted_grid$Y <- shifted_grid$Y + shifts$y[i]
#       
#       # Convert to spatial points
#       shifted_points <- sf::st_as_sf(shifted_grid, coords = c("X", "Y"))
#       
#       # Create circular buffers (the point transect radii)
#       shifted_buffers <- sf::st_buffer(shifted_points, dist = truncation)
#       shifted_survey_area <- sf::st_union(shifted_buffers)
#       
#       # Determine which boxlet centroids fall inside the circular samplers
#       # This is the spatial equivalent of Fewster's 'gbmat' indicator matrix
#       active_boxlets <- sf::st_intersects(boxlet_centroids, shifted_survey_area, sparse = FALSE)[, 1]
#       
#       # Sum the probabilities of active boxlets and adjust for detection
#       Q_b <- sum(p_j[active_boxlets], na.rm = TRUE) * P_a
#       
#       # Expected count for this specific alignment
#       A_b[i] <- N_hat * Q_b
#     }
#     
#     ## Apply the Multinomial Equation
#     
#     # Mean expected count across all grid shifts
#     mean_A <- mean(A_b)
#     
#     # 2D Analytical Variance Equation (from Fewster's mgcv.boxlet.func)
#     var_n <- (1/B) * sum(A_b + (A_b^2) * (1 - 1/N_hat)) - mean_A^2
#     var_n_boxlet <- max(0, var_n) # Ensure non-negative variance before rooting
#     
#     
#     # The Delta Method Integration for Boxlets
#     
#     # For point transects, Encounter Rate (ER) is count / number of points (n/k)
#     k_points <- nrow(segdata)
#     erhat_boxlet <- mean_A / k_points
#     var_ER_boxlet <- var_n_boxlet / (k_points^2)
#     
#     delta_boxlet <- apply_delta(var_ER_boxlet, erhat_boxlet)
#     
#     # Bundle all variances into a comprehensive flat list
#     analytical_variances <- list(
#       var_ER_boxlet = var_ER_boxlet,
#       se_N_boxlet = delta_boxlet$se_N,
#       se_D_boxlet = delta_boxlet$se_D
#     )
#   }
#   
#   return(analytical_variances)
# }

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
    
    # --- 1. GENERATE NEW RANDOM ANIMAL PLACEMENT ---
    cell_counts <- as.vector(rmultinom(n = 1, size = N_total, prob = prob_vec))
    
    base_x <- rep(density_surface$x, cell_counts)
    base_y <- rep(density_surface$y, cell_counts)
    base_sizes <- rep(cell_dimensions, cell_counts)
    
    # Scatter the animals uniformly within the strict boundaries of their assigned boxlets
    animal_x <- base_x + runif(N_total, min = -base_sizes/2, max = base_sizes/2)
    animal_y <- base_y + runif(N_total, min = -base_sizes/2, max = base_sizes/2)
    
    # --- 2. BUILD S4 POPULATION FOR THIS ITERATION ---
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
    
    # --- 3. RUN SURVEY ---
    transects <- generate.transects(design)
    
    if (transect_type == "point") {
      survey <- new("Survey.PT", population = realized_population, transect = transects, rad.truncation = truncation)
    } else {
      survey <- new("Survey.LT", population = realized_population, transect = transects, perp.truncation = truncation)
    }
    
    survey_run <- suppressWarnings(run.survey(survey, region = region))
    obs_data <- survey_run@dist.data
    
    # --- 4. FIT MODEL ---
    if (nrow(obs_data) > 0) {
      tryCatch({
        capture.output(
          m2 <- ds(data = obs_data, transect = transect_type, key = "hn", adjustment = NULL, truncation = truncation, quiet = TRUE)
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