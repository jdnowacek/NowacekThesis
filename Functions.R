#Functions.R

## Load packages

library(tidyverse)
library(dsims)
library(Distance)
library(dssd)
library(dsm)
library(sf)
library(mgcv)
library(qrng) 

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
  
  samplers <- transects@samplers
  
  if (transect_type == "line") {
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
      select(Sample.Label, Effort, x, y, orig_transect) # <--- ADD orig_transect here
    
    # Update obsdata with the new granular Segment Labels
    obsdata$Sample.Label <- new_labels
    obsdata <- obsdata |>
      mutate(object = as.character(object), size = 1) |>
      select(object, Sample.Label, size, distance)
    
  } else {
    # --- POINT TRANSECT APPROACH ---
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
      detection_model = m1, dsm = NULL, N_hat = N_hat,
      se_m1 = se_ds, sigma_hat = sigma_hat,
      obsdata = obsdata, segdata = segdata
    ))
  }
  
  
  ## dsm() models ----- 
  
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
    mutate(N_hat_pred = as.numeric(predictions),
           density = pmax(N_hat_pred / area, .Machine$double.eps))
  
  
  
  ## density, pop from surface -----
  
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
      select(strata, density, x, y, N_hat_pred, area, geometry)
  }
  
  
  
  
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
        mutate(count = replace_na(count, 0)) |> 
        group_by(orig_transect) |> 
        summarize(
          count = sum(count),
          Effort = sum(Effort),
          x = mean(x), # Get the average X-coordinate of the whole line
          .groups = "drop") |> 
        arrange(x) # CRITICAL: Re-sorts the whole lines spatially left-to-right
      
      nspotted <- line_data$count
      lvec <- line_data$Effort
      L <- sum(lvec)
      k <- length(lvec) # 'k' is back to being the true number of lines!
      ntot <- sum(nspotted)
    
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
    
    # extract the squared CV of the detection probability directly from distance model
    cv_Pa_sq <- (as.numeric(m1$dht$individuals$summary$cv.p[1]))^2
    
    # Calculate global Density
    Region_Area <- sum(pred_grid$area, na.rm = TRUE)
    D_hat <- N_hat / Region_Area
    
    # 2. Define Encounter Rates (Observed vs. Model)
    # Traditional estimators use the raw observed encounter rate:
    erhat_obs <- ntot / L
    
    # The Striplet estimator should theoretically use the mean corresponding
    # directly to the striplet model rather than the raw counts:
    erhat_striplet <- mean(Abvec / Lbvec, na.rm = TRUE)
    
    # Delta Method Combiner Function
    # NOTE: Group-size variance is intentionally omitted from this equation 
    # because the simulation forces `size = 1` (individuals) in the obsdata.
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
    
    # Apply to all empirical estimators using the observed encounter rate
    delta_R2 <- apply_delta(var.R2, erhat_obs)
    delta_R3 <- apply_delta(var.R3, erhat_obs)
    delta_S1 <- apply_delta(var.S1, erhat_obs)
    delta_S2 <- apply_delta(var.S2, erhat_obs)
    delta_O1 <- apply_delta(var.O1, erhat_obs)
    delta_O2 <- apply_delta(var.O2, erhat_obs)
    
    # Apply to the striplet estimator using the strictly aligned striplet mean
    delta_striplet <- apply_delta(var_ER_striplet, erhat_striplet)
    
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

# get_bootstrap_invcdf <- function(region,
#                               population_description,
#                               sigma_hat,
#                               transect_type = points_or_lines,
#                               reps = bootstrap_reps,
#                               angle = design_angle,
#                               spacing = design_spacing,
#                               truncation = trunc_dist) {
#   
#   # Pull the density surface from the population description
#   density_surface <- population_description@density@density.surface[[1]]
#   N_total <- as.integer(round(population_description@N))
#   
#   # Calculate expected N for each cell based on density * area
#   density_surface <- density_surface |>
#     mutate(N_expected = density * area)
#   
#   # Sum over Y to get the PDF, then "integrate" to get the CDF over X
#   marg_x <- density_surface |>
#     group_by(x) |>
#     summarize(N_x = sum(N_expected, na.rm = TRUE), .groups = "drop") |>
#     mutate(
#       prob_x = N_x / sum(N_x),
#       cdf_x = cumsum(prob_x)
#     )
#   
#   # Split and calculate Conditional CDFs for Y given X
#   cond_y_list <- density_surface |>
#     group_by(x) |>
#     mutate(
#       prob_y = N_expected / sum(N_expected, na.rm = TRUE),
#       cdf_y = cumsum(prob_y)
#     ) |>
#     split(~x)
#   
#   # Uniform samples from QMC (Sobol)
#   # Generate N_total Sobol points in 2 dimensions
#   U <- sobol(n = N_total, d = 2, randomize = "none")
#   
#   animal_x <- numeric(N_total)
#   animal_y <- numeric(N_total)
#   
#   # Perform Inverse CDF Sampling
#   for (i in 1:N_total) {
#     # Match U1 to the X CDF
#     x_idx <- which(marg_x$cdf_x >= U[i, 1])[1]
#     exact_x <- marg_x$x[x_idx]
#     animal_x[i] <- exact_x
#     
#     # Match U2 to the Y CDF specific to that X column
#     y_dist <- cond_y_list[[as.character(exact_x)]]
#     y_idx <- which(y_dist$cdf_y >= U[i, 2])[1]
#     animal_y[i] <- y_dist$y[y_idx]
#   }
#   
#   # Create a baseline uniform density to act as a placeholder for a valid S4 initialization
#   dummy_density <- dsims::make.density(region = region, x.space = spacing, constant = 1)
#   
#   # Use native dsims tool to build a perfectly configured population description
#   dummy_pop_desc <- dsims::make.population.description(
#     region = region,
#     density = dummy_density,
#     N = N_total,
#     fixed.N = TRUE
#   )
#   
#   detect <- make.detectability(
#     key.function = "hn",
#     scale.param = sigma_hat,
#     truncation = truncation
#   )
#   
#   # Generate a valid population structure natively
#   realized_population <- dsims::generate.population(
#     object = dummy_pop_desc,
#     detectability = detect,
#     region = region
#   )
#   
#   # Extract the natively generated data frame (contains scale.param, individual, etc.)
#   native_pop <- realized_population@population
#   
#   # Overwrite the uniform random X and Y coordinates with your perfect QMC points
#   native_pop$x <- animal_x
#   native_pop$y <- animal_y
#   
#   # Ensure it is a strict base data.frame to satisfy S4 validation
#   realized_population@population <- as.data.frame(native_pop)
#   
#   # Create the survey design (This dictates how the grid shifts randomly)
#   design <- make.design(
#     region        = region,
#     transect.type = transect_type, 
#     design        = "systematic",
#     spacing       = spacing,
#     edge.protocol = "minus",
#     design.angle  = angle,
#     truncation    = truncation
#   )
#   
#   # Empty vector to store variance results
#   N_hat_results <- rep(NA, reps)
#   
#   # Bootstrap loop: Evaluate spatial layout variance ONLY
#   for (b in 1:reps) {
#     
#     # Generate shifting transects for this iteration
#     transects <- generate.transects(design)
#     
#     # Construct the Survey S4 Object based on transect type
#     if (transect_type == "point") {
#       survey <- new("Survey.PT", 
#                     population = realized_population, 
#                     transect = transects, 
#                     rad.truncation = truncation)
#     } else {
#       survey <- new("Survey.LT", 
#                     population = realized_population, 
#                     transect = transects, 
#                     perp.truncation = truncation)
#     }
#     
#     # Run survey to compute Euclidean distances and Detection Probabilities
#     survey_run <- suppressWarnings(run.survey(survey, region = region))
#     obs_data <- survey_run@dist.data
#     
#     # Fit the ds() model safely
#     if (nrow(obs_data) > 0) {
#       tryCatch({
#         capture.output(
#           m2 <- ds(data = obs_data, 
#                    transect = transect_type, 
#                    key = "hn", 
#                    adjustment = NULL, 
#                    truncation = truncation, 
#                    quiet = TRUE)
#         )
#         
#         N_hat_results[b] <- as.numeric(m2$dht$individuals$N$Estimate)
#         
#       }, error = function(e) {
#         # If the model fails to fit (e.g., extremely low counts), gracefully pass NA
#         N_hat_results[b] <- NA
#       })
#     }
#   }
#   
#   # Return data frame identical to standard get_bootstrap output
#   data.frame(
#     replicate = seq_len(reps),
#     N_hat = N_hat_results
#   )
# }

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

# get_bootstrap_disc_density <- function(region,
#                                  population_description,
#                                  sigma_hat,
#                                  transect_type = points_or_lines,
#                                  reps = bootstrap_reps,
#                                  angle = design_angle,
#                                  spacing = design_spacing,
#                                  truncation = trunc_dist) {
#   
#   # Pull the density surface from the population description
#   density_surface <- population_description@density@density.surface[[1]]
#   N_total <- as.integer(round(population_description@N))
#   
#   # Calculate expected N for each cell based on density * area
#   density_surface <- density_surface |>
#     mutate(N_expected = density * area)
#   
#   # Calculate the probability of an animal falling into each specific cell
#   prob_vec <- density_surface$N_expected / sum(density_surface$N_expected, na.rm = TRUE)
#   
#   # Draw exactly N_total animals, distributed across the cells based on those probabilities
#   # rmultinom returns a matrix, so we use as.vector() to flatten it
#   cell_counts <- as.vector(rmultinom(n = 1, size = N_total, prob = prob_vec))
#   
#   # Calculate the width/height of each square cell (Area = L^2, so L = sqrt(Area))
#   cell_dimensions <- sqrt(density_surface$area)
#   
#   # Expand the centroids and dimensions so there is one value per simulated animal
#   # If a cell got 3 animals, its x, y, and size are repeated 3 times in these vectors
#   base_x <- rep(density_surface$x, cell_counts)
#   base_y <- rep(density_surface$y, cell_counts)
#   base_sizes <- rep(cell_dimensions, cell_counts)
#   
#   # Scatter the animals uniformly within the strict boundaries of their assigned boxlets
#   # We offset the centroids by a random number between -(width/2) and +(width/2)
#   animal_x <- base_x + runif(N_total, min = -base_sizes/2, max = base_sizes/2)
#   animal_y <- base_y + runif(N_total, min = -base_sizes/2, max = base_sizes/2)
#   
#   # Create a baseline uniform density to act as a placeholder for a valid S4 initialization
#   dummy_density <- dsims::make.density(region = region, x.space = spacing, constant = 1)
#   
#   # Use native dsims tool to build a perfectly configured population description
#   dummy_pop_desc <- dsims::make.population.description(
#     region = region,
#     density = dummy_density,
#     N = N_total,
#     fixed.N = TRUE
#   )
#   
#   detect <- make.detectability(
#     key.function = "hn",
#     scale.param = sigma_hat,
#     truncation = truncation
#   )
#   
#   # Generate a valid population structure natively
#   realized_population <- dsims::generate.population(
#     object = dummy_pop_desc,
#     detectability = detect,
#     region = region
#   )
#   
#   # Extract the natively generated data frame (contains scale.param, individual, etc.)
#   native_pop <- realized_population@population
#   
#   # Overwrite the uniform random X and Y coordinates with your perfect QMC points
#   native_pop$x <- animal_x
#   native_pop$y <- animal_y
#   
#   # Ensure it is a strict base data.frame to satisfy S4 validation
#   realized_population@population <- as.data.frame(native_pop)
#   
#   # Create the survey design (This dictates how the grid shifts randomly)
#   design <- make.design(
#     region        = region,
#     transect.type = transect_type, 
#     design        = "systematic",
#     spacing       = spacing,
#     edge.protocol = "minus",
#     design.angle  = angle,
#     truncation    = truncation
#   )
#   
#   # Empty vector to store variance results
#   N_hat_results <- rep(NA, reps)
#   
#   # Bootstrap loop: Evaluate spatial layout variance ONLY
#   for (b in 1:reps) {
#     
#     # Generate shifting transects for this iteration
#     transects <- generate.transects(design)
#     
#     # Construct the Survey S4 Object based on transect type
#     if (transect_type == "point") {
#       survey <- new("Survey.PT", 
#                     population = realized_population, 
#                     transect = transects, 
#                     rad.truncation = truncation)
#     } else {
#       survey <- new("Survey.LT", 
#                     population = realized_population, 
#                     transect = transects, 
#                     perp.truncation = truncation)
#     }
#     
#     # Run survey to compute Euclidean distances and Detection Probabilities
#     survey_run <- suppressWarnings(run.survey(survey, region = region))
#     obs_data <- survey_run@dist.data
#     
#     # Fit the ds() model safely
#     if (nrow(obs_data) > 0) {
#       tryCatch({
#         capture.output(
#           m2 <- ds(data = obs_data, 
#                    transect = transect_type, 
#                    key = "hn", 
#                    adjustment = NULL, 
#                    truncation = truncation, 
#                    quiet = TRUE)
#         )
#         
#         N_hat_results[b] <- as.numeric(m2$dht$individuals$N$Estimate)
#         
#       }, error = function(e) {
#         # If the model fails to fit (e.g., extremely low counts), gracefully pass NA
#         N_hat_results[b] <- NA
#       })
#     }
#   }
#   
#   # Return data frame identical to standard get_bootstrap output
#   data.frame(
#     replicate = seq_len(reps),
#     N_hat = N_hat_results
#   )
# }

## EOF