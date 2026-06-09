get_fit_striplet <- function(
    dist_data, 
    transects, 
    region, 
    # transect_type = "line", # Defaulting to line for this specific function 
    spacing = design_spacing, # Required to know how to shift the grid!
    truncation = trunc_dist, 
    x_space = density_grid_spacing, 
    y_space = density_grid_spacing) {
  
  # 1. Fit the basic detection model
  m1 <- ds(
    data = dist_data, 
    transect = transect_type, 
    key = "hn", 
    adjustment = NULL, 
    truncation = truncation, 
    quiet = TRUE
  )
  
  N_hat <- as.numeric(m1$dht$individuals$N$Estimate)
  N_hat <- if (length(N_hat) > 0) as.integer(round(N_hat)) else NA
  
  se_ds <- as.numeric(m1$dht$individuals$N$se)
  se_ds <- if (length(se_ds) > 0) se_ds else NA
  sigma_hat <- as.numeric(exp(coef(m1$ddf)$scale["(Intercept)", "estimate"]))
  
  # 2. Wrangle Observation Data
  obsdata <- dist_data |> 
    filter(!is.na(distance)) |> 
    mutate(
      object = as.character(object), 
      Sample.Label = as.character(Sample.Label), 
      size = 1
    ) |> 
    select(object, Sample.Label, size, distance)
  
  # 3. Wrangle Segment Data & Centroids
  samplers <- transects@samplers
  
  if (transect_type == "line") {
    sf::st_agr(samplers) <- "constant"
    sampler_centroids <- sf::st_centroid(samplers)
    sampler_xy <- sf::st_coordinates(sampler_centroids)
    effort <- as.numeric(sf::st_length(samplers))
  } else {
    sampler_xy <- sf::st_coordinates(samplers)
    effort <- 1
  }
  
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
  
  # Catch empty surveys
  if (nrow(obsdata) == 0) {
    return(list(
      detection_model = m1, dsm = NULL, N_hat = N_hat,
      se_m1 = se_ds, sigma_hat = sigma_hat, obsdata = obsdata, segdata = segdata,
      density_surface = NULL, density = NULL, population_description = NULL,
      analytical_variances = NULL
    ))
  }
  
  # 4. Fit DSM
  if (transect_type == "line") {
    dsm1 <- dsm(count ~ s(x), ddf.obj = m1, segment.data = segdata, 
                observation.data = obsdata, family = quasipoisson(link = "log"), 
                method = "REML", convert.units = 1)
  } else {
    dsm1 <- dsm(count ~ s(x, y), ddf.obj = m1, segment.data = segdata, 
                observation.data = obsdata, family = quasipoisson(link = "log"), 
                method = "REML", convert.units = 1)
  }
  
  # 5. Prediction Grid
  prediction_grid <- make.density(region = region, x.space = x_space, y.space = y_space, constant = 1)
  pred_grid <- prediction_grid@density.surface[[1]] |> mutate(area = as.numeric(sf::st_area(geometry)))
  pred_data <- sf::st_drop_geometry(pred_grid)
  
  pred_grid$N_hat <- predict(dsm1, newdata = pred_data, off.set = pred_grid$area, type = "response")
  
  # 6. Striplet Generation
  if (transect_type == "line") {
    dsm_surface <- pred_grid |> 
      group_by(strata, x) |> 
      summarize(
        N_hat = sum(N_hat, na.rm = TRUE),
        area = sum(area, na.rm = TRUE),
        geometry = sf::st_union(geometry),
        .groups = "drop"
      ) |> 
      mutate(density = pmax(N_hat / area, .Machine$double.eps), y = NA) |> 
      select(strata, density, x, y, geometry)
  } else {
    dsm_surface <- pred_grid |> mutate(density = pmax(N_hat / area, .Machine$double.eps)) |> 
      select(strata, density, x, y, geometry)
  }
  
  gam.density <- make.density(region = region, x.space = x_space, y.space = y_space, density.surface = dsm_surface)
  pop.desc <- make.population.description(region = region, density = gam.density, N = N_hat, fixed.N = TRUE)
  
  # =========================================================================
  # 7. ANALYTICAL VARIANCE BLOCK (Only runs for Line Transects)
  # =========================================================================
  analytical_variances <- NULL
  
  if (transect_type == "line") {
    
    # 7a. Setup Empirical Variables
    obs_counts <- obsdata |> group_by(Sample.Label) |> summarize(count = sum(size), .groups = "drop")
    line_data <- segdata |> left_join(obs_counts, by = "Sample.Label") |> mutate(count = replace_na(count, 0))
    
    nspotted <- line_data$count
    lvec <- line_data$Effort
    L <- sum(lvec)
    k <- length(lvec)
    ntot <- sum(nspotted)
    
    # 7b. Empirical Estimators (R2, R3, S1, S2, O1, O2)
    var.R2 <- (k * sum(lvec^2 * (nspotted/lvec - ntot/L)^2)) / (L^2 * (k - 1))
    var.R3 <- 1 / (L * (k - 1)) * sum(lvec * (nspotted/lvec - ntot/L)^2)
    
    # Stratified (S1, S2)
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
    
    # Overlapping (O1, O2)
    lvec.1 <- lvec[-k]; lvec.2 <- lvec[-1]
    nvec.1 <- nspotted[-k]; nvec.2 <- nspotted[-1]
    ervec.1 <- nvec.1/lvec.1; ervec.2 <- nvec.2/lvec.2
    
    var.O1 <- k / (2 * L^2 * (k - 1)) * sum((nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2)
    var.O2 <- (2 * k) / (L^2 * (k - 1)) * sum(((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1 - ervec.2)^2)
    
    # 7c. The Exact Striplet Variance (Grid-Shift Loop)
    muvec <- dsm_surface$N_hat
    midvec <- dsm_surface$x
    musum <- sum(muvec, na.rm = TRUE)
    
    # Get bounding box to dynamically calculate expected line lengths
    bbox <- sf::st_bbox(region)
    y_length <- as.numeric(bbox["ymax"] - bbox["ymin"])
    x_min <- as.numeric(bbox["xmin"])
    x_max <- as.numeric(bbox["xmax"])
    
    B <- 50 # Number of sub-intervals to approximate the continuous integration
    bvec <- seq(0, spacing, length.out = B)
    Lbvec <- rep(0, B)
    Abvec <- rep(0, B)
    
    for (b_idx in 1:B) {
      b_val <- bvec[b_idx]
      # Deploy theoretical grid for this shift
      lines_grid <- seq(x_min + b_val, x_max, by = spacing)
      Lbvec[b_idx] <- length(lines_grid) * y_length
      
      # Distance from every striplet to the nearest line in THIS specific grid
      min_dist <- sapply(midvec, function(m) min(abs(m - lines_grid)))
      
      # Probability of detection given the fitted half-normal model
      g_b <- exp(- (min_dist^2) / (2 * sigma_hat^2))
      g_b[min_dist > truncation] <- 0 # Absolute truncation
      
      # A(b) = Expected total detections for this grid alignment
      Abvec[b_idx] <- sum(muvec * g_b, na.rm = TRUE)
    }
    
    # Fewster's Equation: Var(n/L) = Var_E(A(b)/L(b)) + E(Var_n|E(A(b)/L(b)))
    var_ER_striplet <- mean((Abvec + (Abvec^2) * (1 - 1/musum)) / (Lbvec^2)) - (mean(Abvec / Lbvec))^2
    
    # Bundle variances into a list
    analytical_variances <- list(
      var_ER_R2 = var.R2,
      var_ER_R3 = var.R3,
      var_ER_S1 = var.S1,
      var_ER_S2 = var.S2,
      var_ER_O1 = var.O1,
      var_ER_O2 = var.O2,
      var_ER_striplet = var_ER_striplet
    )
  }
  
  # Return final integrated list
  list(
    detection_model = m1, dsm = dsm1, N_hat = N_hat, se_m1 = se_ds, 
    sigma_hat = sigma_hat, obsdata = obsdata, segdata = segdata, 
    density_surface = dsm_surface, density = gam.density, 
    population_description = pop.desc, 
    analytical_variances = analytical_variances # New Component!
  )
}