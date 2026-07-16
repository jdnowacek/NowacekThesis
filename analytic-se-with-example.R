source("Functions4.R")

analytic_se <- function(region,
                        dsm_model,
                        ds_results,
                        transect_type,
                        spacing,
                        truncation,
                        integration_spacing,
                        phase_grid = 20) {
  
  # Set up the region and integration surface
  N_hat <- ds_results$N_hat
  region_sf <- region@region
  area_region <- region@area
  bbox <- sf::st_bbox(region_sf)
  
  surface <- make.density(
    region = region,
    x.space = integration_spacing,
    y.space = integration_spacing,
    constant = 1
  )@density.surface[[1]] |>
    mutate(area = as.numeric(sf::st_area(geometry)))
  
  cell_abundance <- c(predict(
    dsm_model,
    newdata = sf::st_drop_geometry(surface),
    off.set = surface$area,
    type = "response"
  ))
  cell_abundance <- cell_abundance * N_hat / sum(cell_abundance)
  
  animal_points <- sf::st_as_sf(
    sf::st_drop_geometry(surface),
    coords = c("x", "y"),
    crs = sf::st_crs(region_sf)
  )
  
  # Reconstruct the fitted mrds distance likelihood
  ddf <- ds_results$ds_model$ddf
  theta <- ddf$par
  ddf_object <- ddf$ds$aux$ddfobj
  point_transect <- transect_type == "point"
  left <- ddf$meta.data$left
  first_profile <- seq_len(nrow(ddf_object$xmat)) == 1
  
  detection_model <- function(parameters) {
    mrds:::assign.par(ddf_object, parameters)
  }
  
  detection <- function(distance, parameters) {
    mrds::detfct(
      distance,
      detection_model(parameters),
      index = 1,
      width = truncation,
      left = left
    )
  }
  
  average_detection <- function(parameters) {
    mrds:::integratepdf(
      detection_model(parameters),
      select = first_profile,
      width = truncation,
      int.range = c(left, truncation),
      standardize = FALSE,
      point = point_transect,
      left = left
    )
  }
  
  log_distance_density <- function(distance, parameters) {
    log(detection(distance, parameters)) -
      log(average_detection(parameters))
  }
  
  P <- average_detection(theta)
  gradient_P <- numDeriv::grad(average_detection, theta)
  shifts <- (seq_len(phase_grid) - 0.5) * spacing / phase_grid
  phases <- if (point_transect) {
    expand.grid(x = shifts, y = shifts)
  } else {
    data.frame(x = shifts, y = 0)
  }
  
  phase_results <- data.frame(
    Q = numeric(nrow(phases)),
    sampled_area = numeric(nrow(phases)),
    v_det = numeric(nrow(phases)),
    v_loc = numeric(nrow(phases)),
    M = numeric(nrow(phases)),
    parameter_shift = numeric(nrow(phases))
  )
  phase_scores <- matrix(0, nrow = nrow(phases), ncol = length(theta))
  
  # Evaluate the surface integrals at every phase
  for (j in seq_len(nrow(phases))) {
    if (!point_transect) {
      line_x <- seq(bbox["xmin"] + phases$x[j], bbox["xmax"], by = spacing)
      samplers <- sf::st_sfc(
        lapply(line_x, function(x) {
          sf::st_linestring(matrix(
            c(x, bbox["ymin"], x, bbox["ymax"]),
            ncol = 2,
            byrow = TRUE
          ))
        }),
        crs = sf::st_crs(region_sf)
      )
      samplers <- sf::st_intersection(samplers, region_sf)
      sampled_area <- 2 * (truncation - left) *
        as.numeric(sum(sf::st_length(samplers)))
    } else {
      samplers <- expand.grid(
        x = seq(bbox["xmin"] + phases$x[j], bbox["xmax"], by = spacing),
        y = seq(bbox["ymin"] + phases$y[j], bbox["ymax"], by = spacing)
      ) |>
        sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(region_sf))
      samplers <- samplers[lengths(sf::st_intersects(samplers, region_sf)) > 0, ]
      sampled_area <- pi * (truncation^2 - left^2) * nrow(samplers)
    }
    
    # Each animal can be encountered independently by every nearby sampler.
    distance <- sf::st_distance(animal_points, samplers)
    encounters <- which(distance >= left & distance <= truncation, arr.ind = TRUE)
    cell <- encounters[, "row"]
    distance <- as.numeric(distance[encounters])
    g <- detection(distance, theta)
    weights <- cell_abundance[cell] * g
    
    score <- numDeriv::jacobian(
      function(parameters) log_distance_density(distance, parameters),
      theta
    )
    expected_log_likelihood <- function(parameters) {
      sum(weights * log_distance_density(distance, parameters)) / N_hat
    }
    
    Q <- sum(weights) / N_hat
    S <- colSums(weights * score) / N_hat
    information <- -numDeriv::hessian(expected_log_likelihood, theta)
    b <- Q / P * solve(t(information), gradient_P)
    influence <- drop(1 - score %*% b)
    centre <- Q - sum(b * S)
    estimator_scale <- area_region / (sampled_area * P)
    conditional_mean <- numeric(length(cell_abundance))
    encounter_mean <- rowsum(g * influence, cell)
    conditional_mean[as.integer(rownames(encounter_mean))] <- encounter_mean[, 1]
    
    phase_results$Q[j] <- Q
    phase_results$sampled_area[j] <- sampled_area
    phase_results$v_det[j] <- estimator_scale^2 * sum(
      cell_abundance[cell] * g * (1 - g) * influence^2
    )
    phase_results$v_loc[j] <- estimator_scale^2 * sum(
      cell_abundance * (conditional_mean - centre)^2
    )
    phase_results$M[j] <- estimator_scale * N_hat * centre
    phase_results$parameter_shift[j] <- sqrt(sum(solve(information, S)^2))
    phase_scores[j, ] <- S
  }
  
  variance_components <- c(
    detection = mean(phase_results$v_det),
    location = mean(phase_results$v_loc),
    phase = mean(phase_results$M^2) - mean(phase_results$M)^2
  )
  
  list(
    se = sqrt(sum(variance_components)),
    variance_components = variance_components,
    mean_phase_abundance = mean(phase_results$M),
    mean_phase_score = colMeans(phase_scores),
    max_parameter_shift = max(phase_results$parameter_shift),
    average_detection = P,
    integration_spacing = integration_spacing,
    integration_cells = nrow(surface),
    cells_per_truncation = truncation / integration_spacing,
    phase_results = phase_results
  )
}


# Small line-transect example ##########

set.seed(42)
trunc_dist <- 5
density_grid_spacing <- 10

region <- make_region(
  lower_x_bound = 0,
  upper_x_bound = 100,
  lower_y_bound = 0,
  upper_y_bound = 100,
  units = "m"
)

density <- make.density(region, x.space = density_grid_spacing, constant = 1)
population <- generate_simulated_pop(
  region,
  density,
  N = 1000,
  scale_param = 2,
  truncation = trunc_dist
)

survey <- survey_data(
  region,
  population$population,
  transect_type = "line",
  spacing = 20,
  truncation = trunc_dist
)

ds_results <- fit_ds(
  region,
  survey$dist_data,
  survey$obsdata,
  survey$segdata,
  transect_type = "line",
  truncation = trunc_dist
)

dsm_model <- dsm(
  count ~ 1,
  ddf.obj = ds_results$ds_model,
  segment.data = survey$segdata,
  observation.data = survey$obsdata,
  family = quasipoisson(link = "log"),
  method = "REML",
  convert.units = 1
)

result <- analytic_se(
  region,
  dsm_model,
  ds_results,
  transect_type = "line",
  spacing = 20,
  truncation = trunc_dist,
  integration_spacing = 1, ## the size of the boxlets (eg 1x1 boxes here)
  phase_grid = 100 ## the number of starting positions used in numerical integration 
)

result$se




# now compare to existing code 
dsm_surface <- density@density.surface[[1]] |>
  mutate(area = as.numeric(sf::st_area(geometry)))
dsm_surface$N_hat_pred <- c(predict(
  dsm_model,
  newdata = sf::st_drop_geometry(dsm_surface),
  off.set = dsm_surface$area,
  type = "response"
))

striplet_boxlet(
  region,
  survey$segdata,
  dsm_model,
  dsm_surface,
  ds_results,
  transect_type = "line",
  spacing = 20,
  truncation = trunc_dist
)$se_N_striplet