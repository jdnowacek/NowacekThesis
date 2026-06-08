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

make_region <- function(lower_x_bound,
                        upper_x_bound,
                        lower_y_bound,
                        upper_y_bound,
                        units) {
  
  # create survey sf object
  simsquare <- sf::st_bbox(c(xmin = lower_x_bound, xmax = upper_x_bound, 
                             ymin = lower_y_bound, ymax = upper_y_bound)) |> 
    sf::st_as_sfc() |> 
    sf::st_as_sf()
  
  # make region with dssd
  region <- dssd::make.region(region.name = "region",
                              shape = simsquare,
                              units = units)
}

## Make theoretical density of animals

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

## Generates true density surface, population description, detectability,
## population iteration

generate_simulated_truth <- function(region,
                                     N = true_N,
                                     x_space = density_grid_spacing,
                                     hotspots = my_hotspots,
                                     scale_param = scale_parameter,
                                     truncation = trunc_dist) {
  
  density_true <- make_hotspot_density(region = region,
                               hotspots = my_hotspots,
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

## Generates survey data from the population specified in generate simulated truth
## returns design, transects, distribution data from the surveys themselves

generate_survey_data <- function(region,
                                 realized_population, # output from (generate_simulated_truth_object)$population
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

get_fit <- function(dist_data,
                    transects,
                    region,
                    transect_type = points_or_lines,
                    truncation = trunc_dist,
                    x_space = density_grid_spacing,
                    y_space = density_grid_spacing) {
  
  # Fit the distance sampling model
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
  
  # Safely extract standard errors (fall back to NA if numeric(0))
  se_ds <- as.numeric(m1$dht$individuals$N$se)
  se_ds <- if (length(se_ds) > 0) se_ds else NA
  
  sigma_hat <- as.numeric(exp(coef(m1$ddf)$scale["(Intercept)", "estimate"]))
  
  
  obsdata <- dist_data |>
    filter(!is.na(distance)) |>
    mutate(
      object = as.character(object),
      Sample.Label = as.character(Sample.Label),
      size = 1
    ) |>
    select(object, Sample.Label, size, distance)
  
  samplers <- transects@samplers
  
  if (transect_type == "line") {
    # Suppress the centroid warning for constant attributes
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
    select(Sample.Label, Effort, x, y)
  
  # Check if there are valid counts for the DSM
  if (nrow(obsdata) == 0) {
    return(list(
      detection_model = m1, dsm = NULL, N_hat = N_hat,
      se_m1 = se_ds, sigma_hat = sigma_hat,
      obsdata = obsdata, segdata = segdata, density_surface = NULL,
      density = NULL, population_description = NULL
    ))
  }
  
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
  
  prediction_grid <- make.density(
    region = region,
    x.space = x_space,
    y.space = y_space,
    constant = 1
  )
  
  pred_grid <- prediction_grid@density.surface[[1]] |>
    mutate(area = as.numeric(sf::st_area(geometry)))
  
  pred_data <- sf::st_drop_geometry(pred_grid)
  
  pred_grid$N_hat <- predict(
    dsm1,
    newdata = pred_data,
    off.set = pred_grid$area,
    type = "response"
  )
  
  if (transect_type == "line") {
    dsm_surface <- pred_grid |>
      group_by(strata, x) |>
      summarize(
        N_hat = sum(N_hat, na.rm = TRUE),
        area = sum(area, na.rm = TRUE),
        geometry = sf::st_union(geometry), # Combines vertical grid squares into a continuous line strip
        .groups = "drop"
      ) |>
      mutate(density = pmax(N_hat / area, .Machine$double.eps),
             y = 0) |>
      select(strata, density, x, y, geometry)
  } else {
    # Keep standard 2D spatial grid ("boxlets")
    dsm_surface <- pred_grid |>
      mutate(density = pmax(N_hat / area, .Machine$double.eps)) |>
      select(strata, density, x, y, geometry)
  }

  gam.density <- make.density(
    region = region,
    x.space = x_space,
    y.space = y_space,
    density.surface = dsm_surface
  )
  
  pop.desc <- make.population.description(
    region = region,
    density = gam.density,
    N = N_hat,
    fixed.N = TRUE
  )
  
  list(
    detection_model = m1,
    dsm = dsm1,
    N_hat = N_hat,
    se_m1 = se_ds,
    sigma_hat = sigma_hat,
    obsdata = obsdata,
    segdata = segdata,
    density_surface = dsm_surface,
    density = gam.density,
    population_description = pop.desc
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