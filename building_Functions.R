## Load packages

library(tidyverse)
library(dsims)
library(Distance)
library(dssd)
library(dsm)
library(sf)
library(mgcv)

### Example Constants

true_N <- 5000

points_or_lines = "line"

trunc_dist <- 750
scale_parameter <- 200
design_spacing <- 750

density_grid_spacing <- 500

design_angle <- 0
bootstrap_reps <- 5


## Make region (common to all analysis steps)

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

## Test

region <- make_region(0, 10000, 0, 10000, "m")

plot(region)

## Make theoretical density of animals

make_density <- function(region,
                         hotspots,
                         x_space = density_grid_spacing) {
  
  # Create basic uniform density
  density_true <- dsims::make.density(region = region,
                                      x.space = x_space,
                                      constant = 1)
  
  # 3. Iteratively add hotspots to create spatial variability
  for (hotspot in hotspots) {
    density_true <- dsims::add.hotspot(object = density_true,
                                       centre = hotspot$centre,
                                       sigma = hotspot$sigma,
                                       amplitude = hotspot$amplitude)
  }
  
  return(density_true)
}

## Test

my_hotspots <- list(
  list(centre = c(2000, 8000), sigma = 800, amplitude = 1.2),
  list(centre = c(5000, 2000), sigma = 6000, amplitude = 0.5),
  list(centre = c(8000, 9000), sigma = 1000, amplitude = 2)
)

density_true <- make_density(region, my_hotspots)

plot(density_true)



## Make survey design

make_design <- function(region,
                        transect_type = points_or_lines,
                        spacing = design_spacing,
                        angle = design_angle,
                        truncation = trunc_dist) {
  
  make.design(
    region        = region,
    transect.type = transect_type,  # Evaluates the literal variable passed ("point" or "line")
    design        = "systematic",
    spacing       = spacing,
    edge.protocol = "minus",
    design.angle  = angle,
    truncation    = truncation
  )
  
}



make_ds_analysis <- function(truncation = trunc_dist) {
  
  make.ds.analysis(
    dfmodel = ~ 1,
    key = "hn",
    truncation = truncation,
    criteria = "AIC"
  )
  
}



make_detectability <- function(scale_param = scale_parameter,
                               truncation = trunc_dist) {
  
  make.detectability(
    key.function = "hn",
    scale.param = scale_param,
    truncation = truncation
  )
}

## Test 

design <- make_design(region, "line", 1000, 0, 750)

analysis <- make_ds_analysis(750)

detectability <- make_detectability(200, 750)


### Functions

generate_simulated_truth <- function(region,
                                     N = true_N,
                                     x_space = density_grid_spacing,
                                     hotspots = my_hotspots,
                                     scale_param = scale_parameter,
                                     truncation = trunc_dist) {
  
  density_true <- make.density(
    region = region,
    x.space = x_space,
    constant = 1
  )
  
  for (hotspot in hotspots) {
    density_true <- add.hotspot(
      object = density_true,
      centre = hotspot$centre,
      sigma = hotspot$sigma,
      amplitude = hotspot$amplitude
    )
  }
  
  # fit.gam <- gam(
  #   density ~ s(x, y),
  #   data = density_true@density.surface[[1]],
  #   family = gaussian(link = "log")
  # )
  # 
  # gam.density_true <- make.density(
  #   region = region,
  #   x.space = x_space,
  #   fitted.model = fit.gam
  # )
  
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
  
  detect_true <- make_detectability(
    scale_param = scale_param,
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
    # gam = fit.gam
  )
}



generate_survey_data <- function(
    region,
    realized_population, # output from (generate_simulated_truth_object)$population
    angle = design_angle,
    transect_type = points_or_lines,
    spacing = design_spacing,
    truncation = trunc_dist) 
{
  
  design <- make_design(
    region = region,
    transect_type = transect_type,
    angle = angle,
    spacing = spacing,
    truncation = truncation
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
  
  dsm_surface <- pred_grid |>
    mutate(density = pmax(N_hat / area, .Machine$double.eps)) |>
    select(strata, density, x, y, geometry)
  
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






get_bootstrap <- function(region,
                          population_description,
                          sigma_hat,
                          transect_type = points_or_lines,
                          reps = bootstrap_reps,
                          angle = design_angle,
                          spacing = design_spacing,
                          truncation = trunc_dist) {
  
  detect <- make_detectability(
    scale_param = sigma_hat,
    truncation = truncation
  )
  
  analyses <- make_ds_analysis(truncation = truncation)
  
  # FIXED: Renamed from point.design to design so make.simulation finds it
  design <- make_design(
    region = region,
    angle = angle,
    transect_type = transect_type,
    spacing = spacing,
    truncation = truncation
  )
  
  sim_obj <- make.simulation(
    reps = reps,
    design = design, # Now correctly points to the object created above
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




### Test Cases

# -------------------------------------------------------------------------
# 1. Setup Common Environment
# -------------------------------------------------------------------------
cat("Setting up region and theoretical truth...\n")

region <- make_region(0, 10000, 0, 10000, "m")
truth <- generate_simulated_truth(region)

# Visualize the baseline reality (Removed 'main' to prevent S4 method errors)
plot(truth$density, region, scale = 0.001)
plot(truth$detectability, truth$population_description)

p_truth <- ggplot() +
  geom_sf(data = region@region) +
  geom_point(data = truth$population@population,
             aes(x = x, y = y), size = 0.8) +
  coord_sf() +
  theme_bw() +
  ggtitle("Theoretical Population Realization")
print(p_truth)


# -------------------------------------------------------------------------
# TEST CASE 1: LINE TRANSECTS
# -------------------------------------------------------------------------
cat("\n--- Running Test Case 1: Line Transects ---\n")

survey_data_lines <- generate_survey_data(
  region = region,
  realized_population = truth$population,
  transect_type = "line"
)

# p_lines <- ggplot() + 
#   geom_sf(data = region@region, fill = NA) + 
#   geom_sf(data = survey_data_lines$transects@samplers, color = "blue") + 
#   geom_point(data = truth$population@population |> 
#                mutate(detected = object %in% 
#                         survey_data_lines$dist_data$object[!is.na(survey_data_lines$dist_data$distance)]), 
#              aes(x = x, y = y, color = detected), size = 0.5, alpha = 0.6) + 
#   scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) + 
#   coord_sf() + 
#   theme_bw() +
#   ggtitle("Line Transect Survey")
# print(p_lines)

fit_lines <- get_fit(
  dist_data = survey_data_lines$dist_data,
  transects = survey_data_lines$transects,
  region = region,
  transect_type = "line"
)

if (!is.null(fit_lines$dsm)) {
  cat("Line Transect DSM Fit Successfully.\n")
  plot(fit_lines$density, region, scale = 0.001)
} else {
  cat("Line DSM failed to fit (likely zero observations).\n")
}


# -------------------------------------------------------------------------
# TEST CASE 2: POINT TRANSECTS
# -------------------------------------------------------------------------
cat("\n--- Running Test Case 2: Point Transects ---\n")

survey_data_points <- generate_survey_data(
  region = region,
  realized_population = truth$population,
  transect_type = "point"
)

# p_points <- ggplot() + 
#   geom_sf(data = region@region, fill = NA) + 
#   geom_sf(data = survey_data_points$transects@samplers, shape = 3, color = "blue") + 
#   geom_point(data = truth$population@population |> 
#                mutate(detected = object %in% 
#                         survey_data_points$dist_data$object[!is.na(survey_data_points$dist_data$distance)]), 
#              aes(x = x, y = y, color = detected), size = 0.5, alpha = 0.6) + 
#   scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) + 
#   coord_sf() + 
#   theme_bw() +
#   ggtitle("Point Transect Survey")
# print(p_points)

fit_points <- get_fit(
  dist_data = survey_data_points$dist_data,
  transects = survey_data_points$transects,
  region = region,
  transect_type = "point"
)

if (!is.null(fit_points$dsm)) {
  cat("Point Transect DSM Fit Successfully.\n")
  plot(fit_points$density, region, scale = 0.001)
} else {
  cat("Point DSM failed to fit (likely zero observations).\n")
}


# -------------------------------------------------------------------------
# TEST CASE 3: BOOTSTRAP PIPELINE CHECK
# -------------------------------------------------------------------------
cat("\n--- Running Test Case 3: Small Bootstrap Pipeline Check ---\n")

if (!is.null(fit_lines$dsm)) {
  bootstrap_check <- get_bootstrap(
    region = region,
    population_description = fit_lines$population_description,
    sigma_hat = fit_lines$sigma_hat,
    transect_type = "line",
    reps = 2
  )
  cat("Bootstrap completed successfully. Sample N_hats:\n")
  print(bootstrap_check$N_hat)
}


