# Load packages. ###
library(tidyverse)
library(dsims)
library(Distance)
library(dssd)
library(dsm)
library(sf)
library(mgcv)

set.seed(1234567890)

### Constants

REGION_NAME <- "simsquare"
REGION_UNITS <- "m"
REGION_BOUNDS <- c(xmin = 0, xmax = 10000, ymin = 0, ymax = 10000)

DENSITY_GRID_SPACING <- 500
HOTSPOTS <- list(
  list(centre = c(2000, 8000), sigma = 800, amplitude = 1.2),
  list(centre = c(5000, 2000), sigma = 6000, amplitude = 0.5),
  list(centre = c(8000, 9000), sigma = 1000, amplitude = 2)
)
TRUE_N <- 1000

TRUNCATION_DISTANCE <- 750
TRUE_DETECTION_SCALE <- 200

POINT_SPACING <- 1500
ANGLE <- 65
BOOTSTRAP_REPS <- 10

### Helper functions

make_point_design <- function(
    region,
    angle = ANGLE,
    spacing = POINT_SPACING,
    truncation = TRUNCATION_DISTANCE) 
    {

  make.design(
    region = region,
    transect.type = "point",
    design = "systematic",
    spacing = spacing,
    edge.protocol = "minus",
    design.angle = angle,
    truncation = truncation
  )
}

make_ds_analysis <- function(truncation = TRUNCATION_DISTANCE) {
  make.ds.analysis(
    dfmodel = ~ 1,
    key = "hn",
    truncation = truncation,
    er.var = "P3",
    criteria = "AIC"
  )
}

make_detectability <- function(
    scale_param,
    truncation = TRUNCATION_DISTANCE) 
    {

  make.detectability(
    key.function = "hn",
    scale.param = scale_param,
    truncation = truncation
  )
}

### Functions

generate_simulated_truth <- function(
    region,
    N = TRUE_N,
    x_space = DENSITY_GRID_SPACING,
    hotspots = HOTSPOTS,
    scale_param = TRUE_DETECTION_SCALE,
    truncation = TRUNCATION_DISTANCE) 
    {

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

  fit.gam <- gam(
    density ~ s(x, y),
    data = density_true@density.surface[[1]],
    family = gaussian(link = "log")
  )

  gam.density_true <- make.density(
    region = region,
    x.space = x_space,
    fitted.model = fit.gam
  )

  density_surface_true <- gam.density_true@density.surface[[1]] |>
    mutate(
      area = as.numeric(sf::st_area(geometry)),
      density = density * N / sum(density * area)
    ) |>
    select(strata, density, x, y, geometry)

  gam.density_true <- make.density(
    region = region,
    x.space = x_space,
    density.surface = density_surface_true
  )

  pop.desc_true <- make.population.description(
    region = region,
    density = gam.density_true,
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
    density = gam.density_true,
    population_description = pop.desc_true,
    population = pop_true,
    detectability = detect_true,
    gam = fit.gam
  )
}

generate_survey_data <- function(
    region,
    realized_population,
    angle = ANGLE,
    spacing = POINT_SPACING,
    truncation = TRUNCATION_DISTANCE) 
    {

  point.design <- make_point_design(
    region = region,
    angle = angle,
    spacing = spacing,
    truncation = truncation
  )

  transects <- generate.transects(point.design)

  survey <- new(
    Class = "Survey.PT",
    population = realized_population,
    transect = transects,
    rad.truncation = truncation
  )

  observed_survey <- run.survey(survey, region = region)

  list(
    design = point.design,
    transects = transects,
    dist_data = observed_survey@dist.data
  )
}

get_fit <- function(
    dist_data,
    transects,
    region,
    truncation = TRUNCATION_DISTANCE,
    x_space = DENSITY_GRID_SPACING,
    y_space = DENSITY_GRID_SPACING) 
    {

  m1 <- ds(
    data = dist_data,
    transect = "point",
    key = "hn",
    adjustment = NULL,
    truncation = truncation
  )

  N_hat <- as.numeric(m1$dht$individuals$N$Estimate)
  N_hat <- as.integer(round(N_hat))

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
  sampler_xy <- sf::st_coordinates(samplers)

  segdata <- samplers |>
    sf::st_drop_geometry() |>
    mutate(
      Sample.Label = as.character(transect),
      Effort = 1,
      x = sampler_xy[, "X"],
      y = sampler_xy[, "Y"]
    ) |>
    select(Sample.Label, Effort, x, y)

  dsm1 <- dsm(
    count ~ s(x, y),
    ddf.obj = m1,
    segment.data = segdata,
    observation.data = obsdata,
    family = quasipoisson(link = "log"),
    method = "REML",
    convert.units = 1
  )

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
    sigma_hat = sigma_hat,
    obsdata = obsdata,
    segdata = segdata,
    density_surface = dsm_surface,
    density = gam.density,
    population_description = pop.desc
  )
}

get_bootstrap <- function(
    region,
    population_description,
    sigma_hat,
    reps = BOOTSTRAP_REPS,
    angle = ANGLE,
    spacing = POINT_SPACING,
    truncation = TRUNCATION_DISTANCE) 
    {

  detect <- make_detectability(
    scale_param = sigma_hat,
    truncation = truncation
  )

  analyses <- make_ds_analysis(truncation = truncation)

  point.design <- make_point_design(
    region = region,
    angle = angle,
    spacing = spacing,
    truncation = truncation
  )

  sim.point <- make.simulation(
    reps = reps,
    design = point.design,
    population.description = population_description,
    detectability = detect,
    ds.analysis = analyses
  )

  sim.point <- run.simulation(sim.point)

  N_hat <- as.numeric(
    sim.point@results[["individuals"]][["N"]][1, "Estimate", seq_len(reps)]
  )

  data.frame(
    replicate = seq_len(reps),
    N_hat = N_hat
  )
}

### Analysis

# set up region, used in both getting truth and in survey
simsquare <- st_bbox(REGION_BOUNDS) |>
  st_as_sfc() |>
  st_as_sf()

region <- make.region(
  region.name = REGION_NAME,
  shape = simsquare,
  units = REGION_UNITS
)

plot(region)

# get truth, this is just for purpose of the simulation
truth <- generate_simulated_truth(region)

plot(truth$density, region, scale = 0.001)
plot(truth$detectability, truth$population_description)
ggplot() + geom_sf(data = region@region) + geom_point(data = truth$population@population, aes(x = x, y = y), size = 0.8) + coord_sf() + theme_bw()

# do the survey 
survey_data <- generate_survey_data(
  region = region,
  realized_population = truth$population
)

plot(region, survey_data$transects) # just the transects
ggplot() + geom_sf(data = region@region, fill = NA) + geom_sf(data = survey_data$transects@samplers, shape = 3, color = "blue") + geom_point(data = truth$population@population |> mutate(detected = individual %in% survey_data$dist_data$individual[!is.na(survey_data$dist_data$distance)]), aes(x = x, y = y, color = detected), size = 0.5, alpha = 0.6) + scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) + coord_sf() + theme_bw()

# get fit, used as truth in bootstrap
fit <- get_fit(
  dist_data = survey_data$dist_data,
  transects = survey_data$transects,
  region = region
)

summary(fit$detection_model)
plot(fit$density, region, scale = 0.001)

# get bootstrap
bootstrap <- get_bootstrap(
  region = region,
  population_description = fit$population_description,
  sigma_hat = fit$sigma_hat
)

sd(bootstrap$N_hat)
