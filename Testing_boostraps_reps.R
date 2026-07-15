# Load necessary parallel packages
library(tidyverse)
library(future)
library(furrr)
library(parallel)

source(here::here("Functions4.R"))

# Setup Constants & Scenario Grid -----

# Region parameters
lower_x <- 0; upper_x <- 5000
lower_y <- 0; upper_y <- 5000

# Simulation constants
true_N <- 1000
density_grid_spacing <- 100
scale_parameter <- 25
trunc_dist <- 80
design_angle <- 0
MAX_BOOT_REPS <- 200
design_spacing <- 500

# Uniform density
my_hotspots <- list(
  list(centre = c(2500, 2500), sigma = 8000, amplitude = 0.1)
)

# Creates the 4 scenarios to test: combinations of transect type and MC method
scenarios <- expand_grid(
  transect_type = c("line", "point"),
  boot_method = c("Standard", "Discrete")
)

# Test Survey ----

test_transect_type <- "line"

test_region <- make_region(
  region_sf = NULL,
  lower_x_bound = lower_x,
  upper_x_bound = upper_x,
  lower_y_bound = lower_y,
  upper_y_bound = upper_y,
  units = "m"
)

test_density_true <- make_hotspot_density(
  region = test_region,
  hotspots = my_hotspots,
  x_space = density_grid_spacing
)

test_sim_truth <- generate_simulated_pop(
  region = test_region,
  density_true = test_density_true,
  N = true_N,
  scale_param = scale_parameter,
  truncation = trunc_dist
)

test_population_description <- test_sim_truth$population_description
test_detectability <- test_sim_truth$detectability

test_survey <- survey_data(
  region = test_region,
  realized_population = test_sim_truth$population,
  # angle = design_angle,
  transect_type = test_transect_type,
  spacing = design_spacing,
  truncation = trunc_dist
)

test_ds <- fit_ds(
  region = test_region,
  dist_data = test_survey$dist_data,
  obsdata = test_survey$obsdata,
  segdata = test_survey$segdata,
  transect_type = test_transect_type,
  truncation = trunc_dist
)

test_design <- test_survey$design

test_analyses <- make.ds.analysis(dfmodel = ~ 1,
                                  key = "hn",
                                  truncation = trunc_dist,
                                  criteria = "AIC")

test_sim <- make.simulation(reps = 3,
                            design = test_design,
                            population.description = test_population_description,
                            detectability = test_detectability,
                            ds.analysis = test_analyses)

# Generate a single instance of a survey: a population, set of transects
# and the resulting distance data
test_survey <- run.survey(test_sim)

plot(test_survey, test_region)

density_sf <- test_density_true@density.surface[[1]]

ggplot(density_sf) +
  geom_sf(aes(fill = density), color = NA) +
  scale_fill_viridis_c(option = "viridis", name = "Density") +
  theme_bw() +
  labs(title = "Simulated Density Surface")

test_ds$sigma_hat
test_ds$N_hat

# Convergence Evaluation Function -----

evaluate_convergence <- function(transect_type, boot_method) {
  
  message(sprintf("Running Scenario: %s Transects with %s Bootstrap...", 
                  toupper(transect_type), toupper(boot_method)))
  
  # Build region
  region <- make_region(
    region_sf = NULL, 
    lower_x_bound = lower_x, 
    upper_x_bound = upper_x, 
    lower_y_bound = lower_y, 
    upper_y_bound = upper_y, 
    units = "m"
  )
  
  # Generate density, population, survey
  density_true <- make_hotspot_density(
    region = region,
    hotspots = my_hotspots,
    x_space = density_grid_spacing
  )
  
  sim_truth <- generate_simulated_pop(
    region = region,
    density_true = density_true,
    N = true_N,
    scale_param = scale_parameter,
    truncation = trunc_dist
  )
  
  sim_survey <- survey_data(
    region = region,
    realized_population = sim_truth$population,
    # angle = design_angle,
    transect_type = transect_type,
    spacing = design_spacing,
    truncation = trunc_dist
  )
  
  # Fit models 
  ds_results <- fit_ds(
    region = region,
    dist_data = sim_survey$dist_data,
    obsdata = sim_survey$obsdata,
    segdata = sim_survey$segdata,
    transect_type = transect_type,
    truncation = trunc_dist
  )
  
  dsm_results <- fit_dsm_and_surface(
    obsdata = sim_survey$obsdata,
    segdata = sim_survey$segdata,
    ds_model = ds_results$ds_model,
    region = region,
    N_hat = ds_results$N_hat,
    transect_type = transect_type,
    x_space = density_grid_spacing,
    y_space = density_grid_spacing
  )
  
  # Run each bootstrap function
  boot_func <- switch(boot_method,
                      "Standard" = get_bootstrap,
                      "Discrete" = get_bootstrap_disc_density
  )
  
  # UPDATE: Added spacing and truncation to align with Functions4.R signature
  boot_res <- boot_func(
    region = region,
    population_description = dsm_results$population_description,
    sigma_hat = ds_results$sigma_hat,
    transect_type = transect_type,
    reps = MAX_BOOT_REPS,
    spacing = design_spacing, # Updated for Functions4.R
    truncation = trunc_dist   # Updated for Functions4.R
  )
  
  # 5. Calculate cumulative standard deviation
  eval_reps <- 50:MAX_BOOT_REPS
  
  cumulative_sd <- purrr::map_dbl(eval_reps, function(i) {
    sd(boot_res$N_hat[1:i], na.rm = TRUE)
  })
  
  tibble(
    transect_type = transect_type,
    boot_method = boot_method,
    iterations = eval_reps,
    cumulative_se = cumulative_sd
  )
}

# Running simulations -----

# Setup Parallel Backend
total_cores <- parallel::detectCores()
target_cores <- max(1, total_cores - 2) # All but 2, minimum of 1

message(sprintf("Starting parallel processing on %d out of %d available cores...", target_cores, total_cores))
future::plan(future::multisession, workers = target_cores)

# Use furrr::future_pmap_dfr to loop over the scenarios dataframe in parallel
convergence_results <- furrr::future_pmap_dfr(
  scenarios,
  evaluate_convergence,
  .options = furrr::furrr_options(seed = TRUE) # CRITICAL: Ensures valid RNG across cores
)

# Plot results ----

# Plot the convergence
convergence_plot <- ggplot(convergence_results, aes(x = iterations, y = cumulative_se, color = boot_method)) +
  geom_line(linewidth = 1) + 
  facet_wrap(~ transect_type, scales = "free_y", labeller = label_both) +
  labs(
    title = "Bootstrap Variance Convergence by Method",
    subtitle = "Comparing Standard, and Discrete Density Bootstraps on a fixed population",
    x = "Number of Bootstrap Iterations",
    y = "Estimated Standard Error (Cumulative)",
    color = "Bootstrap Method"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )

print(convergence_plot)

# Return backend to sequential processing when done
future::plan(future::sequential)