# Load necessary parallel packages
library(tidyverse)
library(future)
library(furrr)
library(parallel)

# Updated to source the new functions file
source(here::here("Functions2.R"))

# Setup Constants & Scenario Grid -----

# region params
lower_x <- 0; upper_x <- 5000
lower_y <- 0; upper_y <- 5000

# Simulation constants
true_N <- 1000
density_grid_spacing <- 500
scale_parameter <- 100
trunc_dist <- 300
design_angle <- 0
MAX_BOOT_REPS <- 1000
design_spacing <- 500

# hotspots_single <- list(
#   list(centre = c(2000, 2000), sigma = 500, amplitude = 5),
#   list(centre = c(8000, 8000), sigma = 500, amplitude = 5),
#   list(centre = c(2000, 8000), sigma = 500, amplitude = 5)
# )

my_hotspots <- list(
  list(centre = c(2500, 2500), sigma = 8000, amplitude = 0.1)
)

# Create a grid of scenarios to test the 3 bootstrap methods on lines and points
scenarios <- expand_grid(
  transect_type = c("line", "point"),
  boot_method = c("Standard", "Discrete")
  # boot_method = c("Standard", "InvCDF", "Discrete")
)

# Convergence Evaluation Function -----

evaluate_convergence <- function(transect_type, boot_method) {
  
  # Note: message() streams better from background workers than cat()
  message(sprintf("Running Scenario: %s Transects with %s Bootstrap...", 
                  toupper(transect_type), toupper(boot_method)))
  
  # 1. Establish Region & Hotspots
  region <- make_region(
    region_sf = NULL, 
    lower_x_bound = lower_x, 
    upper_x_bound = upper_x, 
    lower_y_bound = lower_y, 
    upper_y_bound = upper_y, 
    units = "m"
  )
  
  # 2. Generate Truth & Survey using the single density
  sim_truth <- generate_simulated_pop(
    region = region, 
    N = true_N, 
    x_space = density_grid_spacing, 
    hotspots = my_hotspots, 
    scale_param = scale_parameter, 
    truncation = trunc_dist
  )
  
  sim_survey <- generate_survey_data(
    region = region, 
    realized_population = sim_truth$population,
    angle = design_angle,            
    transect_type = transect_type,
    spacing = design_spacing,
    truncation = trunc_dist
  )
  
  # 3. Fit Models (Replacing the old get_fit wrapper)
  ds_results <- fit_ds_model(
    dist_data = sim_survey$dist_data,
    transect_type = transect_type,
    truncation = trunc_dist
  )
  
  dsm_results <- fit_dsm_model(
    dist_data = sim_survey$dist_data,
    transects = sim_survey$transects,
    region = region,
    detection_model = ds_results$detection_model,
    N_hat = ds_results$N_hat,
    transect_type = transect_type,
    truncation = trunc_dist,
    x_space = density_grid_spacing, 
    y_space = density_grid_spacing  
  )
  
  # analytical_variances <- calculate_variance_estimators(
  #   obsdata = dsm_results$obsdata,
  #   segdata = dsm_results$segdata,
  #   region = region,
  #   detection_model = ds_results$detection_model,
  #   dsm_model = dsm_results$dsm,
  #   dsm_surface = dsm_results$density_surface,
  #   pred_grid = dsm_results$pred_grid,
  #   N_hat = ds_results$N_hat,
  #   sigma_hat = ds_results$sigma_hat,
  #   transect_type = transect_type,
  #   spacing = spacing,
  #   truncation = trunc_dist
  # )
  
  # 4. Dynamically select and run the correct bootstrap function
  boot_func <- switch(boot_method,
                      "Standard" = get_bootstrap,
                      # "InvCDF"   = get_bootstrap_invcdf,
                      "Discrete" = get_bootstrap_disc_density
  )
  
  # Pass the targeted outputs from the new split functions
  boot_res <- boot_func(
    region = region,
    population_description = dsm_results$population_description,
    sigma_hat = ds_results$sigma_hat,
    transect_type = transect_type,
    reps = MAX_BOOT_REPS,
    angle = design_angle,
    spacing = design_spacing,
    truncation = trunc_dist
  )
  
  # 5. Calculate Cumulative Standard Deviation
  eval_reps <- 10:MAX_BOOT_REPS
  
  cumulative_sd <- purrr::map_dbl(eval_reps, function(i) {
    sd(boot_res$N_hat[1:i], na.rm = TRUE)
  })
  
  # Return data formatted for ggplot
  tibble(
    transect_type = transect_type,
    boot_method = boot_method,
    iterations = eval_reps,
    cumulative_se = cumulative_sd
  )
}

# Execution and Visualization -----

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

# Plot the convergence trajectories
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