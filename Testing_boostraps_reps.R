source(here::here("Functions.R"))

# Load necessary parallel packages
library(future)
library(furrr)
library(parallel)

# Setup Constants & Scenario Grid -----

# region params
lower_x <- 0; upper_x <- 10000
lower_y <- 0; upper_y <- 10000

# Simulation constants
true_N <- 1000
density_grid_spacing <- 500
scale_parameter <- 250
trunc_dist <- 500
design_angle <- 0
MAX_BOOT_REPS <- 10000

# Define two different densities
hotspots_mild <- list(
  list(centre = c(5000, 5000), sigma = 3000, amplitude = 1)
)

hotspots_extreme <- list(
  list(centre = c(2000, 2000), sigma = 500, amplitude = 5),
  list(centre = c(8000, 8000), sigma = 500, amplitude = 5),
  list(centre = c(2000, 8000), sigma = 500, amplitude = 5)
)

# Create a grid of scenarios to test
scenarios <- expand_grid(
  transect_type = c("line", "point"),
  hotspot_name = c("Mild", "Extreme")
)

# Convergence Evaluation Function -----

evaluate_convergence <- function(transect_type, hotspot_name) {
  
  # Note: cat() output might not print to console sequentially in parallel processing
  message(sprintf("Running Scenario: %s Transects with %s Hotspots...", 
                  toupper(transect_type), toupper(hotspot_name)))
  
  # 1. Establish Region & Hotspots
  region <- make_region(lower_x, upper_x, lower_y, upper_y, units = "m")
  my_hotspots <- if(hotspot_name == "Mild") hotspots_mild else hotspots_extreme
  
  # 2. Generate Truth & Survey
  sim_truth <- generate_simulated_pop(
    region = region, 
    N = true_N, 
    hotspots = my_hotspots, 
    scale_param = scale_parameter, 
    truncation = trunc_dist
  )
  
  # Adjust spacing so lines and points have reasonable effort
  spacing <- if(transect_type == "line") 1000 else 500 
  
  sim_survey <- generate_survey_data(
    region = region, 
    realized_population = sim_truth$population,
    transect_type = transect_type,
    spacing = spacing,
    truncation = trunc_dist
  )
  
  # 3. Fit Model & Density Surface
  fit_obj <- get_fit(
    dist_data = sim_survey$dist_data,
    transects = sim_survey$transects,
    region = region,
    spacing = spacing,
    transect_type = transect_type,
    truncation = trunc_dist
  )
  
  # 4. Run Maximum Bootstraps
  boot_res <- get_bootstrap(
    region = region,
    population_description = fit_obj$population_description,
    sigma_hat = fit_obj$sigma_hat,
    transect_type = transect_type,
    reps = MAX_BOOT_REPS,
    spacing = spacing,
    truncation = trunc_dist
  )
  
  # 5. Calculate Cumulative Standard Deviation
  # We start at 5 reps because SD is too volatile/undefined for < 5
  eval_reps <- 5:MAX_BOOT_REPS
  
  cumulative_sd <- purrr::map_dbl(eval_reps, function(i) {
    sd(boot_res$N_hat[1:i], na.rm = TRUE)
  })
  
  # Return data formatted for ggplot
  tibble(
    transect_type = transect_type,
    hotspot_name = hotspot_name,
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
  list(scenarios$transect_type, scenarios$hotspot_name),
  evaluate_convergence,
  .options = furrr::furrr_options(seed = TRUE) # CRITICAL: Ensures valid RNG across cores
)

# Plot the convergence trajectories
library(ggplot2)
convergence_plot <- ggplot(convergence_results, aes(x = iterations, y = cumulative_se, color = hotspot_name)) +
  geom_line(linewidth = 1) + # Changed size to linewidth to avoid ggplot2 deprecation warnings
  facet_wrap(~ transect_type, scales = "free_y", labeller = label_both) +
  labs(
    title = "Bootstrap Variance Convergence",
    subtitle = "Standard Error of N_hat stabilizing over iterations",
    x = "Number of Bootstrap Iterations",
    y = "Estimated Standard Error (Cumulative)",
    color = "Distribution Setup"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )

print(convergence_plot)

# Return backend to sequential processing when done
future::plan(future::sequential)