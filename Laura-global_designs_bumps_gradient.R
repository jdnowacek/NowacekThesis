library(dsims)
# Effort allocation formulas

# Universal detectability
detect <- make.detectability(key.function = "hn",
                             scale.param = 10,
                             truncation = 20)

# Analyse using O2 variance estimator
analyses <- make.ds.analysis(dfmodel = ~1,
                             key = c("hn"),
                             truncation = 20,
                             er.var = "O2")

# Analyse using O2 variance estimator
analyses.R2 <- make.ds.analysis(dfmodel = ~1,
                                key = c("hn"),
                                truncation = 20,
                                er.var = "R2")


# Global region
outer1 <- matrix(c(0,0,0,1000,5000,1000,5000,0,0,0),ncol=2, byrow=TRUE)
pol1 <- sf::st_polygon(list(outer1))
sfc <- sf::st_sfc(pol1)
mp1 <- sf::st_sf(geom = sfc)
region.global <- make.region(region.name = "region", 
                             shape = mp1)

# Basic density object
density <- make.density(region = region.global,
                        x.space = 50,
                        constant = 1)

# Global variable distribution example
density.var <- add.hotspot(density, c(600, 200), 200, 1)
density.var <- add.hotspot(density.var, c(1000, 600), 300, 1)
density.var <- add.hotspot(density.var, c(2400, 300), 200, 0.5)
density.var <- add.hotspot(density.var, c(3500, 400), 1000, -0.25)
density.var <- add.hotspot(density.var, c(3800, 100), 800, 0.6)
density.var <- add.hotspot(density.var, c(1000, 600), 350, 1)
density.var <- add.hotspot(density.var, c(2600, 600), 600, 0.5)
density.var <- add.hotspot(density.var, c(4000, 900), 250, 1)
plot(density.var)

# Global gradient distribution example
density.grad <- density
density.grad@density.surface[[1]]$density <- density.grad@density.surface[[1]]$x*8e-04
plot(density.grad)

# Make population descriptions
pop.desc.var <- make.population.description(region = region.global,
                                            density = density.var,
                                            N = c(2500))

# Make designs
design.global <- make.design(region = region.global,
                             line.length = 20000,
                             truncation = 20)

transects <- generate.transects(design.global)
plot(region.global, transects)

# Make simulations
sim.global.var <- make.simulation(reps = 10,
                                  design = design.global,
                                  pop.desc.var,
                                  detectability = detect,
                                  ds.analysis = analyses)

sim.global.var.R2 <- make.simulation(reps = 10,
                                     design = design.global,
                                     pop.desc.var,
                                     detectability = detect,
                                     ds.analysis = analyses.R2)

# Check simulation setup
eg.survey <- run.survey(sim.global.var)
plot(eg.survey)

eg.survey@transect@cov.area/region.global@area

# Run simulations
Sys.time()
sim.global.var <- run.simulation(sim.global.var)
Sys.time()
save(sim.global.var, file = "sim_global_var_Jul26.robj")

# Run simulations
Sys.time()
sim.global.var.R2 <- run.simulation(sim.global.var.R2)
Sys.time()
save(sim.global.var.R2, file = "sim_global_var_R2_Jul26.robj")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots - these are in the read_results.R file

load("sim_global_var_Jul26.robj")
load("sim_global_var_R2_Jul26.robj")

library(dsims)
summary(sim.global.var)
# Using the O2 estimator the 95%CI coverage is 96% where as it is 99% for the R2 estimator, the true standard error (sd of means) is around 200 and we can see that the O2 is a better estimator than the R2 estimator but it is still over-estimating the variability. I wonder how well the striplet estimation method would do at estimating the variability for this scenario.

# Estimates of Abundance (N)
# 
# Truth mean.Estimate percent.bias   RMSE CI.coverage.prob mean.se sd.of.means
# 1  2500       2501.01         0.04 205.94             0.96  227.01      206.14

summary(sim.global.var.R2)
# Estimates of Abundance (N)
# 
# Truth mean.Estimate percent.bias   RMSE CI.coverage.prob mean.se sd.of.means
# 1  2500       2499.57        -0.02 194.93             0.99  249.76      195.12

