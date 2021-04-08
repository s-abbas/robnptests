# # Simulations for asymptotic M-estimators ----
#
# # Load package ----
# library(robTests)
#
# # Normal distribution ----
# set.seed(4815)
# sim.data <- replicate(1000, stats::rnorm(200))
#
# # p-values for different sample sizes
# sample.grid <- expand.grid(m = c(10, 30, 50, 100), n = c(10, 30, 50, 100))
# sample.grid <- sample.grid[sample.grid[, 1] <= sample.grid[, 2], ]
# estimators <- c("huber", "hampel", "bisquare")
#
# # Dataframe to save results.normal
# results.normal <- data.frame(m = rep(sample.grid[, "m"], 3), n = rep(sample.grid[, "n"], 3), estimator = rep(estimators, each = 10), p.value = numeric(30), std.error = numeric(30))
#
# # Simulation study
# for (i in 1:nrow(results.normal)) {
#   # Sample sizes
#   m <- results.normal[i, "m"]
#   n <- results.normal[i, "n"]
#
#   # p-values
#   p.values <- apply(sim.data, 2, function(z) m_test(x = z[1:m], y = z[(m + 1):(m + n)], psi = results.normal[i, "estimator"], method = "asymptotic")$p.value)
#
#   # Save p-value and standard error
#   results.normal[i, "p.value"] <- mean(p.values < 0.05)
#   results.normal[i, "std.error"] <- sd(p.values < 0.05)/sqrt(nrow(sim.data))
# }
#
# results.normal <- results.normal[order(results.normal[, "m"]), ]
#
# # Chi2 distribution ----
# set.seed(4815)
# sim.data <- replicate(1000, stats::rchisq(200, df = 3))
#
# # p-values for different sample sizes
# sample.grid <- expand.grid(m = c(10, 30, 50, 100), n = c(10, 30, 50, 100))
# sample.grid <- sample.grid[sample.grid[, 1] <= sample.grid[, 2], ]
# estimators <- c("huber", "hampel", "bisquare")
#
# # Dataframe to save results.chi
# results.chi <- data.frame(m = rep(sample.grid[, "m"], 3), n = rep(sample.grid[, "n"], 3), estimator = rep(estimators, each = 10), p.value = numeric(30), std.error = numeric(30))
#
# # Simulation study
# for (i in 1:nrow(results.chi)) {
#   # Sample sizes
#   m <- results.chi[i, "m"]
#   n <- results.chi[i, "n"]
#
#   # p-values
#   p.values <- apply(sim.data, 2, function(z) m_test(x = z[1:m], y = z[(m + 1):(m + n)], psi = results.chi[i, "estimator"], method = "asymptotic")$p.value)
#
#   # Save p-value and standard error
#   results.chi[i, "p.value"] <- mean(p.values < 0.05)
#   results.chi[i, "std.error"] <- sd(p.values < 0.05)/sqrt(nrow(sim.data))
# }
#
# results.chi <- results.chi[order(results.chi[, "m"]), ]
#
# # Save results ----
# usethis::use_data(results.normal, results.chi, internal = TRUE, overwrite = TRUE)
