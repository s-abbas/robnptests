# Simulations for asymptotic M-estimators ----

# Load package ----
library(robTests)

# N(0, 1)-distribution ----
set.seed(168)
sim.data <- replicate(1000, stats::rnorm(200))

# Test size for different sample sizes
sample.grid <- expand.grid(m = c(10, 30, 50, 100), n = c(10, 30, 50, 100))
sample.grid <- sample.grid[sample.grid[, 1] <= sample.grid[, 2], ]
estimators <- c("huber", "hampel", "bisquare")

# Dataframe to save results.normal
results.normal <- data.frame(m = rep(sample.grid[, "m"], 3), n = rep(sample.grid[, "n"], 3), estimator = rep(estimators, each = 10), size = numeric(30), std.error = numeric(30))

# Simulation study
for (i in 1:nrow(results.normal)) {
  # Sample sizes
  m <- results.normal[i, "m"]
  n <- results.normal[i, "n"]

  # Size of the tests
  p.values <- apply(sim.data, 2, function(z) m_test(x = z[1:m], y = z[(m + 1):(m + n)], psi = results.normal[i, "estimator"], method = "asymptotic")$p.value)

  # Save test size and standard error
  results.normal[i, "size"] <- mean(p.values < 0.05)
  results.normal[i, "std.error"] <- sd(p.values < 0.05)/sqrt(nrow(sim.data))
}

# Sort results by sample size
results.normal <- results.normal[order(results.normal[, "m"]), ]

# t(2)-distribution ----
set.seed(168)
sim.data <- replicate(1000, stats::rt(200, df = 2))

# Test size for different sample sizes
sample.grid <- expand.grid(m = c(10, 30, 50, 100), n = c(10, 30, 50, 100))
sample.grid <- sample.grid[sample.grid[, 1] <= sample.grid[, 2], ]
estimators <- c("huber", "hampel", "bisquare")

# Dataframe to save results.t
results.t <- data.frame(m = rep(sample.grid[, "m"], 3), n = rep(sample.grid[, "n"], 3), estimator = rep(estimators, each = 10), size = numeric(30), std.error = numeric(30))

# Simulation study
for (i in 1:nrow(results.t)) {
  # Sample sizes
  m <- results.t[i, "m"]
  n <- results.t[i, "n"]

  # Size of the test
  p.values <- apply(sim.data, 2, function(z) m_test(x = z[1:m], y = z[(m + 1):(m + n)], psi = results.t[i, "estimator"], method = "asymptotic")$p.value)

  # Save p-value and standard error
  results.t[i, "size"] <- mean(p.values < 0.05)
  results.t[i, "std.error"] <- sd(p.values < 0.05)/sqrt(nrow(sim.data))
}

# Sort results by sample size
results.t <- results.t[order(results.t[, "m"]), ]

# Chi2(3)-distribution ----
set.seed(168)
sim.data <- replicate(1000, stats::rchisq(200, df = 3))

# Test size for different sample sizes
sample.grid <- expand.grid(m = c(10, 30, 50, 100), n = c(10, 30, 50, 100))
sample.grid <- sample.grid[sample.grid[, 1] <= sample.grid[, 2], ]
estimators <- c("huber", "hampel", "bisquare")

# Dataframe to save results.chi
results.chi <- data.frame(m = rep(sample.grid[, "m"], 3), n = rep(sample.grid[, "n"], 3), estimator = rep(estimators, each = 10), size = numeric(30), std.error = numeric(30))

# Simulation study
for (i in 1:nrow(results.chi)) {
  # Sample sizes
  m <- results.chi[i, "m"]
  n <- results.chi[i, "n"]

  # Size of the tests
  p.values <- apply(sim.data, 2, function(z) m_test(x = z[1:m], y = z[(m + 1):(m + n)], psi = results.chi[i, "estimator"], method = "asymptotic")$p.value)

  # Save test size and standard error
  results.chi[i, "size"] <- mean(p.values < 0.05)
  results.chi[i, "std.error"] <- sd(p.values < 0.05)/sqrt(nrow(sim.data))
}

# Sort results by sample size
results.chi <- results.chi[order(results.chi[, "m"]), ]

# Save results ----
usethis::use_data(results.normal, results.chi, results.t, internal = TRUE, overwrite = TRUE)
