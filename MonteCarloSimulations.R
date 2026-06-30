# Monte Carlo Simulation for the Maximum Estimator of a Uniform Distribution

# Exact Theoretical Expressions -----------------------------------------------

# Theoretical formulas

# Bias function of the maximum estimator
bias_exact <- function(theta, n) {
  -theta / (n + 1)
}

# Variance function of the maximum estimator
var_exact <- function(theta, n) {
  n * theta^2 / ((n + 1)^2 * (n + 2))
}

# MSE function of the maximum estimator
mse_exact <- function(theta, n) {
  2 * theta^2 / ((n + 1) * (n + 2))
}


# Monte Carlo Simulation -------------------------------------------------------

# Set random seed for reproducibility
set.seed(123)

# Baseline setting: theta = 1, sample size n = 10, and B = 10000 replications.
theta <- 1
n <- 10
B <- 10000

# Generate B Monte Carlo samples.
# For each replication, draw n observations from U(0, theta)
# and store the sample maximum as the estimator.
t_hat <- replicate(B, {
  y <- runif(n, min = 0, max = theta)
  max(y)
})

# Monte Carlo estimates of bias, variance, and MSE
mc_bias <- mean(t_hat - theta)
mc_var <- var(t_hat)
mc_mse <- mean((t_hat - theta)^2)

mc_bias
mc_var
mc_mse

# Theoretical values
bias_exact(theta, n)
var_exact(theta, n)
mse_exact(theta, n)

# The Monte Carlo estimates closely match the corresponding theoretical values.
# This indicates that the simulation reproduces the analytical properties of
# the maximum estimator under the U(0, theta) model.


# Effect of Sample Size on Bias ------------------------------------------------

# Apply the Monte Carlo simulation for different sample sizes.
n_values <- c(5, 10, 20, 50, 100, 250)

bias_data <- data.frame()

for (n in n_values) {
  t_hat <- replicate(B, {
    y <- runif(n, min = 0, max = theta)
    max(y)
  })
  
  # Bias values for the B simulated maximum estimators:
  # max(Y_1) - theta, max(Y_2) - theta, ..., max(Y_B) - theta
  bias_values <- t_hat - theta
  
  bias_data <- rbind(
    bias_data,
    data.frame(n = factor(n), bias = bias_values)
  )
}

# Exact bias values for each sample size
exact_bias_df <- data.frame(
  n = factor(n_values),
  exact_bias = bias_exact(theta, n_values)
)

# Boxplot of simulated bias distributions
boxplot(
  bias ~ as.factor(n),
  data = bias_data,
  xlab = "Sample size n",
  ylab = "Bias = max(Y) - theta",
  main = "Bias distribution for different n",
  names = n_values
)

# Add exact bias values to the boxplot
exact_bias_vals <- bias_exact(theta, n_values)

points(
  x = 1:length(n_values),
  y = exact_bias_vals,
  col = "red",
  pch = 19,
  cex = 1.5
)

# The distribution of the bias becomes increasingly concentrated around zero
# as the sample size increases. The simulated bias closely follows the
# theoretical expectation, showing convergence toward the exact bias.


# Bias Approximation Error -----------------------------------------------------

# Construct 95% confidence intervals using Monte Carlo standard errors.
# alpha = 0.05 gives the critical value z_(1-alpha/2), approximately 1.96.
alpha <- 0.05
z_alpha <- qnorm(1 - alpha / 2)

bias_results <- data.frame()

for (n in n_values) {
  t_hat <- replicate(B, {
    y <- runif(n, 0, theta)
    max(y)
  })
  
  # Monte Carlo bias
  mc_bias <- mean(t_hat) - theta
  
  # Exact theoretical bias
  ex_bias <- bias_exact(theta, n)
  
  # Difference between simulated and theoretical bias
  diff_bias <- mc_bias - ex_bias
  
  # Monte Carlo standard error
  se_bias <- sd(t_hat - theta) / sqrt(B)
  
  bias_results <- rbind(
    bias_results,
    data.frame(
      n = n,
      mc_bias = mc_bias,
      ex_bias = ex_bias,
      diff_bias = diff_bias,
      se_bias = se_bias
    )
  )
}

# Difference between Monte Carlo bias and exact bias
plot(
  bias_results$n,
  bias_results$diff_bias,
  pch = 19,
  type = "b",
  xlab = "n",
  ylab = "MC bias - exact bias",
  main = "Difference between MC bias and exact bias"
)

# 95% confidence intervals: estimate ± 1.96 * SE
arrows(
  bias_results$n,
  bias_results$diff_bias - z_alpha * bias_results$se_bias,
  bias_results$n,
  bias_results$diff_bias + z_alpha * bias_results$se_bias,
  angle = 90,
  code = 3,
  length = 0.05
)

# The observed differences remain close to zero and lie within the simulation
# error bounds. This indicates that the Monte Carlo procedure accurately
# approximates the theoretical bias of the estimator.


# Variance Approximation Error -------------------------------------------------

# Compare the Monte Carlo estimate of the variance with the theoretical variance.
var_results <- data.frame()

for (n in n_values) {
  t_hat <- replicate(B, {
    y <- runif(n, 0, theta)
    max(y)
  })
  
  # Monte Carlo variance
  mc_var <- var(t_hat)
  
  # Exact theoretical variance
  ex_var <- var_exact(theta, n)
  
  # Difference between simulated and theoretical variance
  diff_var <- mc_var - ex_var
  
  # Standard error for the variance approximation
  dev2 <- (t_hat - mean(t_hat))^2
  se_var <- sd(dev2) / sqrt(B)
  
  var_results <- rbind(
    var_results,
    data.frame(
      n = n,
      diff_var = diff_var,
      se_var = se_var
    )
  )
}

# Difference between Monte Carlo variance and exact variance
plot(
  var_results$n,
  var_results$diff_var,
  type = "b",
  pch = 19,
  xlab = "n",
  ylab = "MC var - exact var",
  main = "Difference between MC variance and exact variance"
)

# 95% confidence intervals
arrows(
  var_results$n,
  var_results$diff_var - z_alpha * var_results$se_var,
  var_results$n,
  var_results$diff_var + z_alpha * var_results$se_var,
  angle = 90,
  code = 3,
  length = 0.05
)

# As n increases, the theoretical variance of the maximum estimator decreases.
# The Monte Carlo estimates remain close to the theoretical variance, and the
# observed differences fall within the simulation error bounds.


# Mean Squared Error Approximation --------------------------------------------

# Compare the Monte Carlo estimate of MSE with the theoretical MSE.
mse_results <- data.frame()

for (n in n_values) {
  t_hat <- replicate(B, {
    y <- runif(n, 0, theta)
    max(y)
  })
  
  # Monte Carlo MSE
  mc_mse <- mean((t_hat - theta)^2)
  
  # Exact theoretical MSE
  ex_mse <- mse_exact(theta, n)
  
  # Difference between simulated and theoretical MSE
  diff_mse <- mc_mse - ex_mse
  
  # Standard error for the MSE approximation
  sq_err <- (t_hat - theta)^2
  se_mse <- sd(sq_err) / sqrt(B)
  
  mse_results <- rbind(
    mse_results,
    data.frame(
      n = n,
      diff_mse = diff_mse,
      se_mse = se_mse
    )
  )
}

# Difference between Monte Carlo MSE and exact MSE
plot(
  mse_results$n,
  mse_results$diff_mse,
  type = "b",
  pch = 19,
  xlab = "n",
  ylab = "MC MSE - exact MSE",
  main = "Difference between MC MSE and exact MSE"
)

# 95% confidence intervals
arrows(
  mse_results$n,
  mse_results$diff_mse - z_alpha * mse_results$se_mse,
  mse_results$n,
  mse_results$diff_mse + z_alpha * mse_results$se_mse,
  angle = 90,
  code = 3,
  length = 0.05
)

# The observed MSE differences remain small and lie within the simulation
# error bounds across the sample sizes considered.


# Sensitivity to the Parameter theta -------------------------------------------

# Repeat the simulation for different theta values while keeping n fixed.
theta_values <- c(0.25, 0.5, 1, 2, 5)
n <- 10

theta_results <- data.frame()

for (theta in theta_values) {
  t_hat <- replicate(B, {
    y <- runif(n, 0, theta)
    max(y)
  })
  
  mc_bias <- mean(t_hat - theta)
  ex_bias <- bias_exact(theta, n)
  
  mc_var <- var(t_hat)
  ex_var <- var_exact(theta, n)
  
  mc_mse <- mean((t_hat - theta)^2)
  ex_mse <- mse_exact(theta, n)
  
  theta_results <- rbind(
    theta_results,
    data.frame(
      theta = theta,
      mc_bias = mc_bias,
      ex_bias = ex_bias,
      mc_var = mc_var,
      ex_var = ex_var,
      mc_mse = mc_mse,
      ex_mse = ex_mse
    )
  )
}

theta_results

# The simulation remains accurate across different theta values. Although bias,
# variance, and MSE change with theta, the differences between theoretical and
# simulated values remain small.