# Define the equations
equation <- function(params, Data, R) {
  lambda1 <- params[1]
  lambda2 <- params[2]
  alpha <- params[3]
  pi1 <- params[4]
  pi2 <- 1 - pi1
  
  # Calculate the negative log-likelihood based on the equations
  eq1 <- 8 / lambda1 - sum(8 * (lambda1 - 1) * Data) - sum((15 * (2 * (R + 1) - 1)) * pi1 * exp(-Data * alpha * lambda1) * Data) / (alpha * (pi1 * exp(-Data * alpha * lambda1) + pi2 * exp(-Data * alpha * lambda2)))
  eq2 <- 7 / lambda2 - sum(7 * (lambda2 - 1) * Data) - sum((15 * (2 * (R + 1) - 1)) * pi2 * exp(-Data * alpha * lambda2) * Data) / (alpha * (pi1 * exp(-Data * alpha * lambda1) + pi2 * exp(-Data * alpha * lambda2)))
  eq3 <- 8 / lambda1 - 7 / (1 - pi1) + sum((15 * (2 * (R + 1) - 1) * exp(-Data * alpha * (lambda1 + lambda2))) / (alpha * (pi1 * exp(-Data * alpha * lambda1) + pi2 * exp(-Data * alpha * lambda2))))
  eq4 <- 7 / lambda2 - 7 / (1 - pi1) + sum((15 * (2 * (R + 1) - 1) * exp(-Data * alpha * (lambda1 + lambda2))) / (alpha * (pi1 * exp(-Data * alpha * lambda1) + pi2 * exp(-Data * alpha * lambda2))))
  
  # Return the negative sum of the equations
  return(-sum(c(eq1, eq2, eq3, eq4)))
}

# Sample data
Data <- c(79, 178, 203, 272, 276, 285, 350, 356, 392, 471, 503, 584, 622, 663, 966)
R <- c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0)

# Initial guess for parameters
init_params <- c(lambda1 = 0.4, lambda2 = 0.4, alpha = 1, pi1 = 0.4)

# Optimize the negative log-likelihood function
result <- optim(init_params, equation, Data = Data, R = R)

# Estimated parameters
lambda1_est <- result$par[1]
lambda2_est <- result$par[2]
alpha_est <- result$par[3]
pi1_est <- result$par[4]
pi2_est <- 1 - pi1_est

# Calculate the log-likelihood value
log_likelihood <- 8 * log(pi1_est) + 7 * log(1 - pi1_est) + 8 * log(lambda1_est) + 7 * log(lambda2_est) +
  sum(-8 * (lambda1_est - 1) * Data - 8 * Data - 7 * (lambda2_est - 1) * Data - 7 * Data +
        15 * (2 * (R + 1) - 1) * (1 / alpha_est) * log(pi1_est * exp(-Data * alpha_est * lambda1_est) +
                                                         pi2_est * exp(-Data * alpha_est * lambda2_est)))

# Calculate AIC value
AIC_value <- 8 - 2 * log_likelihood

# Print the estimated parameters and AIC value
cat("Estimated lambda1:", lambda1_est, "\n")
cat("Estimated lambda2:", lambda2_est, "\n")
cat("Estimated alpha:", alpha_est, "\n")
cat("Estimated pi1:", pi1_est, "\n")
cat("Estimated pi2:", pi2_est, "\n")
cat("Log-Likelihood Value:", log_likelihood, "\n")
cat("AIC Value:", AIC_value, "\n")
