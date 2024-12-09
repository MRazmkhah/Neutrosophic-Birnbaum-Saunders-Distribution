# ------------------------------------------------------------------------
# Neutrosophic Parameter Estimation and Inferences for Birnbaum-Saunders Model (Real Data 2)
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------
library(GoFKernel)

# ------------------------------------------------------------------------
# Input Data
# ------------------------------------------------------------------------
# Original data points
xl <- c(304.12, 355.34, 310.93, 309.47, 309.12, 292.80, 327.49,
        280.33, 259.99, 238.19, 229.98, 226.77, 223.57, 233.12,
        216.37, 208.16, 206.30, 193.44, 177.31, 157.88, 153.18,
        143.93, 132.78, 127.87, 118.28, 116.75, 116.96, 114.21,
        106.86)

xu <- c(307.82, 355.34, 310.93, 309.47, 312.10, 292.80, 327.49,
        280.33, 259.99, 242.45, 229.98, 226.77, 223.57, 233.12,
        216.37, 208.16, 209.14, 193.44, 177.31, 157.88, 153.18,
        143.93, 132.78, 127.87, 118.28, 116.75, 116.96, 114.21,
        110.62)

# Add small offsets to avoid optimization issues
L <- xl - 1e-9
U <- xu + 1e-9
x0 <- (L + U) / 2
n <- length(xl)  # Sample size

# ------------------------------------------------------------------------
# Parameter Estimation Functions
# ------------------------------------------------------------------------

# Function to optimize beta (ML estimate for beta)
MLb <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- sum(x) / n
  
  DL_b <- function(bb) { # Derivative of likelihood with respect to beta
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }
  
  uniroot(DL_b, lower = r, upper = s)$root
}

# Function to optimize alpha (ML estimate for alpha)
MLa <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- sum(x) / n
  
  beta_hat <- uniroot(function(bb) {
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }, lower = r, upper = s)$root
  
  sqrt(s / beta_hat + beta_hat / r - 2)
}

# ------------------------------------------------------------------------
# Calculate Bounds for Alpha and Beta
# ------------------------------------------------------------------------
ML_b_L <- optim(par = x0, fn = MLb, lower = L, upper = U, method = "L-BFGS-B")$value
ML_b_U <- optim(par = x0, fn = MLb, lower = L, upper = U, method = "L-BFGS-B", control = list(fnscale = -1))$value
ML_a_L <- optim(par = x0, fn = MLa, lower = L, upper = U, method = "L-BFGS-B")$value
ML_a_U <- optim(par = x0, fn = MLa, lower = L, upper = U, method = "L-BFGS-B", control = list(fnscale = -1))$value

# ------------------------------------------------------------------------
# Log-Likelihood Function
# ------------------------------------------------------------------------
loglikelihood <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- sum(x) / n
  
  beta_hat <- uniroot(function(bb) {
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }, lower = r, upper = s)$root
  
  alpha_hat <- sqrt(s / beta_hat + beta_hat / r - 2)
  
  -n * log(2) - n * log(alpha_hat) - n * log(beta_hat) - n * log(2 * pi) / 2 + 
    sum(log(sqrt(beta_hat / x) + (beta_hat / x)^1.5)) -
    (1 / (2 * alpha_hat^2)) * sum(x / beta_hat + beta_hat / x - 2)
}

loglike_min <- optim(par = x0, fn = loglikelihood, lower = L, upper = U, method = "L-BFGS-B")$value
loglike_max <- optim(par = x0, fn = loglikelihood, lower = L, upper = U, method = "L-BFGS-B", control = list(fnscale = -1))$value

# ------------------------------------------------------------------------
# AIC and BIC Calculations
# ------------------------------------------------------------------------
AIC <- function(y) {
  2 * 2 - 2 * y
}

BIC <- function(y) {
  4 * log(n) - 2 * y
}

AIC_min <- AIC(loglike_max)
AIC_max <- AIC(loglike_min)
BIC_min <- BIC(loglike_max)
BIC_max <- BIC(loglike_min)

# ------------------------------------------------------------------------
# Modified KS Test Function
# ------------------------------------------------------------------------
F_t <- function(t, a, b) {
  xx <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
  pnorm(xx, mean = 0, sd = 1)
}

Modified_KS <- function(t,a_hat,b_hat){
  n <- length(t)
  t_sort <- sort(t) #Step 1 in Algorithm 1 
  V_hat <- F_t(t_sort,a_hat,b_hat) #Step 2
  
  Y_hat <- qnorm(V_hat) #Step 3
  Z <- (Y_hat - mean(Y_hat))/sqrt(sum((Y_hat-mean(Y_hat))^2)/(n-1))
  U_hat <- pnorm(Z) #Step 4
  
  F_n <- ecdf(t_sort)
  ECDF <- F_n(t_sort)
  
  D <- NULL
  for (j in 1:n){
    if (j==1){
      D[j] <- abs(U_hat[j]-ECDF[j])
    }
    else{
      
      D[j] <- max(abs(ECDF[j]-U_hat[j]),abs(U_hat[j]-ECDF[j-1]))
      
    }
  }
  
  Max_D <- max(D)
  KS_star <- (sqrt(n)-0.01+0.85/sqrt(n))*Max_D
  
  c(Max_D , KS_star)
}

DDl <- Modified_KS(xl, ML_a_L, ML_b_L)
DDu <- Modified_KS(xu, ML_a_U, ML_b_U)

# ------------------------------------------------------------------------
# P-Value Calculation Using Monte Carlo Simulation
# ------------------------------------------------------------------------
KS_Montecarlo <- function(data, a, b, iterations = 5000) {
  KS_data <- Modified_KS(data, a, b)[2]
  
  F_tt <- function(t) {
    xx <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
    pnorm(xx, mean = 0, sd = 1)
  }
  
  F_t_inv <- inverse(F_tt, lower = 0, upper = 1e5)
  
  RV <- function(n) {
    sapply(runif(n), F_t_inv)
  }
  
  KS_Monte <- replicate(iterations, {
    sample <- RV(length(data))
    Modified_KS(sample, a, b)[2]
  })
  
  mean(KS_Monte > KS_data)
}

set.seed(12345)
p_value_u <- KS_Montecarlo(xl, ML_a_L, ML_b_L)
p_value_l <- KS_Montecarlo(xu, ML_a_U, ML_b_U)

# ------------------------------------------------------------------------
# Final Results
# ------------------------------------------------------------------------
results <- paste0(
  "Estimated alpha = [", round(ML_a_L, 4), ", ", round(ML_a_U, 4), "] ",
  "& Estimated beta = [", round(ML_b_L, 4), ", ", round(ML_b_U, 4), "] ",
  "& Log-likelihood = [", round(loglike_min, 4), ", ", round(loglike_max, 4), "] ",
  "& AIC = [", round(AIC_min, 4), ", ", round(AIC_max, 4), "] ",
  "& BIC = [", round(BIC_min, 4), ", ", round(BIC_max, 4), "]",
  "& Modified_KS = [", round(DDl[2], 4), ",", round(DDu[2], 4), "]" ,
  "& p-value = [", round(p_value_l, 4), ",", round(p_value_u, 4), "]"
)

print(results)


