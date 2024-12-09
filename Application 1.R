# ------------------------------------------------------------------------
# Neutrosophic Parameter Estimation and Inferences for Birnbaum-Saunders Model (Real Dataset 1)
# ------------------------------------------------------------------------

# Load Required Libraries
library(GoFKernel)

# ------------------------------------------------------------------------
# Setup for Real Data
# ------------------------------------------------------------------------
e <- 0.001  # Indeterminacy noise
real_data_BS <- c(70, 90, 96, 97, 99, 100, 103, 104, 104, 105, 107, 108, 108, 108, 109, 109, 112, 112, 113, 114, 114, 114,
                  116, 119, 120, 120, 120, 121, 121, 123, 124, 124, 124, 124, 124, 128, 128, 129, 129, 130, 130, 130, 131, 131,
                  131, 131, 131, 132, 132, 132, 133, 134, 134, 134, 134, 134, 136, 136, 137, 138, 138, 138, 139, 139, 141, 141,
                  142, 142, 142, 142, 142, 142, 144, 144, 145, 146, 148, 148, 149, 151, 151, 152, 155, 156, 157, 157, 157, 157,
                  158, 159, 162, 163, 163, 164, 166, 166, 168, 170, 174, 196, 212)

n <- length(real_data_BS)

# Define Neutrosophic Dataset
L <- real_data_BS - e  # Lower bound
U <- real_data_BS + e  # Upper bound
x0 <- (L + U) / 2      # Initial optimization point

# ------------------------------------------------------------------------
# Beta Optimization Function
# ------------------------------------------------------------------------
MLb <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- sum(x) / n
  DL_b <- function(bb) {
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }
  h <- uniroot(Vectorize(DL_b), lower = r, upper = s)$root
  return(h)
}

ML_b_L <- optim(par = x0, fn = MLb, lower = L, upper = U, hessian = FALSE, method = "L-BFGS-B")$value
ML_b_U <- optim(par = x0, fn = MLb, lower = L, upper = U, hessian = FALSE, method = "L-BFGS-B", control = list(fnscale = -1))$value

# ------------------------------------------------------------------------
# Alpha Optimization Function
# ------------------------------------------------------------------------
MLa <- function(x) {
  r <- n / sum(1 / x)
  s <- sum(x) / n
  DL_b <- function(bb) {
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }
  b0 <- uniroot(Vectorize(DL_b), lower = r, upper = s)$root
  sqrt(s / b0 + b0 / r - 2)
}

ML_a_L <- optim(par = x0, fn = MLa, lower = L, upper = U, hessian = FALSE, method = "L-BFGS-B")$value
ML_a_U <- optim(par = x0, fn = MLa, lower = L, upper = U, hessian = FALSE, method = "L-BFGS-B", control = list(fnscale = -1))$value

# ------------------------------------------------------------------------
# Log-Likelihood Function for Real Data
# ------------------------------------------------------------------------
loglikelihood <- function(x) {
  r <- n / sum(1 / x)
  s <- sum(x) / n
  DL_b <- function(bb) {
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }
  beta_hat <- uniroot(Vectorize(DL_b), lower = r, upper = s)$root
  alpha_hat <- sqrt(s / beta_hat + beta_hat / r - 2)
  -n * log(2) - n * log(alpha_hat) - n * log(beta_hat) - n * log(2 * pi) / 2 +
    sum(log(sqrt(beta_hat / x) + (beta_hat / x)^1.5)) -
    (1 / (2 * alpha_hat^2)) * sum(x / beta_hat + beta_hat / x - 2)
}

loglike_min <- optim(par = x0, fn = loglikelihood, lower = L, upper = U, hessian = FALSE, method = "L-BFGS-B")$value
loglike_max <- optim(par = x0, fn = loglikelihood, lower = L, upper = U, hessian = FALSE, method = "L-BFGS-B", control = list(fnscale = -1))$value


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

DDl <- Modified_KS(L, ML_a_L, ML_b_L)
DDu <- Modified_KS(U, ML_a_U, ML_b_U)

# ------------------------------------------------------------------------
# Output Results for Real Data
# ------------------------------------------------------------------------
set.seed(88888888)
result <- paste0(
  "Estimated alpha = [", round(ML_a_L, 4), ", ", round(ML_a_U, 4), "] ",
  "& Estimated beta = [", round(ML_b_L, 4), ", ", round(ML_b_U, 4), "] ",
  "& Log-likelihood = [", round(loglike_min, 4), ", ", round(loglike_max, 4), "] ",
  "& AIC = [", round(AIC_min, 4), ", ", round(AIC_max, 4), "] ",
  "& BIC = [", round(BIC_min, 4), ", ", round(BIC_max, 4), "]",
  "& Modified_KS = [", round(DDl[2], 4), ",", round(DDu[2], 4), "]" 
)

print(result)







