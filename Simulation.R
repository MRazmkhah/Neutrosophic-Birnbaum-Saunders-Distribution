# ------------------------------------------------------------
# Simulation Study for Estimation of Birnbaum-Saunders Parameters under Uncertainty
# ------------------------------------------------------------

# Load Required Libraries
library(optimParallel)
library(GoFKernel)

# ------------------------------------------------------------
# Simulation Setup
# ------------------------------------------------------------
N <- 1000  # Number of iterations
n <- 30    # Sample size (must be greater than 1)
a <- 0.25  # True alpha (shape parameter)
b <- 1     # True beta (scale parameter)
e <- 0.1   # Uncertainty noise

# ------------------------------------------------------------
# Function for Parameter Estimation
# ------------------------------------------------------------
f2 <- function(n) {
  # Vectors to store results
  A_L <- c()  # Lower bound of alpha
  A_U <- c()  # Upper bound of alpha
  B_L <- c()  # Lower bound of beta
  B_U <- c()  # Upper bound of beta
  
  # BS CDF
  F_t <- function(t) {
    xx <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
    pnorm(xx, mean = 0, sd = 1)
  }
  
  # Inverse BS CDF
  F_t.inv <- inverse(F_t, lower = 0.001, upper = 500)
  
  # Parallel Cluster Setup
  cl <- makeCluster(4)
  setDefaultCluster(cl = cl)
  
  i <- 0
  repeat {
    i <- i + 1
    
    # Generate Random Variables from BS Distribution
    RV <- function(n) {
      Q <- runif(n)  # Generate uniform quantiles
      lapply(Q, FUN = F_t.inv)  # Apply inverse CDF
    }
    xl <- unlist(RV(n))  # Lower boundary of simulated data
    NN <- runif(n, min = -e, max = e)
    xu <- xl + NN  # Upper boundary of simulated data
    x0 <- (xl + xu) / 2  # Initial point for optimization
    
    # Adjust bounds
    MM <- matrix(c(xl, xu), n, 2)
    L <- apply(MM, 1, min, na.rm = TRUE)
    U <- apply(MM, 1, max, na.rm = TRUE)
    
    # Beta Optimization Function
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
    
    # Optimize Beta
    B_L[i] <- optimParallel(
      par = x0,
      fn = MLb,
      lower = L,
      upper = U,
      hessian = FALSE,
      method = "L-BFGS-B"
    )$value
    B_U[i] <- optimParallel(
      par = x0,
      fn = MLb,
      lower = L,
      upper = U,
      hessian = FALSE,
      method = "L-BFGS-B",
      control = list(fnscale = -1)
    )$value
    
    # Alpha Optimization Function
    MLa <- function(x) {
      n <- length(x)
      r <- n / sum(1 / x)
      s <- sum(x) / n
      DL_b <- function(bb) {
        K <- n / sum(1 / (bb + x))
        bb^2 - bb * (2 * r + K) + r * (s + K)
      }
      b0 <- uniroot(Vectorize(DL_b), lower = r, upper = s)$root
      sqrt(s / b0 + b0 / r - 2)
    }
    
    # Optimize Alpha
    A_L[i] <- optimParallel(
      par = x0,
      fn = MLa,
      lower = L,
      upper = U,
      hessian = FALSE,
      method = "L-BFGS-B"
    )$value
    A_U[i] <- optimParallel(
      par = x0,
      fn = MLa,
      lower = L,
      upper = U,
      hessian = FALSE,
      method = "L-BFGS-B",
      control = list(fnscale = -1)
    )$value
    
    # Stop if required number of iterations is reached
    if (i == N) {
      print("Loop Done")
      break
    }
  }
  
  # Stop Parallel Cluster
  setDefaultCluster(cl = NULL)
  stopCluster(cl)
  
  # Metrics Calculation
  Ave_Bias_a_L <- mean(A_L) - a
  Ave_Bias_a_U <- mean(A_U) - a
  MSE_a_L <- sqrt(mean((A_L - a)^2))
  MSE_a_U <- sqrt(mean((A_U - a)^2))
  Ave_Bias_b_L <- mean(B_L) - b
  Ave_Bias_b_U <- mean(B_U) - b
  MSE_b_L <- sqrt(mean((B_L - b)^2))
  MSE_b_U <- sqrt(mean((B_U - b)^2))
  
  # Results
  c <- round(
    c(mean(A_L), mean(A_U), mean(B_L), mean(B_U),
      Ave_Bias_a_L, Ave_Bias_a_U, Ave_Bias_b_L, Ave_Bias_b_U,
      MSE_a_L, MSE_a_U, MSE_b_L, MSE_b_U), 4
  )
  return(c)
}

# Vectorize the Function for Different Sample Sizes
f3 <- Vectorize(f2)

# ------------------------------------------------------------
# Execution for Different Sample Sizes
# ------------------------------------------------------------
E <- lapply(n, FUN = f3)
for (u in 1:length(n)) {
  d <- E[[u]]
  P <- paste0(
    "NAE of alpha = [", d[1], ",", d[2], "] & ",
    "NAE of beta = [", d[3], ",", d[4], "] & ",
    "NAB of alpha = [", d[5], ",", d[6], "] & ",
    "NAB of beta = [", d[7], ",", d[8], "] & ",
    "NSME of alpha = [", d[9], ",", d[10], "] & ",
    "NSME of beta = [", d[11], ",", d[12], "]"
  )
  print(P)
}
