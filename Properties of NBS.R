# Statistical properties
EX <- function(x) {
  a <- x[1]
  b <- x[2]
  b / 2 * (a^2 + 2)
}

VarX <- function(x) {
  a <- x[1]
  b <- x[2]
  (b^2 / 4) * (5 * a^4 + 4 * a^2)
}

CV <- function(x) {
  a <- x[1]
  b <- x[2]
  sqrt(5 * a^4 + 4 * a^2) / (a^2 + 2)
}

sk <- function(x) {
  a <- x[1]
  b <- x[2]
  (44 * a^3 + 24 * a) / (5 * a^2 + 4)^1.5
}

ku <- function(x) {
  a <- x[1]
  b <- x[2]
  3 + (558 * a^4 + 240 * a^2) / (5 * a^2 + 4)^2
}

# Define matrices for parameter boundaries
A <- matrix(c(
  0.1, 0.35, 0.5, 0.75, 1, 1.5, 2, 3,
  0.1, 0.35, 0.5, 0.75, 1, 1.5, 2, 3,
  0.1, 0.35, 0.5, 0.75, 1, 1.5, 2, 3
), ncol = 2, byrow = TRUE)

B <- matrix(c(
  0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1,
  1, 2, 1, 2, 1, 2, 1, 2,
  2, 3, 2, 3, 2, 3, 2, 3
), ncol = 2, byrow = TRUE)

# Initialize results matrix
D <- matrix(NA, nrow = 12, ncol = 1)

# Loop through rows of A and B to compute statistical properties
for (j in 1:nrow(A)) {
  a_L <- A[j, 1]
  a_U <- A[j, 2]
  b_L <- B[j, 1]
  b_U <- B[j, 2]
  
  # Midpoints for optimization
  a0 <- (a_L + a_U) / 2
  b0 <- (b_L + b_U) / 2
  
  # Mean
  mu_L <- round(optim(par = c(a0, b0), fn = EX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B")$value, 3)
  mu_U <- round(optim(par = c(a0, b0), fn = EX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)
  
  # Variance
  var_L <- round(optim(par = c(a0, b0), fn = VarX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B")$value, 3)
  var_U <- round(optim(par = c(a0, b0), fn = VarX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)
  
  # Coefficient of Variation
  CV_L <- round(optim(par = c(a0, b0), fn = CV, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B")$value, 3)
  CV_U <- round(optim(par = c(a0, b0), fn = CV, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)
  
  # Skewness
  skew_L <- round(optim(par = c(a0, b0), fn = sk, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                        method = "L-BFGS-B")$value, 3)
  skew_U <- round(optim(par = c(a0, b0), fn = sk, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                        method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)
  
  # Kurtosis
  kur_L <- round(optim(par = c(a0, b0), fn = ku, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B")$value, 3)
  kur_U <- round(optim(par = c(a0, b0), fn = ku, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)
  
  # Populate results
  D[j, ] <- paste0("[", b_L, ", ", b_U, "] & [", a_L, ", ", a_U, "] & ",
                   "[", mu_L, ", ", mu_U, "] & [", var_L, ", ", var_U, "] & ",
                   "[", CV_L, ", ", CV_U, "] & [", skew_L, ", ", skew_U, "] & ",
                   "[", kur_L, ", ", kur_U, "]")
}

# Print the results
print(D)


