# H_N(t): Neutrosophic Hazard Function for Birnbaum-Saunders
H_N <- function(t, a, b) { 
  u <- (1 / sqrt(2 * pi)) * exp(-1 / (2 * (a^2)) * (t / b + b / t - 2)) * 
    (t^(-3 / 2) * (t + b)) / (2 * a * sqrt(b))
  x <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
  hazard <- u / (1 - pnorm(x, mean = 0, sd = 1))
  return(hazard)
}

# Define values of BS variable
t <- seq(0, 2, 0.01)

# Define neutrosophic parameter bounds
a_L <- 0.1  # Lower bound of alpha
a_U <- 0.35 # Upper bound of alpha
b_L <- 1    # Lower bound of beta
b_U <- 1    # Upper bound of beta

# Open a graphical window
windows(3, 3)

# Compute the initial hazard function
y1 <- H_N(t, a_L, b_L)

# Find the maximum hazard value for scaling the plot
ff1 <- function(t) { H_N(t, a_L, b_L) }
ans1 <- optimize(ff1, interval = c(0, 2), maximum = TRUE)

ff2 <- function(t) { H_N(t, a_U, b_U) }
ans2 <- optimize(ff2, interval = c(0, 2), maximum = TRUE)

L <- max(ans1$objective, ans2$objective)

# Plot the hazard function
plot(t, y1, type = "l", ylab = "Hazard Function", xlab = expression(paste(t[N])), ylim = c(0, L))
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)

# Add polygons for shaded regions
for (b in seq(b_L, b_U, 0.05)) {
  for (a in seq(a_L, a_U, 0.05)) {
    y2 <- H_N(t, a, b)
    polygon(c(t, rev(t)), c(y2, rev(y1)), col = "lightblue", border = NA)
    y1 <- y2
  }
}

# Final polygons for boundary cases
y3 <- H_N(t, a_U, b_U)
polygon(c(t, rev(t)), c(y3, rev(y3)), lwd = 1, col = "lightblue")

y4 <- H_N(t, a_L, b_L)
polygon(c(t, rev(t)), c(y4, rev(y4)), lwd = 1, col = "lightblue")
