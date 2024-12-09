# F_N(t): Neutrosophic Cumulative Distribution Function (CDF) of Birnbaum-Saunders
F_N <- function(t, a, b) {
  x <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
  cdf <- pnorm(x, mean = 0, sd = 1)
  return(cdf)
}

# Define values of BS variable
t <- seq(0, 3, 0.01)

# Define neutrosophic parameters directly
a_L <- 0.1
a_U <- 0.35
b_L <- 1
b_U <- 1

# Open a graphical window
windows(3, 3)

# Plot the CDF for the defined values of alpha and beta
y1 <- F_N(t, a_L, b_L)

# Find the maximum values for plotting
ff1 <- function(t) { F_N(t, a_L, b_L) }
ans1 <- optimize(ff1, interval = c(0, 3), maximum = TRUE)

ff2 <- function(t) { F_N(t, a_U, b_U) }
ans2 <- optimize(ff2, interval = c(0, 3), maximum = TRUE)

L <- max(ans1$objective, ans2$objective)

# Create the plot
plot(t, y1, type = "l", ylab = "CDF", xlab = expression(paste(t[N])), ylim = c(0, L))
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)

# Add polygons for shaded regions
for (b in seq(b_L, b_U, 0.05)) {
  for (a in seq(a_L, a_U, 0.05)) {
    y2 <- F_N(t, a, b)
    polygon(c(t, rev(t)), c(y2, rev(y1)), col = "lightblue", border = NA)
    y1 <- y2
  }
}

# Final polygons for boundary cases
y3 <- F_N(t, a_U, b_U)
polygon(c(t, rev(t)), c(y3, rev(y3)), lwd = 1, col = "lightblue")

y4 <- F_N(t, a_L, b_L)
polygon(c(t, rev(t)), c(y4, rev(y4)), lwd = 1, col = "lightblue")
