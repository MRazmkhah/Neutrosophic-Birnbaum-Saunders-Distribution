# Probability Density Function (PDF) of the Neutrosophic Birnbaum-Saunders Distribution
f_N <- function(t, a, b) { 
  # Compute the PDF
  pdf <- (1 / sqrt(2 * pi)) * 
    exp(-1 / (2 * (a^2)) * (t / b + b / t - 2)) * 
    (t^(-3/2) * (t + b)) / (2 * a * sqrt(b))
  return(pdf)
}

# Define BS variable values and neutrosophic parameters
t <- seq(0, 3, 0.01)  # BS variable values
a_L <- 0.25           # Lower bound of neutrosophic alpha
a_U <- 0.25           # Upper bound of neutrosophic alpha
b_L <- 0.5            # Lower bound of neutrosophic beta
b_U <- 1              # Upper bound of neutrosophic beta

# Open a graphical window
windows(3, 3)

# Plot the first PDF with alpha = a_L and beta = b_L
y1 <- f_N(t, a_L, b_L)
plot(t, y1, type = "l", 
     ylab = "PDF", xlab = expression(paste(t[N])), 
     lwd = 1, ylim = c(0, 5))
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)

# Add interval-based neutrosophic parameters and fill areas
for (b in seq(b_L, b_U, 0.01)) {
  for (a in seq(a_L, a_U, 0.01)) {
    # Compute the next PDF
    y2 <- f_N(t, a, b)
    # Add shaded areas between successive curves
    polygon(c(t, rev(t)), c(y2, rev(y1)), col = "#6BD7AF", border = NA)
    y1 <- y2  # Update y1 for the next iteration
  }
}

# Final updates to the plot with boundaries
y3 <- f_N(t, a_U, b_U)
polygon(c(t, rev(t)), c(y3, rev(y1)), lwd = 1, col = "#6BD7AF")
y4 <- f_N(t, a_L, b_L)
polygon(c(t, rev(t)), c(y4, rev(y4)), lwd = 1, col = "#6BD7AF")


