# Install and load the skew-normal package if not already installed
if (!require("sn")) install.packages("sn")
library(sn)

# Function to compute the correct xi for a given desired mean
get_xi_for_mean <- function(mean_desired, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  xi <- mean_desired - omega * sqrt(2/pi) * delta
  return(xi)
}

# Function to compute omega based on life expectancy
get_omega <- function(life_expectancy, max_age = 100, k = 6) {
  return(k * (max_age / life_expectancy))
}

# Desired mean (life expectancy)
life_expectancy <- c(60, 75, 90)
max_age <- 100
k <- 6
alpha <- -3

# Age range for plotting
x <- seq(20, 110, length.out = 500)

# Colors for plotting
colors <- c("blue", "red", "green")

# Open PDF device
pdf("skew_test.pdf", width = 8, height = 5.5)

plot(NULL, xlim = range(x), ylim = c(0, 0.18),
     xlab = "Age", ylab = "Density",
     main = paste0("Distribution of life expectancy for mean expectancy: ", paste(life_expectancy, collapse = ", ")))

legend_labels <- character(0)

for (i in seq_along(life_expectancy)) {
  omega_i <- get_omega(life_expectancy[i], max_age, k)
  xi_i <- get_xi_for_mean(life_expectancy[i], omega_i, alpha)
  y <- dsn(x, xi = xi_i, omega = omega_i, alpha = alpha)
  lines(x, y, col = colors[i], lwd = 2, lty = 1) # <-- all solid lines
  legend_labels <- c(legend_labels, paste0("mean=", life_expectancy[i], ", omega=", round(omega_i,2), ", alpha=", alpha))
}

legend("topright", legend = legend_labels, col = colors[1:length(life_expectancy)],
       lwd = 2, lty = 1, cex=0.8) # <-- all solid lines

dev.off()
