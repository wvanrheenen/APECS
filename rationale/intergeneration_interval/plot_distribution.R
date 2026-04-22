# Load required libraries (add these if not in libraries_simPed.R)
library(ggplot2)
library(wesanderson)
library(dplyr)
library(scales)

source("../../src/libraries_simPed.R")
set.seed(123)

generate_gen_yr <- function(mean_gen_yr, min_gen_yr = (mean_gen_yr-5), max_gen_yr = (mean_gen_yr+5)) {
  gen_yr <- rnorm(1, mean = mean_gen_yr, sd = (max_gen_yr - min_gen_yr) / 6)
  return((max(min_gen_yr, min(max_gen_yr, gen_yr))))
}

# Generate simulated data
n_sim <- 2500
mean_gen_yr <- 30
simulated_intervals <- data.frame(
  yob = 1950,  # Fixed cohort for demonstration
  simulated_interval = replicate(n_sim, generate_gen_yr(mean_gen_yr))
)

# Calculate summary statistics for plotting
interval_estimates <- simulated_intervals %>%
  group_by(yob) %>%
  summarise(
    simulated_mean = mean(simulated_interval),
    simulated_lower = quantile(simulated_interval, 0.025),
    simulated_upper = quantile(simulated_interval, 0.975),
    .groups = 'drop'
  )

# Calculate and print 1 SD of the inter-generation interval distribution
interval_sd <- sd(simulated_intervals$simulated_interval, na.rm = TRUE)
cat(sprintf("Standard deviation of inter-generation intervals: %.2f years\n", interval_sd))
cat(sprintf("Mean: %.2f, SD: %.2f (95%% range: %.0f-%.0f)\n", 
            mean(simulated_intervals$simulated_interval),
            interval_sd,
            quantile(simulated_intervals$simulated_interval, 0.025),
            quantile(simulated_intervals$simulated_interval, 0.975)))

# Summary table
summary_stats <- simulated_intervals %>%
  summarise(
    mean = mean(simulated_interval),
    sd = sd(simulated_interval),
    min = min(simulated_interval),
    max = max(simulated_interval),
    n = n()
  )
print("Full summary:")
print(summary_stats)

# # Save the data
# Define Darjeeling1 palette (MISSING in your code)
darjeeling_cols <- wes_palette("Darjeeling1", 5)


# Save plot to PNG (same style as your life expectancy plot)
p <- ggplot(simulated_intervals, aes(x = simulated_interval)) +
  geom_density(linewidth = 1.2, fill = darjeeling_cols[1], alpha = 0.3, color = darjeeling_cols[1]) +
  labs(title = "Inter-generation Interval Distribution",
       x = "Inter-generation Interval (years)", 
       y = "Probability") +
  scale_x_continuous(limits = c(22, 38)) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 5), face = "bold"),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.margin = margin(b = 2),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = unit(c(5,5,5,5), "pt")
  )

# Save EXACTLY like your disease plot (3x3 PDF)
ggsave("intergen_interval_density.pdf", p, 
       width = 3, height = 3, units = "in", dpi = 300)

cat("Inter-generation density plot saved as 9x3 PDF (exact disease plot match)\n")
cat("Darjeeling1 palette used - density line with probability y-axis\n")

# Save raw simulation data
write.csv(simulated_intervals, "intergen_interval_raw_data.csv", row.names = FALSE)
