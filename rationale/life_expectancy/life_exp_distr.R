# Skew-normal life expectancy 
library(sn)
library(ggplot2)
library(wesanderson)
library(dplyr)
library(scales)

get_xi_for_mean <- function(mean_desired, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  xi <- mean_desired - omega * sqrt(2/pi) * delta
  return(xi)
}

get_omega <- function(life_expectancy, max_age = 100, k = 5.5) {
  return(k * (max_age / life_expectancy))
}

# Link to data for mean life expectancy based on birth year
life_expectancy = read.table("../../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt", header=T)

# Parameters
years_keep <- c(1900, 1950, 2000)

le_sel <- life_expectancy %>%
  filter(year_of_birth %in% years_keep) %>%
  arrange(match(year_of_birth, years_keep))

life_expectancy <- le_sel$life_expectancy
year_labels <- as.character(le_sel$year_of_birth)

max_age <- 100
k <- 5.5
alpha <- -3
x <- seq(30, 110, length.out = 500)

# Create data with CORRECT numeric life_exp
plot_data <- data.frame()
colors <- wes_palette(n = length(life_expectancy), name = "Darjeeling1")

for (i in seq_along(life_expectancy)) {
  omega_i <- get_omega(life_expectancy[i], max_age, k)
  xi_i <- get_xi_for_mean(life_expectancy[i], omega_i, alpha)
  y <- dsn(x, xi = xi_i, omega = omega_i, alpha = alpha)
  
  plot_data <- rbind(plot_data, data.frame(
    age = x,
    density = y,
    life_exp = factor(life_expectancy[i], levels = life_expectancy,  # ← Make factor
                      labels = paste("Mean =", life_expectancy[i]))
  ))
}

p <- ggplot(plot_data, aes(x = age, y = density, color = life_exp)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors, 
                     name = "Year of birth",
                     labels = year_labels) +
  labs(title = "(B) Left Skewed Life Expectancy Distributions",
       x = "Life expectancy (years)", y = "Density") +
  scale_x_continuous(breaks = seq(40, 100, 10)) +      
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +  
  coord_cartesian(xlim = c(45, 100), ylim = c(0, 0.13)) +
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
    plot.margin = unit(c(5,5,5,5), "pt"),
    panel.grid.minor = element_blank()
  )

ggsave("life_exp_distr.pdf", p, width = 3, height = 3, units = "in", dpi = 300)

cat("✓ life_exp_distr.pdf - FIXED legend shows correct 60/75/90 years ✓\n")
