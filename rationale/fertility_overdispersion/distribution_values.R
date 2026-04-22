# Script to plot theoretical Poisson vs Negative Binomial offspring distributions
# Fixed: no rowwise/reframe - simple loop approach

library(ggplot2)
library(dplyr)
library(tidyr)

max_offspring <- 20
means <- c(1, 2, 3, 4)
sizes <- c(1, 2, 3, 4, 5)

# Create all combinations using nested loops
plot_data <- data.frame()

for (m in means) {
  for (s in c(sizes, Inf)) {
    offspring <- 0:max_offspring
    if (is.infinite(s)) {
      # Poisson
      prob <- dpois(offspring, lambda = m)
      dist_type <- "Poisson"
    } else {
      # Negative Binomial
      prob <- dnbinom(offspring, size = s, mu = m)
      dist_type <- "Negative Binomial"
    }
    temp_df <- data.frame(
      offspring = offspring,
      mean = m,
      size = s,
      dist = dist_type,
      prob = prob
    )
    plot_data <- rbind(plot_data, temp_df)
  }
}

plot_data <- plot_data %>%
  mutate(
    size_label = ifelse(is.infinite(size), "inf (Poisson)", as.character(size))
  ) %>%
  filter(prob > 1e-4)

# Plot
p <- ggplot(plot_data, aes(x = offspring, y = prob, 
                          color = factor(size_label), 
                          linetype = dist)) +
  geom_line(size = 1.2, alpha = 0.9) +
  facet_wrap(~ paste("Fertility Rate =", mean), scales = "free_y", ncol = 2) +
  scale_color_viridis_d(name = "NB size (theta)", option = "plasma") +
  labs(title = "Theoretical Offspring Distributions: Poisson vs Negative Binomial",
       subtitle = "Higher theta → closer to Poisson.",
       x = "Number of Offspring", y = "Probability") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
ggsave("fertility_distributions_theoretical.png", p, width = 12, height = 8, dpi = 300)
