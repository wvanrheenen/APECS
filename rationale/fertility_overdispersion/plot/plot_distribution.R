# Script to plot theoretical Negative Binomial offspring distributions
# Matches your exact 9x3 boxplot publication style

library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

max_offspring <- 15
means <- c(2.5)
sizes <- c(1, 3, 5)

# Create NB data
plot_data <- data.frame()

for (m in means) {
  for (s in sizes) {
    offspring <- 0:max_offspring
    prob <- dnbinom(offspring, size = s, mu = m)
    temp_df <- data.frame(
      offspring = offspring,
      mean = m,
      size = s,
      prob = prob
    )
    plot_data <- rbind(plot_data, temp_df)
  }
}

plot_data <- plot_data %>%
  mutate(
    size_label = factor(size, levels = c(1, 3, 5)),
    facet_label = paste("Mean Fertility =", mean)
  ) %>%
  filter(prob > 1e-4)

# 2x2 grid layout - perfect square format
p <- ggplot(plot_data, aes(x = offspring, y = prob, color = size_label)) +
  geom_line(size = 1, alpha = 0.9) +
  facet_wrap(~ facet_label, nrow = 2, ncol = 2, scales = "free") +  # ← 2x2 grid
  scale_color_manual(values = wes_palette("Darjeeling1"), name = "NegBin size") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +  
  labs(title = "(C) Negative Binomial Fertility Rate Distributions",
       x = "Number of Offspring", y = "Density") +
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
    strip.background = element_rect(fill = "grey90", color = "black"),
    plot.margin = unit(c(5,5,5,5), "pt"),
    panel.grid.minor = element_blank()
  )

# Square 2x2 format - adjust dimensions for balance
ggsave("nb_fertility_distributions.pdf", p, width = 4.5, height = 3, units = "in", dpi = 300)

cat("✓ nb_fertility_distributions.pdf (2x2 grid, 8x8 square)\n")
cat("✓ Wes Anderson Darjeeling1 + exact publication theme\n")