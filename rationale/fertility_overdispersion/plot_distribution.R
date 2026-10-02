library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)
library(RColorBrewer)

# Link to data for mean fert rate based on birth year
fert_rate = read.table("../../data/fertility_rate/Gapminder/GM_fertility_rate_Netherlands_1800_2100.txt", header=T)

# Parameters
years_keep <- c(1900, 1950, 2000)

fr_sel <- fert_rate %>%
  filter(year %in% years_keep) %>%
  arrange(match(year, years_keep))

fertility_rate <- fr_sel$mean_fertility
year_labels <- as.character(fr_sel$year)

max_offspring <- 15
sizes <- 5

# Create NB data
plot_data <- data.frame()
# First three colours from the four-colour ColorBrewer RdBu palette
colors <- RColorBrewer::brewer.pal(4, "RdBu")[c(1, 3, 4)]

for (i in seq_along(fertility_rate)) {
  m <- fertility_rate[i]
  offspring <- 0:max_offspring
  prob <- dnbinom(offspring, size = sizes, mu = m)

  temp_df <- data.frame(
    offspring = offspring,
    mean = m,
    year = year_labels[i],
    prob = prob
  )

  plot_data <- rbind(plot_data, temp_df)
}

plot_data <- plot_data %>%
  mutate(
    year_label = factor(year, levels = year_labels)
  ) %>%
  filter(prob > 1e-4)

# Plot lines for each fertility year
p <- ggplot(plot_data, aes(x = offspring, y = prob, color = year_label)) +
  geom_line(size = 1, alpha = 0.9) +
  scale_color_manual(
    values = colors,
    name = "Year of birth",
    labels = year_labels
  ) +
  scale_y_continuous(
    breaks = c(0, 0.1, 0.2, 0.3),
    labels = c("0", "0.1", "0.2", "0.3")
  ) +
  labs(
    title = "(B) Fertility Rate Distributions",
    x = "Number of Offspring",
    y = "Density"
  ) +
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

ggsave("nb_fertility_distributions.pdf", p, width = 3, height = 3, units = "in", dpi = 300)

cat("✓ nb_fertility_distributions.pdf - FIXED legend shows correct years ✓\n")