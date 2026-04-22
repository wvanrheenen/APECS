#!/usr/bin/env Rscript
options(warn = 2)

# Snakemake integration
if ("snakemake" %in% ls()) {
  input_tsv <- snakemake@input[[1]]
  output_pdf <- snakemake@output[[1]]
} else {
  input_tsv <- "monogenic_inheritance.tsv"
  output_pdf <- "monogenic_prob.pdf"
}

library(dplyr)
library(ggplot2)
library(wesanderson)
library(readr)
library(scales)

# Read data
df <- read_tsv(input_tsv)

# Filter and prepare (match Snakefile exactly)
df <- df %>% 
  filter(penetrance %in% c("0.2", "0.5", "0.8")) %>%
  mutate(
    penetrance = factor(penetrance, levels = c("0.2", "0.5", "0.8")),
    lambda = factor(lambda, levels = c("1.5", "3")),
    generation = factor(generation),
    simulated_proportion = as.numeric(simulated_proportion),
    expected_proportion = as.numeric(expected_proportion)
  )

# Expected segments
expected_segments <- df %>%
  group_by(penetrance, lambda, generation) %>%
  summarise(expected_proportion = first(expected_proportion), .groups = "drop") %>%
  mutate(
    generation_num = as.numeric(as.character(generation)),
    x = generation_num - 0.5,
    xend = generation_num + 0.5,
    y = expected_proportion,
    yend = expected_proportion
  )

# Labels
penetrance_labels <- setNames(
  paste0("Penetrance: ", percent(as.numeric(levels(df$penetrance)), accuracy = 0.1)),
  levels(df$penetrance)
)

# Darjeeling1 palette (matching disease_onset.R)
darjeeling_cols <- wes_palette("Darjeeling1", 2)

# Plot with disease_onset.R matching theme
p <- ggplot(df, aes(x = generation, y = simulated_proportion, fill = lambda)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, 
               outlier.shape = 16, outlier.size = 1) +
  # geom_segment(
  #   data = expected_segments,
  #   aes(x = x, xend = xend, y = y, yend = yend),
  #   color = "black", linewidth = 0.8, linetype = "dashed",
  #   inherit.aes = FALSE
  # ) +
  facet_wrap(~ penetrance, scales = "free_x", ncol = 3,
             labeller = labeller(penetrance = penetrance_labels)) +
  labs(
    title = "(A) Monogenic Inheritance Validation",
    x = "Generations after founder", 
    y = "Proportion affected offspring",
    fill = "Fertility rate"
  ) +
  scale_fill_manual(values = darjeeling_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, NA)) +
  coord_cartesian(xlim = c(0.5, max(as.numeric(levels(df$generation))) + 0.5)) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.margin = margin(b = 2),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(output_pdf, p, width = 9, height = 3, units = "in", dpi = 300)

cat("Monogenic inheritance plot saved as 9x3 PDF\n")
cat("Darjeeling1 palette: λ=1.5 =", darjeeling_cols[1], ", λ=2 =", darjeeling_cols[2], "\n")
