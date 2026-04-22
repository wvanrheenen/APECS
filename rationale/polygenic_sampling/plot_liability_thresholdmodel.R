#!/usr/bin/env Rscript

# Liability threshold model visualization - Phenotypic liability only (p3)
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(wesanderson)
})

set.seed(123)

# Parameters
h2 <- 0.5       # heritability
K <- 0.05       # lifetime prevalence (e.g., 5%)
threshold <- qnorm(1 - K)  # liability threshold from prevalence
n <- 25000      # sample size for simulation

# Simulate genetic and environmental components
G <- rnorm(n, mean=0, sd=sqrt(h2))
E <- rnorm(n, mean=0, sd=sqrt(1 - h2))
P <- G + E  # Phenotypic liability

# Assign phenotype status based on threshold
Phenotype <- ifelse(P > threshold, "Affected", "Unaffected")

# Create data frame
data <- data.frame(P = P, Phenotype = Phenotype)

# Single plot: Phenotypic liability distribution with threshold (original p3 style)
p3 <- ggplot(data, aes(x = P, fill = Phenotype)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), geom = "area",
                xlim = c(threshold, 4), fill = "red", alpha = 0.2) +
  scale_fill_manual(values = c("Unaffected" = "gray", "Affected" = "red")) +
  scale_y_continuous(limits = c(0, 900), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Liability Threshold Model: Phenotypic Liability P ~ N(0,1)",
    subtitle = paste0("P = G + E | h² = ", h2, " | Prevalence K = ", percent(K), " | Threshold LT = ", round(threshold, 2)),
    x = "Phenotypic liability", 
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  ) +
  annotate("text", x = threshold + 1.2, y = 540, 
           label = sprintf("Liability\nThreshold \nLT = %.2f", threshold),
           color = "black", angle = 0, vjust = -0.5, size = 4) +
  annotate("text", x = threshold + 1.5, y = 135, 
           label = paste0("Lifetime risk\nK = ", scales::percent(K)), 
           color = "black", size = 4)

# Save single plot
ggsave("liability_threshold_model.pdf", p3, width = 8, height = 6, dpi = 300, 
       device = "pdf", bg = "white")
