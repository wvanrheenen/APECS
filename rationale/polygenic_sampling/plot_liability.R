#!/usr/bin/env Rscript

# Liability threshold model visualization Rscript
#
# Visualizes phenotype liability P = G + E,
# where G ~ N(0, h2), E ~ N(0, 1-h2), and P ~ N(0, 1).
# Shows how exceeding the threshold (based on K=prevalence)
# determines case/control phenotypes.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
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

# Prepare data frame
data <- data.frame(G=G, E=E, P=P, Phenotype=Phenotype)

p1 <- ggplot(data, aes(x=G)) +
  geom_histogram(bins=100, fill=wes_palette("Darjeeling1")[2], alpha=0.7) +
  labs(title="Genetic component G", subtitle="~ N(0,h2)", x="G", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

p2 <- ggplot(data, aes(x=E)) +
  geom_histogram(bins=100, fill=wes_palette("Darjeeling1")[3], alpha=0.7) +
  labs(title="Environmental component E", subtitle="~ N(0,(1-h2))", x="E", y="Count") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )


# Plot phenotypic liability distribution with threshold
# Base histogram plot of P colored by Phenotype
p3_base <- ggplot(data, aes(x=P, fill=Phenotype)) +
  geom_histogram(bins=100, alpha=0.7, position="identity") 

max_count <- max(ggplot_build(p3_base)$data[[1]]$count)

p3 <- p3_base +
  geom_vline(xintercept=threshold, color="red", linetype="dashed") +
  stat_function(fun = dnorm, args = list(mean=0, sd=1), geom = "area",
                xlim = c(threshold, 4), fill = "red", alpha = 0.2) +
  scale_fill_manual(values=c("Unaffected"="gray", "Affected"="red")) +
  labs(title="Phenotype status by Phenotypic liability P",
       subtitle = "Phenotypic liability = Genetic (G) + Non-genetic (E)",
       x="P", y="Count") +
  theme_minimal() +
  theme(legend.position = "bottom",
  plot.title = element_text(size = 12, hjust = 0.5),        
  plot.subtitle = element_text(size = 10, hjust = 0.5)  
    ) +
  annotate("text", x=threshold + 1.2, y=max_count * 0.6, label=sprintf("Liability\nThreshold \nLT = %.2f", threshold),
           color="black", angle=0, vjust=-0.5, size=4) +
  annotate("text", x=threshold + 1.5, y=max_count * 0.15, label=paste0("Lifetime risk\nK = ", scales::percent(K)), color="black", size=4)


# Arrange plots in grid
# Create blank spacer grob
blank <- rectGrob(gp = gpar(col = NA, fill = NA))

# List of grobs: p1, p2 for left col; blank spacer in middle; p3 for full right col
grobs_list <- list(p1, p2, blank, p3)

# Layout matrix:
# Columns: 1 = left plots; 
#          2 = blank spacer;
#          3 = full-height right (p3)
layout_mat <- rbind(
  c(1, 3, 4),
  c(2, 3, 4)
)

# Arrange plots with widths setting spacer narrower
pdf("polygenic_model.pdf", width=6, height=7)
grid.arrange(grobs = grobs_list,
             layout_matrix = layout_mat,
             widths = c(1, 0.2, 2.1))
dev.off()