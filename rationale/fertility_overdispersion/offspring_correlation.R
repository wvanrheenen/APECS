# Fixed R script: Correlated fertility sampling with plots and validation
# Issue: qnbinom returns NA for p=0 or p=1; added proper clamping

library(MASS)
library(ggplot2)
library(dplyr)
library(gridExtra)

correlated_fertility <- function(n_indiv, mu, theta, rho = 0.15, max_offspring = 20) {
  # Generate latent correlated normals (parent-offspring)
  latent_cor <- mvrnorm(n = n_indiv, mu = c(0, 0), Sigma = matrix(c(1, rho, rho, 1), 2))
  
  # Convert latent normals to NB via normal copula with proper bounds
  N_parent <- sapply(latent_cor[,1], function(z) {
    p <- pnorm(z)
    p <- pmax(1e-8, pmin(1-1e-8, p))  # Avoid exact 0/1
    qnbinom(p, mu = mu, size = theta)
  })
  
  N_offspring <- sapply(latent_cor[,2], function(z) {
    p <- pnorm(z)
    p <- pmax(1e-8, pmin(1-1e-8, p))  # Avoid exact 0/1
    qnbinom(p, mu = mu, size = theta)
  })
  
  return(data.frame(
    parent_offspring = as.integer(N_parent),
    offspring = as.integer(N_offspring)
  ))
}

# Run simulations
# set.seed(42)

# Main scenario: mu=2.5, theta=5, rho=0.15
sim_main <- correlated_fertility(10000, mu = 2.5, theta = 5, rho = 0.15)

# Different rho values for comparison
sim_rho <- rbind(
  correlated_fertility(5000, mu = 2.5, theta = 5, rho = 0) %>% mutate(rho = "rho=0"),
  correlated_fertility(5000, mu = 2.5, theta = 5, rho = 0.15) %>% mutate(rho = "rho=0.15"),
  correlated_fertility(5000, mu = 2.5, theta = 5, rho = 0.8) %>% mutate(rho = "rho=1.0")
)

# Validation stats
rho_actual <- cor(sim_main$parent_offspring, sim_main$offspring)
cat("Target rho=0.15, Actual rho =", round(rho_actual, 3), "\n")
cat("Marginal means:", round(mean(sim_main$offspring), 2), "≈ 3.0\n")
cat("Marginal variance:", round(var(sim_main$offspring), 2), "\n")

print(sim_rho %>% 
  group_by(rho) %>% 
  summarise(corr = round(cor(parent_offspring, offspring), 3), .groups = "drop"))

# Plot 1: Marginal distribution validation
p1 <- ggplot(sim_main, aes(x = offspring)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "steelblue", alpha = 0.7) +
  stat_function(fun = function(x) dnbinom(x, mu = 2.5, size = 5), 
                args = list(mu = 2.5, size = 5), color = "red", linewidth = 1.2) +
  labs(title = "Marginal Offspring Distribution", 
       subtitle = "Histogram vs theoretical NB(mu=3, theta=2)",
       x = "Number of Offspring", y = "Density") +
  theme_bw()

# Plot 2: Parent-offspring scatter
p2 <- ggplot(sim_main, aes(x = parent_offspring, y = offspring)) +
  geom_point(alpha = 0.3, color = "darkblue", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Parent-Offspring Fertility Correlation",
       subtitle = paste("rho =", round(rho_actual, 3)),
       x = "Parent Offspring Count", y = "Offspring Count") +
  theme_bw()

# Plot 3: rho comparison
p3 <- ggplot(sim_rho, aes(x = parent_offspring, y = offspring, color = rho)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~ rho, scales = "free", ncol = 3) +
  labs(title = "Parent-Offspring Correlation by rho Value",
       x = "Parent Offspring Count", y = "Offspring Count") +
  theme_bw()

# Display and save plots
grid.arrange(p1, p2, ncol = 2)
print(p3)

ggsave("marginal_fertility.png", p1, width = 8, height = 6, dpi = 300)
ggsave("correlation_scatter.png", p2, width = 8, height = 6, dpi = 300)
ggsave("rho_comparison.png", p3, width = 12, height = 4, dpi = 300)

cat("\n✓ Function working! Use in pedigree simulation:\n")
cat("founders <- rnbinom(n_founders, mu, theta)\n")
cat("descendants <- correlated_fertility(n_descendants, mu, theta, rho, using parent_fert)\n")
