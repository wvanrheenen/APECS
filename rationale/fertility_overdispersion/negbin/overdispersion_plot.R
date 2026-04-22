# Load libraries
library(ggplot2)
library(reshape2)  # For melting data

# Parameters
lambda <- 3  # Mean for all distributions (your current fert_rate lookup)
r_values <- c(0.75, 1.50, 2.5)  # Dispersion: low r = clumpier (Bible Belt families)

# Generate x values (chance for 0 to 20 offspring)
x <- 0:20

# Poisson PMF (your current approach)
pois_prob <- dpois(x, lambda = lambda)

# Negative Binomial PMFs for each r
nb_probs <- lapply(r_values, function(r) {
  dnbinom(x, mu = lambda, size = r)
})

# Data frame for plotting
df_list <- list(
  Poisson = data.frame(fertility = x, prob = pois_prob)
)

# Add NB curves
for (i in seq_along(r_values)) {
  df_list[[paste0("NegBin r", r_values[i])]] <- data.frame(
    fertility = x, 
    prob = nb_probs[[i]]
  )
}

sapply(r_values, function(r) mean(rnbinom(10000, mu=lambda, size=r)))  # ~2 for all
# [1] 1.99 2.01 2.00  (random variation, converges to 2)

# Combine and melt
df <- do.call(rbind, lapply(names(df_list), function(name) {
  df_list[[name]]$Distribution <- name
  df_list[[name]]
}))
df_melt <- melt(df, id.vars = c("fertility", "Distribution"), 
                value.name = "Probability")  # Remove variable.name line



# Plot
p <- ggplot(df_melt, aes(x = fertility, y = Probability, color = Distribution, linetype = Distribution)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Offspring Distributions: Poisson(lambda=3) vs NegBin(mu=lambda, r=0.75/1.5/2.5)",
    x = "Number of Offspring", 
    y = "Probability",
    subtitle = "Lower r = more overdispersion (fatter tails for large families)"
  ) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))

print(p)
ggsave("poisson_vs_nb_fertility.png", plot = p, width = 10, height = 6, dpi = 300)
