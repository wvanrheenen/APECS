### Polygenic inheritance rationale plotting script ###

# Load libraries
options(warn = 2)
library(readr)
library(dplyr)
library(ggplot2)

# Read the polygenic output
df <- read_tsv("polygenic_K_h2.out")

# Make sure columns are appropriate
df <- df %>%
  mutate(
    h2_als = factor(h2_als),       # treat heritability as categorical
    K_als = as.factor(K_als),      # input prevalence categories for x-axis
    rg     = factor(rg)            # facetting variable
  )

# # Plot: input K on x, simulated K on y, boxplots grouped by h², facetted by rg
# p <- ggplot(df, aes(x = K_als, y = K_empirical, fill = h2_als)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   # Add theoretical prevalence values as reference points (black diamonds)
#   geom_point(aes(y = as.numeric(as.character(K_als))),
#              color = "black", shape = 18, size = 2,
#              position = position_dodge(width = 0.8)) +
#   labs(
#     title = "Simulated vs Expected Prevalence of ALS",
#     subtitle = "Polygenic liability-threshold simulations across genetic correlations",
#     x = "Input K_ALS (prevalence)",
#     y = "Empirical prevalence from simulations",
#     fill = "Heritability (h²)"
#   ) +
#   facet_wrap(~ rg, labeller = label_both) +   # facet on genetic correlation
#   theme_minimal()

# print(p)

# # Plot: input h² on x, empirical h² on y, boxplots grouped by K_ALS, facetted by rg
# p <- ggplot(df, aes(x = h2_als, y = h2_empirical, fill = K_als)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   # Add theoretical input h2 as black diamonds for reference
#   geom_point(aes(y = as.numeric(as.character(h2_als))),
#              color = "black", shape = 18, size = 2,
#              position = position_dodge(width = 0.8)) +
#   labs(
#     title    = "Simulated vs Expected Heritability of ALS",
#     subtitle = "Polygenic liability-threshold simulations across genetic correlations",
#     x        = "Input heritability (h²_ALS)",
#     y        = "Empirical heritability from simulations",
#     fill     = "Input prevalence (K_ALS)"
#   ) +
#   facet_wrap(~ rg, labeller = label_both) +   # facet on genetic correlation
#   theme_minimal()

# print(p)

# Plot facetted by lambda for simulated vs expected prevalence
p_lambda <- ggplot(df, aes(x = K_als, y = K_empirical, fill = h2_als)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(aes(y = as.numeric(as.character(K_als))),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.8)) +
  labs(
    title = "Simulated vs Expected Prevalence of ALS",
    subtitle = "Polygenic liability-threshold simulations across lambda values",
    x = "Input K_ALS (prevalence)",
    y = "Empirical prevalence from simulations",
    fill = "Heritability (h²)"
  ) +
  facet_wrap(~ lambda, scales = "free_x", labeller = labeller(lambda = function(x) paste0("Fertility rate: ", x))) +
  theme_minimal()
print(p_lambda)

# Plot facetted by lambda for simulated vs expected heritability
p_lambda <- ggplot(df, aes(x = h2_als, y = h2_empirical, fill = K_als)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(aes(y = as.numeric(as.character(h2_als))),
             color = "black", shape = 18, size = 2,
             position = position_dodge(width = 0.8)) +
  labs(
    title = "Simulated vs Expected Heritability of ALS",
    subtitle = "Polygenic liability-threshold simulations across lambda values",
    x = "Input heritability (h²_ALS)",
    y = "Empirical heritability from simulations",
    fill = "Input prevalence (K_ALS)"
  ) +
  facet_wrap(~ lambda, scales = "free_x", labeller = labeller(lambda = function(x) paste0("Fertility rate: ", x))) +
  theme_minimal()
print(p_lambda)

# # Plot facetted by lambda for simulated vs expected prevalence
# p_gen <- ggplot(df, aes(x = K_als, y = K_empirical, fill = h2_als)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   geom_point(aes(y = as.numeric(as.character(K_als))),
#              color = "black", shape = 18, size = 2,
#              position = position_dodge(width = 0.8)) +
#   labs(
#     title = "Simulated vs Expected Prevalence of ALS",
#     subtitle = "Polygenic liability-threshold simulations across generations",
#     x = "Input K_ALS (prevalence)",
#     y = "Empirical prevalence from simulations",
#     fill = "Heritability (h²)"
#   ) +
#   facet_wrap(~ gen, labeller = label_both) +   # facet on lambda
#   theme_minimal()
# print(p_gen)

# # Plot facetted by gen (number of generations) for simulated vs expected heritability
# p_gen <- ggplot(df, aes(x = h2_als, y = h2_empirical, fill = K_als)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   geom_point(aes(y = as.numeric(as.character(h2_als))),
#              color = "black", shape = 18, size = 2,
#              position = position_dodge(width = 0.8)) +
#   labs(
#     title = "Simulated vs Expected Heritability of ALS",
#     subtitle = "Polygenic liability-threshold simulations across generations",
#     x = "Input heritability (h²_ALS)",
#     y = "Empirical heritability from simulations",
#     fill = "Input prevalence (K_ALS)"
#   ) +
#   facet_wrap(~ gen, labeller = label_both) +   # facet on number of generations
#   theme_minimal()
# print(p_gen)