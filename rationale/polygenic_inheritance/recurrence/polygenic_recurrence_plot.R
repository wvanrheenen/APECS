### Polygenic inheritance rationale plotting script ###

# Load libraries
options(warn = 2)
library(readr)
library(dplyr)
library(ggplot2)
library(wesanderson)

# Read the polygenic output
df <- read_tsv("polygenic_recurrence.out")

# Make sure appropriate columns are factor or numeric
df <- df %>%
  mutate(
    lambda = factor(lambda),
    K_als = factor(K_als),       # disease prevalence for faceting
    h2_als = as.numeric(h2_als), # heritability on x-axis
    gen = factor(gen)            # if needed for coloring, else ignore
  )

# A helper function to plot for a given lambda type
plot_lambda <- function(data, lambda_emp_col, lambda_theor_col, lambda_name) {
  ggplot(data, aes(x = factor(h2_als))) +
    geom_boxplot(aes(y = !!rlang::sym(lambda_emp_col), fill = lambda),
                 position = position_dodge(width = 0.8), alpha = 0.6) +
    geom_point(aes(y = !!rlang::sym(lambda_theor_col), group = lambda),
               color = "black", shape = 18, size = 2,
               position = position_dodge(width = 0.8)) +
    facet_wrap(~ K_als, scales = "free_x", labeller = labeller(K_als = function(x) paste0("K: ", x))) +
    scale_fill_manual(values = wes_palette(n=length(unique(data$lambda)), name="Darjeeling1")) +
    labs(
      title = paste("Simulated vs Theoretical Recurrence Risk:", lambda_name),
      subtitle = "Theoretical risk marked by black diamond",
      x = "Heritability (h2)",
      y = paste("Theoretical and Simulated Recurrence Risk"),
      fill = "Fertility rate"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Plot lambda1
pdf("polygenic_recurrence_risk.pdf", width = 7, height = 5)
p_lambda1 <- plot_lambda(df, 
                        lambda_emp_col = "lambda1_empirical", 
                        lambda_theor_col = "lambda1_theory",
                        lambda_name = "First-degree relatives")
print(p_lambda1)                        
dev.off()

# # Plot lambda2
# p_lambda2 <- plot_lambda(df, 
#                         lambda_emp_col = "lambda2_empirical", 
#                         lambda_theor_col = "lambda2_theory",
#                         lambda_name = "Second-degree relatives")
# print(p_lambda2)


## Correlation
cor_lambda1_spearman <- cor(df$lambda1_empirical, df$lambda1_theory, method = "spearman", use = "complete.obs")
cor_lambda2_spearman <- cor(df$lambda2_empirical, df$lambda2_theory, method = "spearman", use = "complete.obs")

cor_lambda1_spearman
cor_lambda2_spearman