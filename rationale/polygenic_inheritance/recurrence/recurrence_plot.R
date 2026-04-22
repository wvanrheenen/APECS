library(readr)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(rlang)

df <- read_tsv("polygenic_recurrence.out") %>%
  mutate(
    lambda = factor(lambda),
    K_als = factor(K_als),
    h2_als = as.numeric(h2_als),
    gen = factor(gen)
  )

# Create a summary data frame for dashed lines, similar to expected_segments in your monogenic plot
expected_segments <- df %>%
  arrange(lambda, h2_als) %>%
  group_by(lambda, K_als, h2_als) %>%
  summarize(
    y = unique(lambda1_theory),
    .groups = "drop"
  ) %>%
  mutate(
    x = as.numeric(factor(h2_als)) - 0.5,
    xend = as.numeric(factor(h2_als)) + 0.5,
    yend = y
  )

k_als_labeller <- function(x) paste0("K: ", x)

plot_lambda <- function(data, lambda_emp_col, lambda_theor_col, lambda_name) {
  ggplot(data, aes(x = factor(h2_als))) +
    geom_boxplot(aes(y = !!sym(lambda_emp_col), fill = lambda),
                 position = position_dodge(width = 0.8), alpha = 0.7) +
    geom_segment(data = expected_segments,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black",
                 linetype = "dashed",
                 inherit.aes = FALSE) +
    facet_wrap(~ K_als, scales = "free_x", labeller = labeller(K_als = k_als_labeller)) +
    scale_fill_manual(values = wes_palette("Darjeeling1", n = length(unique(data$lambda))),
                      name = "Fertility rate") +
    labs(
      title = paste("Simulated vs Theoretical Recurrence Risk:", lambda_name),
      subtitle = "Theoretical risk shown as dashed line",
      x = "Heritability (h2)",
      y = "Recurrence Risk"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(size = 14, face = "bold")
    )
}

png("polygenic_recurrence.png", width = 10, height = 4, units = "in", res = 300)
p_lambda1 <- plot_lambda(df, 
                         lambda_emp_col = "lambda1_empirical", 
                         lambda_theor_col = "lambda1_theory",
                         lambda_name = "First-degree relatives")
print(p_lambda1)
dev.off()
