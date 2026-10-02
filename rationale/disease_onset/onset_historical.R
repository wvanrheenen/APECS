#' make sure warnings are treated as errors
options(warn = 1)
source("../../src/libraries_simPed.R")
library(RColorBrewer)

# Read age_of_onset data
disease_onset <- read.csv("../../data/age_of_onset/age_of_onset.csv", header = TRUE)

# Normalize each disease onset probability to sum to 1 independently
disease_onset_norm <- disease_onset %>%
  mutate(across(-age, ~ .x / sum(.x))) 

# Pivot to long format for ggplot
disease_long <- disease_onset_norm %>%
  pivot_longer(cols = -age, names_to = "disease", values_to = "probability") 

new_order <- c(
  "C9_ALS",
  "Polygenic_ALS",
  "C9_FTD",
  "Polygenic_FTD",
  "C9_Dementia",
  "Polygenic_Dementia"
)

new_labels <- gsub("^C9_", "Monogenic ", new_order)
new_labels <- gsub("^Polygenic_", "Polygenic ", new_labels)

disease_long$disease <- factor(disease_long$disease,
                               levels = new_order,
                               labels = new_labels)
                               
# Six-class ColorBrewer RdBu palette
rdbu_cols <- RColorBrewer::brewer.pal(8, "RdBu")

# Explicit disease mapping:
# ALS: darkest red / darkest blue
# FTD: medium red / medium blue
# Dementia: lightest red / lightest blue
disease_cols <- c(
  "Monogenic ALS"       = rdbu_cols[1],
  "Polygenic ALS"       = rdbu_cols[8],
  "Monogenic FTD"       = rdbu_cols[2],
  "Polygenic FTD"       = rdbu_cols[7],
  "Monogenic Dementia"  = rdbu_cols[3],
  "Polygenic Dementia"  = rdbu_cols[6]
)


p <- ggplot(disease_long, aes(x = age, y = probability, color = disease)) +
  geom_line(linewidth = 1) +
  labs(title = "(A) Historical Disease Onset",
       x = "Age at Onset (years)", y = "Density") +
  scale_color_manual(
    values = disease_cols,
    name = "Disease"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face="bold", margin = margin(b = 5)),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt"))

ggsave("historical_disease_onset_9x3.pdf", p, 
       width = 9, height = 3, units = "in", dpi = 300)

cat("Historical disease onset plot saved as PDF\n")
cat("RdBu-6 palette used for 6 disease-onset distributions\n")
