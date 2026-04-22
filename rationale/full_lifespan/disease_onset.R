## disease_onset.R
## This script checks if simulated Mendelian/polygenic ALS/FTD/Dementia onset ages 
## match input parameter distributions from disease_onset.csv (density format)

library(dplyr)
library(ggplot2)
library(wesanderson)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
sim_file <- args[1]
onset_file <- args[2]
set_name <- args[3]
plot_dir <- args[4]

# Read in the simulation results and reference distribution  
simulations <- read.csv(sim_file, stringsAsFactors = FALSE)
disease_onset <- read.csv(onset_file, stringsAsFactors = FALSE)

# Create plots directory
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Darjeeling1 palette
darjeeling_cols <- wes_palette("Darjeeling1", 2)


## ===== MENDLIAN ALS vs C9_ALS =====
als_mendel_cases <- simulations %>%
  filter(mendel_ALS_Y_common == 1, !is.na(age_mendel_ALS_Y_common)) %>%
  pull(age_mendel_ALS_Y_common) %>%
  c(simulations %>%
      filter(mendel_ALS_Y_patho == 1, !is.na(age_mendel_ALS_Y_patho)) %>%
      pull(age_mendel_ALS_Y_patho))

cat("Mendelian ALS cases: n =", length(als_mendel_cases), "\n")
print(summary(als_mendel_cases))

c9_als_ref <- data.frame(
  age = disease_onset$age,
  density_c9_als = disease_onset$C9_ALS / 100
)

p1 <- ggplot() +
  stat_density(data = data.frame(age = als_mendel_cases), 
               aes(x = age, color = "Simulated"), 
               geom = "line", position = "identity") +
  geom_line(data = c9_als_ref, 
            aes(x = age, y = density_c9_als, color = "Historical") 
            ) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Source") +
  labs(title = "(B) ALS Monogenic Onset",
       x = "Age at Onset (years)", y = "Density") +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(file.path(plot_dir, "als_mendel_vs_c9als.pdf"), p1, 
       width = 3, height = 3, units = "in", dpi = 300)
       
## ===== POLYGENIC ALS vs Polygenic_ALS =====
als_poly_cases <- simulations %>%
  filter(polygenicY_ALS == 1, !is.na(age_polygenicY_ALS)) %>%
  pull(age_polygenicY_ALS)

cat("\nPolygenic ALS cases: n =", length(als_poly_cases), "\n")
print(summary(als_poly_cases))

poly_als_ref <- data.frame(
  age = disease_onset$age,
  density_poly_als = disease_onset$Polygenic_ALS / 100
)

p2 <- ggplot() +
  stat_density(data = data.frame(age = als_poly_cases), 
               aes(x = age, color = "Simulated"), 
               geom = "line", position = "identity") +               
  geom_line(data = poly_als_ref, 
            aes(x = age, y = density_poly_als, color = "Historical") 
            ) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Source") +
  labs(title = "(E) ALS Polygenic Onset",
       x = "Age at Onset (years)", y = "Density") +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(file.path(plot_dir, "als_polygenic_vs_polyals.pdf"), p2, 
       width = 3, height = 3, units = "in", dpi = 300)

## ===== MENDLIAN FTD vs C9_FTD =====
ftd_mendel_cases <- simulations %>%
  filter(mendel_FTD_Y_common == 1, !is.na(age_mendel_FTD_Y_common)) %>%
  pull(age_mendel_FTD_Y_common) %>%
  c(simulations %>%
      filter(mendel_FTD_Y_patho == 1, !is.na(age_mendel_FTD_Y_patho)) %>%
      pull(age_mendel_FTD_Y_patho),
    simulations %>%
      filter(mendel_FTD_Y_ftd == 1, !is.na(age_mendel_FTD_Y_ftd)) %>%
      pull(age_mendel_FTD_Y_ftd))

cat("\nMendelian FTD cases: n =", length(ftd_mendel_cases), "\n")
print(summary(ftd_mendel_cases))

c9_ftd_ref <- data.frame(
  age = disease_onset$age,
  density_c9_ftd = disease_onset$C9_FTD / 100
)

p3 <- ggplot() +
  stat_density(data = data.frame(age = ftd_mendel_cases), 
               aes(x = age, color = "Simulated"), 
               geom = "line", position = "identity") +
  geom_line(data = c9_ftd_ref, 
            aes(x = age, y = density_c9_ftd, color = "Historical") 
            ) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Source") +
  labs(title = "(C) FTD Monogenic Onset",
       x = "Age at Onset (years)", y = "Density") +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(file.path(plot_dir, "ftd_mendel_vs_c9ftd.pdf"), p3, 
       width = 3, height = 3, units = "in", dpi = 300)

## ===== POLYGENIC FTD vs Polygenic_FTD =====
ftd_poly_cases <- simulations %>%
  filter(polygenicY_FTD == 1, !is.na(age_polygenicY_FTD)) %>%
  pull(age_polygenicY_FTD)

cat("\nPolygenic FTD cases: n =", length(ftd_poly_cases), "\n")
print(summary(ftd_poly_cases))

poly_ftd_ref <- data.frame(
  age = disease_onset$age,
  density_poly_ftd = disease_onset$Polygenic_FTD / 100
)

p4 <- ggplot() +
  stat_density(data = data.frame(age = ftd_poly_cases), 
               aes(x = age, color = "Simulated"), 
               geom = "line", position = "identity") +   
  geom_line(data = poly_ftd_ref, 
            aes(x = age, y = density_poly_ftd, color = "Historical") 
            ) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Source") +
  labs(title = "(F) FTD Polygenic Onset",
       x = "Age at Onset (years)", y = "Density") +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(file.path(plot_dir, "ftd_polygenic_vs_polyftd.pdf"), p4, 
       width = 3, height = 3, units = "in", dpi = 300)

## ===== MENDLIAN DEMENTIA vs C9_Dementia =====
dem_mendel_cases <- simulations %>%
  filter(mendel_dem_Y_common == 1, !is.na(age_mendel_dem_Y_common)) %>%
  pull(age_mendel_dem_Y_common) %>%
  c(simulations %>%
      filter(mendel_dem_Y_patho == 1, !is.na(age_mendel_dem_Y_patho)) %>%
      pull(age_mendel_dem_Y_patho))

cat("\nMendelian Dementia cases: n =", length(dem_mendel_cases), "\n")
print(summary(dem_mendel_cases))

c9_dem_ref <- data.frame(
  age = disease_onset$age,
  density_c9_dem = disease_onset$C9_Dementia / 100
)

p5 <- ggplot() +
  stat_density(data = data.frame(age = dem_mendel_cases), 
               aes(x = age, color = "Simulated"), 
               geom = "line", position = "identity") +
  geom_line(data = c9_dem_ref, 
            aes(x = age, y = density_c9_dem, color = "Historical") 
            ) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Source") +
  labs(title = "(D) Dementia Monogenic Onset",
       x = "Age at Onset (years)", y = "Density") +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(file.path(plot_dir, "dem_mendel_vs_c9dem.pdf"), p5, 
       width = 3, height = 3, units = "in", dpi = 300)

## ===== POLYGENIC DEMENTIA vs Polygenic_Dementia =====
dem_poly_cases <- simulations %>%
  filter(polygenicY_dementia == 1, !is.na(age_polygenicY_dementia)) %>%
  pull(age_polygenicY_dementia)

cat("\nPolygenic Dementia cases: n =", length(dem_poly_cases), "\n")
print(summary(dem_poly_cases))

poly_dem_ref <- data.frame(
  age = disease_onset$age,
  density_poly_dem = disease_onset$Polygenic_Dementia / 100
)

p6 <- ggplot() +
  stat_density(data = data.frame(age = dem_poly_cases), 
               aes(x = age, color = "Simulated"), 
               geom = "line", position = "identity") +   
  geom_line(data = poly_dem_ref, 
            aes(x = age, y = density_poly_dem, color = "Historical") 
            ) +
  scale_color_manual(values = wes_palette(n = 2, name = "Darjeeling1"),
                     name = "Source") +
  labs(title = "(G) Dementia Polygenic Onset",
       x = "Age at Onset (years)", y = "Density") +
  coord_cartesian(xlim = c(20, 100)) +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin = margin(b = 2),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(5,5,5,5), "pt")
  )

ggsave(file.path(plot_dir, "dem_polygenic_vs_polydem.pdf"), p6, 
       width = 3, height = 3, units = "in", dpi = 300)

# Print comprehensive summary statistics table
cat("\n=== COMPREHENSIVE SUMMARY STATISTICS ===\n")
summary_stats <- data.frame(
  Metric = c("Mean", "Median", "SD", "Min", "Max", "N"),
  Mendelian_ALS = c(round(mean(als_mendel_cases), 2), round(median(als_mendel_cases), 2), 
                    round(sd(als_mendel_cases), 2), round(min(als_mendel_cases), 2), 
                    round(max(als_mendel_cases), 2), length(als_mendel_cases)),
  Polygenic_ALS = c(round(mean(als_poly_cases), 2), round(median(als_poly_cases), 2), 
                    round(sd(als_poly_cases), 2), round(min(als_poly_cases), 2), 
                    round(max(als_poly_cases), 2), length(als_poly_cases)),
  Mendelian_FTD = c(round(mean(ftd_mendel_cases), 2), round(median(ftd_mendel_cases), 2), 
                    round(sd(ftd_mendel_cases), 2), round(min(ftd_mendel_cases), 2), 
                    round(max(ftd_mendel_cases), 2), length(ftd_mendel_cases)),
  Polygenic_FTD = c(round(mean(ftd_poly_cases), 2), round(median(ftd_poly_cases), 2), 
                    round(sd(ftd_poly_cases), 2), round(min(ftd_poly_cases), 2), 
                    round(max(ftd_poly_cases), 2), length(ftd_poly_cases)),
  Mendelian_Dem = c(round(mean(dem_mendel_cases), 2), round(median(dem_mendel_cases), 2), 
                    round(sd(dem_mendel_cases), 2), round(min(dem_mendel_cases), 2), 
                    round(max(dem_mendel_cases), 2), length(dem_mendel_cases)),
  Polygenic_Dem = c(round(mean(dem_poly_cases), 2), round(median(dem_poly_cases), 2), 
                    round(sd(dem_poly_cases), 2), round(min(dem_poly_cases), 2), 
                    round(max(dem_poly_cases), 2), length(dem_poly_cases))
)
print(summary_stats)

cat("\nAll 6 plots saved to 'plots/' directory\n")
cat("✓ ALS Mendelian vs C9_ALS\n")
cat("✓ ALS Polygenic vs Polygenic_ALS\n") 
cat("✓ FTD Mendelian vs C9_FTD\n")
cat("✓ FTD Polygenic vs Polygenic_FTD\n")
cat("✓ Dementia Mendelian vs C9_Dementia\n")
cat("✓ Dementia Polygenic vs Polygenic_Dementia\n")
cat("Darjeeling1 palette: Simulated =", darjeeling_cols[1], ", Historical =", darjeeling_cols[2], "\n")
