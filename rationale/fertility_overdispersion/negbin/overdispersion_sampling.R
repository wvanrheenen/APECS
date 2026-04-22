#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

set.seed(123)

# Parameters
k <- 3
fert_rate <- data.frame(year = c(1940, 1970, 2000, 2030), mean_fertility = c(4, 2.5, 1.5, 1))
gen_years <- c(1940, 1970, 2000, 2030)
r_nb <- 2.5 
rho <- 0.15 
theta <- 1 / r_nb 
n_peds <- 1500

# Functions (unchanged)
init_ped <- function(k, r_nb, theta, gen_years){
  data.frame(gen = 0, id = "C0", pid = NA_character_, mid = NA_character_,
             sex = 1, yob = gen_years[1], fert_propensity = rgamma(1, shape=r_nb, scale=theta))
}

add_gen <- function(df_ped, k, r_nb, rho, theta, gen_years){
  g <- 0
  while(g < k){
    I1s <- filter(df_ped, gen == g)
    if(nrow(I1s) > 0) {
      for(i in 1:nrow(I1s)){
        parent_propensity <- I1s$fert_propensity[i]
        child_propensity <- rho * parent_propensity + (1-rho) * rgamma(1, shape=r_nb, scale=theta)
        
        I2 <- data.frame(gen = g, id = gsub("C", "P", I1s$id[i]), pid = NA_character_, 
                        mid = NA_character_, sex = abs(I1s$sex[i] - 1),
                        yob = I1s$yob[i], fert_propensity = parent_propensity)
        
        base_lambda <- fert_rate$mean_fertility[fert_rate$year == gen_years[g+1]]
        lambda_fert <- base_lambda * child_propensity
        n_II <- rnbinom(1, size=r_nb, mu=lambda_fert)
        
        # while(n_II == 0 & g == 0){
        #   n_II <- rnbinom(1, size=r_nb, mu=lambda_fert)
        # }
        
        if(n_II > 0){
          IIs <- data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
          colnames(IIs) <- colnames(df_ped)
          for(j in 1:n_II){
            IIs$gen[j] <- g + 1
            IIs$id[j] <- paste(I1s$id[i], j, sep="_")
            IIs$pid[j] <- ifelse(I1s$sex[i] == 0, I1s$id[i], I2$id)
            IIs$mid[j] <- ifelse(I1s$sex[i] == 1, I1s$id[i], I2$id)
            IIs$sex[j] <- sample(c(0,1), 1)
            IIs$yob[j] <- gen_years[g+2]
            IIs$fert_propensity[j] <- child_propensity
          }
          df_ped <- bind_rows(df_ped, I2, IIs)
        } else {
          df_ped <- bind_rows(df_ped, I2)
        }
      }
      g <- g + 1
    } else {
      g <- g - 1
    }
  }
  return(df_ped)
}

# SIMULATE PEDIGREES
all_peds <- data.frame()
for(i in 1:n_peds){
  ped <- init_ped(k, r_nb, theta, gen_years)
  ped <- add_gen(ped, k, r_nb, rho, theta, gen_years)
  ped$ped_id <- i
  all_peds <- bind_rows(all_peds, ped)
}

cat("Simulated", nrow(all_peds), "individuals\n\n")
write.table(all_peds, "all_pedigrees.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# FERTILITY ESTIMATES (one per mother)
fert_estimates <- data.frame(gen=numeric(), yob_children=numeric(), N_offspring=numeric(), ped_id=numeric())
max_gen <- max(all_peds$gen)

for(ped in unique(all_peds$ped_id)){
  ped_data <- filter(all_peds, ped_id == ped)
  for(i in 1:max_gen){
    women <- filter(ped_data, gen == i-1, sex == 1)
    for(j in 1:nrow(women)){
      II <- filter(ped_data, mid == women$id[j])
      if(nrow(II) > 0) {
        fert_estimates <- bind_rows(fert_estimates, data.frame(
          gen = i-1, yob_children = round(mean(II$yob)), 
          N_offspring = nrow(II), ped_id = ped))
      } else {
        fert_estimates <- bind_rows(fert_estimates, data.frame(
          gen = i-1, yob_children = women$yob[j] + 30, 
          N_offspring = 0, ped_id = ped))
      }
    }
  }
}
fwrite(fert_estimates, "./fertility_estimates_per_pedigree.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


# FAMILY PROPENSITY (one per family per gen)
family_propensity <- all_peds %>%
  group_by(ped_id, gen) %>%
  summarise(family_propensity = mean(fert_propensity), .groups = 'drop')

# FAMILY FERTILITY (mean per family per gen)
family_fert_mean <- fert_estimates %>%
  group_by(ped_id, gen) %>%
  summarise(mean_N_offspring = mean(N_offspring), .groups = 'drop')

# SUMMARY STATISTICS
fert_stats <- fert_estimates %>%
  group_by(gen) %>%
  summarise(
    n_individuals = n(),
    mean_fertility = round(mean(N_offspring), 3),
    .groups = 'drop') %>%
  mutate(target_fertility = c(4.0, 2.5, 1.5)[gen+1])  # G0=4.0, G1=2.5, G2=1.5

propensity_stats <- family_propensity %>%
  group_by(gen) %>%
  summarise(
    n_families = n(),
    mean_propensity = round(mean(family_propensity), 3),
    target_propensity = 1.0,
    .groups = 'drop')

# CORRELATIONS (G0-G1 family lineage)
family_lineage <- family_fert_mean %>%
  left_join(family_propensity, by = c("ped_id", "gen"))

# G0→G1
g0 <- family_lineage %>% filter(gen == 0)
g1 <- family_lineage %>% filter(gen == 1)
g0_g1 <- merge(g0, g1, by = "ped_id", suffixes = c("_g0", "_g1"))
cor_prop_g0g1 <- cor(g0_g1$family_propensity_g0, g0_g1$family_propensity_g1, use = "complete.obs")
cor_fert_g0g1 <- cor(g0_g1$mean_N_offspring_g0, g0_g1$mean_N_offspring_g1, use = "complete.obs")

# G1→G2  
g2 <- family_lineage %>% filter(gen == 2)
g1_g2 <- merge(g1, g2, by = "ped_id", suffixes = c("_g1", "_g2"))
cor_prop_g1g2 <- cor(g1_g2$family_propensity_g1, g1_g2$family_propensity_g2, use = "complete.obs")
cor_fert_g1g2 <- cor(g1_g2$mean_N_offspring_g1, g1_g2$mean_N_offspring_g2, use = "complete.obs")

# FINAL SUMMARY OUTPUT
cat("=== SIMULATION VALIDATION ===\n\n")

cat("INPUT PARAMETERS:\n")
cat(sprintf("  Heritability (ρ):     %.3f\n", rho))
cat(sprintf("  Propensity target:    %.3f\n", 1.0))
cat(sprintf("  Fertility targets:    G0=%.1f, G1=%.1f, G2=%1f\n", 4.0, 2.5, 1.5))
cat("\n")

cat("OUTPUT - MEAN FERTILITY (per mother):\n")
print(fert_stats, n = Inf)
cat("\n")

cat("OUTPUT - MEAN PROPENSITY (per family):\n")
print(propensity_stats, n = Inf)
cat("\n")

cat("OUTPUT - G0→G1 CORRELATIONS (family lineage):\n")
cat(sprintf("  Propensity: %.3f (target ρ=%.3f)\n", cor_prop_g0g1, rho))
cat(sprintf("  Fertility:  %.3f\n", cor_fert_g0g1))
cat(sprintf("  N families: %d\n", nrow(g0_g1)))

cat("OUTPUT - G1→G2 CORRELATIONS (family lineage):\n")
cat(sprintf("  Propensity: %.3f (target ρ=%.3f)\n", cor_prop_g1g2, rho))
cat(sprintf("  Fertility:  %.3f\n", cor_fert_g1g2))
cat(sprintf("  N families: %d\n", nrow(g1_g2)))


# ============================
# === PLOTTING SECTION (4 GENS)
# ============================

# --- 1. Family Propensity Distribution (All 4 Generations) ---
family_propensity <- all_peds %>%
  group_by(ped_id, gen) %>%
  summarise(family_propensity = mean(fert_propensity), .groups = 'drop')

# Theoretical Gamma density (mean=1.0)
x_prop <- seq(0, 3, length.out = 100)
gamma_theory <- data.frame(
  propensity = x_prop,
  density = dgamma(x_prop, shape = r_nb, scale = theta)
)

gen_stats <- family_propensity %>%
  group_by(gen) %>%
  summarise(mean_prop = mean(family_propensity), .groups = 'drop')

p_propensity <- ggplot(family_propensity, aes(x = family_propensity)) +
  geom_histogram(aes(y = after_stat(density), fill = factor(gen)), 
                 bins = 30, alpha = 0.7, color = "black", linewidth = 0.3) +
  geom_line(data = gamma_theory, aes(x = propensity, y = density),
            color = "black", linewidth = 0.5) +
  geom_vline(data = gen_stats, aes(xintercept = mean_prop),
             color = "black", linetype = "dotted", linewidth = 1, alpha = 0.8) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "black", linewidth = 1) +
  facet_wrap(~gen, scales = "free_y",
             labeller = labeller(gen = c("0" = "G0 (1940)", "1" = "G1 (1970)", 
                                        "2" = "G2 (2000)", "3" = "G3 (2030)"))) +
  scale_fill_manual(values = c("0" = "red", "1" = "steelblue", "2" = "forestgreen", "3" = "orange"),
                    name = "Generation", labels = c("G0", "G1", "G2", "G3")) +
  labs(title = "Family Propensity: Observed vs Theoretical Gamma (4 Generations)",
       x = "Family Fertility Propensity", y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black"))

print(p_propensity)
ggsave("family_propensity_4gen.png", p_propensity, width = 16, height = 4, dpi = 300)

# --- 2. Observed vs Theoretical Offspring (G0, G1, G2 mothers) ---
x_max <- 20
fert_targets <- c(4.0, 2.5, 1.5)  # G0, G1, G2 mothers

# Theoretical NB and Poisson for 3 mother generations
theoretical <- data.frame(
  offspring = rep(0:x_max, 3),
  gen = rep(0:2, each = x_max + 1),
  density = unlist(lapply(fert_targets, function(mu) dnbinom(0:x_max, mu = mu, size = r_nb))),
  type = "NegBin"
)

poisson <- data.frame(
  offspring = rep(0:x_max, 3),
  gen = rep(0:2, each = x_max + 1),
  density = unlist(lapply(fert_targets, function(mu) dpois(0:x_max, lambda = mu))),
  type = "Poisson"
)

all_theoretical <- bind_rows(theoretical, poisson)

p_combined <- ggplot(fert_estimates %>% filter(gen <= 2), aes(x = N_offspring)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, alpha = 0.6, 
                 fill = "steelblue", color = "black") +
  geom_line(data = all_theoretical, aes(x = offspring, y = density, color = type), 
            linewidth = 1.2) +
  facet_wrap(~gen, scales = "free_y",
             labeller = labeller(gen = c("0" = "G0 mothers (1940, target=4.0)",
                                        "1" = "G1 mothers (1970, target=2.5)",
                                        "2" = "G2 mothers (2000, target=1.5)"))) +
  labs(title = "Observed vs Theoretical: NegBin (r=2.5) vs Poisson (3 Mother Generations)",
       x = "Number of Offspring", y = "Density", color = "Distribution") +
  scale_color_manual(values = c("NegBin" = "red", "Poisson" = "darkgreen")) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_combined)
ggsave("observed_vs_theoretical_4gen.png", p_combined, width = 15, height = 5, dpi = 300)

cat("\nPlots saved:\n- family_propensity_4gen.png\n- observed_vs_theoretical_4gen.png\n")
