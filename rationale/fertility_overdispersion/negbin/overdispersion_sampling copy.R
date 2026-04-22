library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
set.seed(123)

## We apply 'propensity' per pedigree: the trend to get the same amount of offspring within a family, mirrorring
## Real live cultural clusters that produce more/less offspring than people in their time

# Parameters from study
k <- 2
fert_rate <- data.frame(year = c(1940, 1970, 2000), mean_fertility = c(3, 1.5, 1))  # Gen-specific
gen_years <- c(1940, 1970, 2000)  # G0=1920 (high fert), G1=1950, G2=1980 (low fert)
r_nb <- 2.5 
rho <- 0.15 
theta <- 1 / r_nb 
n_peds <- 5000

#' ### Function to initiate pedigree with fertility propensity
init_ped = function(k, r_nb, theta, gen_years){  # Add gen_years parameter
  core_ped = data.frame(gen = 0,
                        id = "C0",
                        pid = as.character(NA),
                        mid = as.character(NA),
                        sex = 1,
                        yob = gen_years[1],  # 1940 for G0
                        fert_propensity = rgamma(1, shape=r_nb, scale=theta))
  return(core_ped)
}

#' ### Function to add next generation with CORRELATED fertility
add_gen = function(df_ped, lambda, k, r_nb, rho, theta, gen_years){
  g = 0
  while(g < k){
    I1s = filter(df_ped, gen == g)
    if(nrow(I1s) > 0) {
      for(i in 1:nrow(I1s)){
        # NEW: Child inherits parent's fertility propensity
        parent_propensity <- I1s$fert_propensity[i]
        child_propensity <- rho * parent_propensity + (1-rho) * rgamma(1, shape=r_nb, scale=theta)
        
        I2 = data.frame(gen = g, 
                        id = gsub("C", "P", I1s$id[i]),
                        pid = NA, mid = NA,
                        sex = abs(I1s$sex[i] - 1),
                        yob = I1s$yob[i],  # ← ADD THIS
                        fert_propensity = parent_propensity)  # Partner shares family trait
        
        # NEW: NegBin with family-specific rate
        # KEY FIX: Time-specific base rate by generation
        base_lambda <- fert_rate$mean_fertility[fert_rate$year == gen_years[g+1]]  # Child's era
        lambda_fert <- base_lambda * child_propensity  # Scale by family trait
        
        n_II = rnbinom(1, size=r_nb, mu=lambda_fert)
        
        # Ensure gen 0 has kids
        while(n_II == 0 & g == 0){
          n_II = rnbinom(1, size=r_nb, mu=lambda_fert)
        }
        
        IIs = data.frame(matrix(NA, nrow=n_II, ncol=ncol(df_ped)))
        colnames(IIs) = colnames(df_ped)
        j = 0
        while(j < n_II){
          j = j+1
          IIs$gen[j] = g + 1
          IIs$id[j] = paste(I1s$id[i], j, sep="_")
          IIs$pid[j] = ifelse(I1s$sex[i] == 0, I1s$id[i], I2$id)
          IIs$mid[j] = ifelse(I1s$sex[i] == 1, I1s$id[i], I2$id)
          IIs$sex[j] = sample(c(0,1), 1)
          IIs$yob[j] = gen_years[g+2]  # G1 kids=1970 (index 2), G2 kids=2000 (index 3)
          IIs$fert_propensity[j] = child_propensity  # ALL siblings inherit same family trait
        }
        
        if(n_II > 0){
          df_ped = bind_rows(df_ped, I2, IIs)
        } else {
          df_ped = bind_rows(df_ped, I2)
        }
      }
      g = g+1
    } else {
      g = g-1
    }
  }
  return(df_ped)
}

# SIMULATE 5000 three-generation pedigrees
all_peds <- data.frame()
for(i in 1:n_peds){
  ped <- init_ped(k, r_nb, theta, gen_years)
  ped <- add_gen(ped, lambda, k, r_nb, rho, theta, gen_years)
  ped$ped_id <- i
  all_peds <- bind_rows(all_peds, ped)
}
print(paste("Simulated", nrow(all_peds), "individuals"))

# Save all_peds to TSV file (tab-separated)
write.table(all_peds, "all_pedigrees.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## Plot distribution of offspring per individual
# Estimate fertility per individual (your exact method)
fert_estimates <- data.frame(gen=numeric(), yob_children=numeric(), N_offspring=numeric(), ped_id=numeric())

max_gen <- max(all_peds$gen)
for(ped in unique(all_peds$ped_id)){
  ped_data <- filter(all_peds, ped_id == ped)
  
  for(i in 1:max_gen){
    women <- filter(ped_data, gen == i-1 & sex == 1)  # Mothers only
    for(j in 1:nrow(women)){
      II <- filter(ped_data, mid == women$id[j])  # Her kids
      if(nrow(II) > 0) {
        fert_estimates_i <- data.frame(
          gen = i-1,
          yob_children = round(mean(II$yob)),
          N_offspring = nrow(II),
          ped_id = ped
        )
      } else {
        fert_estimates_i <- data.frame(
          gen = i-1,
          yob_children = women$yob[j] + 30,  # mean_gen_yr ~30
          N_offspring = 0,
          ped_id = ped
        )
      }
      fert_estimates <- bind_rows(fert_estimates, fert_estimates_i)
    }
  }
}

fwrite(fert_estimates, "./fertility_estimates_per_pedigree.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## Calculate statistics on fert rate
fert_stats <- fert_estimates %>%
  group_by(gen) %>%
  summarise(
    individuals = n(),
    mean_fertility = round(mean(N_offspring), 2),
    var_fertility = round(var(N_offspring), 2),
    sd_fertility = round(sd(N_offspring), 2),
    overdispersion = round(var_fertility / mean_fertility, 2),
    prop_childless = round(mean(N_offspring == 0) * 100, 1),
    .groups = 'drop'
  ) %>%
  mutate(
    target_mean = case_when(gen == 0 ~ 3.0, gen == 1 ~ 1.5),
    mean_error = round((mean_fertility - target_mean) / target_mean * 100, 1)
  )

print(fert_stats)

## theoretical distribution of offspring per generation (depending on lambda)
# Plot observed histograms + theoretical NB curves for each generation
x_max <- 20  # Max offspring to plot

# Create theoretical NB densities for each generation
theoretical_g01 <- data.frame(
  offspring = rep(0:x_max, 2),  # Only 2 generations
  gen = rep(c(0,1), each = x_max+1),
  density = c(
    dnbinom(0:x_max, mu = 3.0, size = r_nb),    # Gen 0 only
    dnbinom(0:x_max, mu = 1.5, size = r_nb)     # Gen 1 only
  )
)

# Add Poisson theoretical curves for comparison
poisson_g01 <- data.frame(
  offspring = rep(0:x_max, 2),
  gen = rep(c(0,1), each = x_max+1),
  density = c(
    dpois(0:x_max, lambda = 3.0),    # Poisson Gen 0
    dpois(0:x_max, lambda = 1.5)     # Poisson Gen 1
  ),
  type = "Poisson"
)

theoretical_g01$type <- "NegBin"  # Label existing curves

# Combine both distributions
all_theoretical <- bind_rows(poisson_g01, theoretical_g01)

# Enhanced plot with Poisson + NegBin
p_combined <- ggplot(fert_estimates, aes(x = N_offspring)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, alpha = 0.6, 
                 fill = "steelblue", color = "black") +
  geom_line(data = all_theoretical, aes(x = offspring, y = density, color = type), 
            linewidth = 1.2) +
  facet_wrap(~gen, labeller = labeller(gen = c("0" = "Gen 0 (1940, mean fert=3.0)", 
                                               "1" = "Gen 1 (1970, mean fert=1.5)")), 
             scales = "free_y") +
  labs(title = "Observed vs Theoretical: NegBin (r=2.5) vs Poisson", 
       x = "Number of Offspring", y = "Density",
       color = "Distribution") +
  scale_color_manual(values = c("NegBin" = "red", "Poisson" = "darkgreen")) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_combined)
ggsave("observed_vs_poisson_negbin_g01.png", p_combined, width = 12, height = 5, dpi = 300)

### Calculate mean propensity and distribution
# ONE propensity per pedigree per generation (family average)
propensity_stats <- all_peds %>%
  group_by(ped_id, gen) %>%
  summarise(
    family_propensity = mean(fert_propensity),  # Average within pedigree
    n_individuals = n(),
    .groups = 'drop'
  ) %>%
  group_by(gen) %>%
  summarise(
    n_families = n(),
    mean_propensity = round(mean(family_propensity), 3),
    sd_propensity = round(sd(family_propensity), 3),
    .groups = 'drop'
  )
print(propensity_stats)

family_propensity <- all_peds %>%
  group_by(ped_id, gen) %>%
  summarise(family_propensity = mean(fert_propensity), .groups = 'drop')

# Theoretical Gamma density (shape=2.5, scale=0.4, mean=1.0)
x_prop <- seq(0, 3, length.out = 100)
gamma_theory <- data.frame(
  propensity = x_prop,
  density = dgamma(x_prop, shape = r_nb, scale = theta),
  type = "Theoretical Gamma"
)

gen_stats <- family_propensity %>%
  group_by(gen) %>%
  summarise(mean_prop = round(mean(family_propensity), 3), .groups = 'drop')

# Enhanced plot
p_propensity <- ggplot(family_propensity, aes(x = family_propensity)) +
  # Colored histograms by generation
  geom_histogram(aes(y = after_stat(density), fill = factor(gen)), 
                 bins = 30, alpha = 0.7, color = "black", linewidth = 0.3) +
  # Theoretical Gamma curve
  geom_line(data = gamma_theory, aes(x = propensity, y = density), 
            color = "black", linewidth = 0.5) +
  # Mean propensity annotations
  geom_text(data = gen_stats, aes(x = Inf, y = Inf, label = paste0("Mean: ", mean_prop)),
            hjust = 1.1, vjust = 1.8, size = 4, fontface = "bold") +
  # Mean lines
  geom_vline(data = gen_stats, aes(xintercept = mean_prop), 
             color = "black", linetype = "dotted", linewidth = 1, alpha = 0.8) +
  # Population mean
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "black", linewidth = 1) +
  # Facet with custom labels + colors
  facet_wrap(~gen, scales = "free_y", 
             labeller = labeller(gen = c("0" = "Gen 0", "1" = "Gen 1", "2" = "Gen 2"))) +
  # Colors: Red=Gen0, Blue=Gen1, Green=Gen2
  scale_fill_manual(values = c("0" = "red", "1" = "steelblue", "2" = "forestgreen"),
                    name = "Generation", labels = c("Gen 0", "Gen 1", "Gen 2")) +
  labs(title = "Family Propensity: Observed vs Theoretical Gamma", 
       x = "Family Fertility Propensity", y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black"))

print(p_propensity)
ggsave("family_propensity_vs_gamma.png", p_propensity, width = 12, height = 4, dpi = 300)


## Check correlation of parent/offspring fertility rate
family_fert_mean <- fert_estimates %>%
  group_by(ped_id, gen) %>%
  summarise(
    mean_N_offspring = mean(N_offspring),  # Mean fertility of all mothers in family
    n_mothers = n(),
    .groups = 'drop'
  )

# 2. Merge with family propensity
family_lineage <- family_fert_mean %>%
  left_join(family_propensity, by = c("ped_id", "gen"))

# 3. G0→G1 family-lineage correlation (grandmother → mother generation)
g0_lineage <- family_lineage %>% filter(gen == 0)
g1_lineage <- family_lineage %>% filter(gen == 1)
g0_g1_lineage <- merge(g0_lineage, g1_lineage, by = "ped_id", all = TRUE,
                       suffixes = c("_g0", "_g1"))

# 4. Correlations down the family line
cor_prop_lineage <- cor(g0_g1_lineage$family_propensity_g0, 
                       g0_g1_lineage$family_propensity_g1, use = "complete.obs")
cor_fert_lineage <- cor(g0_g1_lineage$mean_N_offspring_g0, 
                       g0_g1_lineage$mean_N_offspring_g1, use = "complete.obs")

cat(sprintf("✅ G0→G1 FAMILY LINEAGE correlations:\n"))
cat(sprintf("   Propensity: %.3f (target ρ=0.15)\n", cor_prop_lineage))
cat(sprintf("   Mean Fertility: %.3f\n", cor_fert_lineage))
cat(sprintf("   N family lines: %d\n", nrow(g0_g1_lineage)))

