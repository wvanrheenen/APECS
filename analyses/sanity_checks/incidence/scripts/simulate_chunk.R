# scripts/simulate_chunk.R

args <- commandArgs(trailingOnly = TRUE)
chunk_id <- as.integer(args[1])
n_indiv <- as.integer(args[2])
current_year <- as.integer(args[3])
outfile <- args[4]


suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(tidyverse))
suppressMessages(library(MASS))
suppressMessages(library(mvnfast))
suppressMessages(library(kinship2))
suppressMessages(library(pedtools))
suppressMessages(library(ribd))
suppressMessages(library(RColorBrewer))
suppressMessages(library(this.path))
suppressMessages(library(wesanderson))
suppressMessages(library(broom))
suppressMessages(library(fmsb))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


# Load your functions and data
init_ped = function(life_expectancy, current_year, sim_nr){
  core_ped = data.frame(gen = 0, # generation
                        id = sim_nr, # individual id - founder
                        pid = as.character(NA), # parental id
                        mid = as.character(NA), # maternal id
                        sex = 1, # sex, 0 = male, 1 = female
                        a1_common = 1,
                        a2_common = 0,
                        a1_patho = 0,
                        a2_patho = 0) # first allele at disease locus
  # simulate life expectancy
  core_ped$yob = sample(1900:2025, 1)
  core_ped$life_expectancy = round(rgamma(1, shape=150, scale=(life_expectancy$life_expectancy[life_expectancy$year_of_birth == (core_ped$yob)])/150))
  # simulate age
  core_ped$age = current_year - core_ped$yob
  if(core_ped$age > core_ped$life_expectancy){
    core_ped$status = "dead"
    core_ped$age_censored = core_ped$life_expectancy
  } else {
    core_ped$status = "alive"
    core_ped$age_censored = core_ped$age
  }
  return(core_ped)
}

## Create phenotype
# Calculate age-specific incidence as difference of cumulative incidence

#' ### Function to add genetic values for Mendelian and polygenic inheritance
add_pheno = function(df_ped, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, rg) {

  # Mendelian ALS (autosomal dominant)
  df_ped$age_penetrance_ALS_common = sapply(df_ped$age_censored, function(age) {
    penetrance_ALS_common$cum_incidence_c9_ALS[which.min(abs(penetrance_ALS_common$age - age))]
  })
  df_ped$age_penetrance_ALS_patho = sapply(df_ped$age_censored, function(age) {
    penetrance_ALS_patho$cum_incidence_patho_ALS[which.min(abs(penetrance_ALS_patho$age - age))]
  })

  df_ped = df_ped %>%
  mutate(
    mendel_ALS_y1_common = ifelse(a1_common == 1, rbinom(n(), 1, age_penetrance_ALS_common), 0),
    mendel_ALS_y2_common = ifelse(a2_common == 1, rbinom(n(), 1, age_penetrance_ALS_common), 0),
    mendel_ALS_y1_patho = ifelse(a1_patho == 1, rbinom(n(), 1, age_penetrance_ALS_patho), 0),
    mendel_ALS_y2_patho = ifelse(a2_patho == 1, rbinom(n(), 1, age_penetrance_ALS_patho), 0), 
    mendel_ALS_Y_common = ifelse(mendel_ALS_y1_common + mendel_ALS_y2_common > 0, 1, 0),
    mendel_ALS_Y_patho = ifelse(mendel_ALS_y1_patho + mendel_ALS_y2_patho > 0, 1, 0),   
    mendel_ALS_Y = ifelse(mendel_ALS_y1_common + mendel_ALS_y2_common + mendel_ALS_y1_patho + mendel_ALS_y2_patho > 0, 1, 0)
  )

  # polygenic model:
  df_ped$G_ALS = NA
  
  # Sample G from multivariate normal for founders
  founders = which(is.na(df_ped$pid)) 
  df_ped$G_ALS[founders] = rnorm(length(founders), mean=0, sd=sqrt(h2_ALS)) 

  # sample non-genetic value E:
  # Without genetic correlation
  df_ped$E_ALS = rnorm(nrow(df_ped), 0, sqrt(1-h2_ALS))

  # define phenotype ALS
  df_ped$P_ALS = df_ped$G_ALS + df_ped$E_ALS
  df_ped$age_related_K_ALS = sapply(df_ped$age_censored, function(age) {
  K_ALS$cum_risk_adjusted[which.min(abs(K_ALS$age - age))]
  })
  df_ped$LT_ALS = -qnorm(df_ped$age_related_K_ALS, 0, 1) # liability threshold
  df_ped$polygenicY_ALS = ifelse(df_ped$P_ALS > df_ped$LT_ALS, 1, 0)
  
  # final phenotype
  df_ped$Y_ALS = ifelse(df_ped$polygenicY_ALS + df_ped$mendel_ALS_Y > 0, 1, 0)
  return(df_ped)
}

## Extract years of interest
years = function(df) {
  df %>%
    mutate(
      yod = ifelse(status == "dead", yob + age_censored, NA_real_),
      year_censored = yob + age_censored,
      onset_ALS = case_when(
        status == "dead" & Y_ALS == 1 ~ yod - 3,
        status == "alive" & Y_ALS == 1 ~ year_censored - 2,
        TRUE ~ NA_real_
      )
    )
}

expand_individual_years = function(indiv, current_year) {
  yob = indiv$yob
  max_year = ifelse(indiv$status == "dead", indiv$yod, current_year)
  tibble(
    year = yob:max_year,
    age = year - yob,
    alive = 1,
    death = ifelse(year == indiv$yod, 1, 0),
    als_onset = ifelse(year == indiv$onset_ALS, 1, 0)
  ) %>%
    replace_na(list(death = 0, als_onset = 0))
}

create_age_year_population <- function(indiv, current_year) {
  # Expand individual's timeline
  df = expand_individual_years(indiv, current_year)
  # Create population count matrix
  df %>%
    group_by(age, year) %>%
    summarise(population = sum(alive), .groups = "drop") %>%
    pivot_wider(
      names_from = year,
      values_from = population,
      values_fill = 0
    ) %>%
    complete(age = full_seq(age, 1), fill = list(population = 0)) %>%
    arrange(age)
}

reorder_year_columns <- function(df) {
  # Extract all columns except 'age'
  year_cols <- setdiff(names(df), "age")
  # Convert year column names to numeric for proper ordering
  ordered_cols <- year_cols[order(as.numeric(year_cols))]
  # Select 'age' first, then ordered year columns
  dplyr::select(df, age, dplyr::all_of(ordered_cols))
}

current_year = current_year
life_expectancy = read.table(file.path("/hpc/hers_en/pbeele/simPed/data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
penetrance_ALS_common = read.table(file.path("/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = read.table(file.path("/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_ALS = 0.5 # additive polygenic heritability
K_ALS = read.table(file.path("/hpc/hers_en/pbeele/simPed/data/lifetime_disease_risk/lifetime_ALS_risk.txt" ), header=T) # derived from Levison 2025
rg = 0.8 

# Simulate n_indiv individuals
individuals <- purrr::map_df(1:n_indiv, function(i) {
  sim_nr <- paste0("chunk", chunk_id, "_sim", i)
  init_ped(life_expectancy, current_year, sim_nr) %>%
    add_pheno(penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, rg) %>%
    years()
})

expanded_years <- individuals %>%
  dplyr::group_split(id) %>%
  purrr::map_df(~expand_individual_years(.x, current_year))

saveRDS(expanded_years, outfile)