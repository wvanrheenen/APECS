# This is a script to check if the lifetime risk, as derived from Levison in theory, matches the lifetime risk on developing ALS in our simulations
# In order to do so, we simulate 1000s of individuals, not pedigrees
# For these people, we know their life-expectancy, based on their year of birth, which includes death due to ALS 
# We give these people an age-specific risk to develop ALS
# We assume they all die 3 years after their ALS diagnosis (which isn't true, but serves the sanity check) 
# Since we know year of birth, year of diagnosis, and year of death, we can now calculate (cumulative) incidence rates
# We can also calculate age-specific mortality rates because we know how old someone was at time of death and how many people lived passed that age

# Load libraries and functions
source("../../../src/libraries_simPed.R")
source("../../../src/functions_simPed.R")

## Create individual
#' ### Function to initiate pedigree with one founder with a mutation
init_ped = function(life_expectancy, current_year, sim_nr){
  core_ped = data.frame(gen = 0, # generation
                        id = sim_nr, # individual id - founder
                        pid = as.character(NA), # parental id
                        mid = as.character(NA), # maternal id
                        sex = 1, # sex, 0 = male, 1 = female
                        a1_common = 0,
                        a2_common = 0,
                        a1_patho = 0,
                        a2_patho = 0) # first allele at disease locus
  # simulate life expectancy
  core_ped$yob = sample(1925:1975, 1)
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
        TRUE ~ NA_real_
      )
    )
}



## Input parameters 
life_expectancy = read.table(file.path("../../../data/life_expectancy_at_fifteen/OWID_life_expectancy_Netherlands_1835_2085.txt"), header=T)
current_year = 2025
penetrance_ALS_common = read.table(file.path("../../../data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
penetrance_ALS_patho = read.table(file.path("../../../data/lifetime_disease_risk/lifetime_ALS_c9_risk.txt" ), header=T) # derived from Van Wijk 2024
h2_ALS = 0.5 # additive polygenic heritability
K_ALS = read.table(file.path("../../../data/lifetime_disease_risk/lifetime_ALS_risk.txt" ), header=T) # derived from Levison 2025
rg = 0.8 

## Interactive use of functions:
individual = init_ped(life_expectancy=life_expectancy, current_year=current_year, sim_nr=1)
individual = add_pheno(individual, penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, rg)
individual = years(individual) 
print(individual)

## PART 2: Calculate population statistics
# Helper function to expand individual timelines
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

# restate current_year
current_year <- 2085

# Create two individuals with unique IDs
individuals = map_df(1:100, function(i) {
  init_ped(life_expectancy, current_year, sim_nr = paste0("sim", i)) %>%
    add_pheno(penetrance_ALS_common, penetrance_ALS_patho, h2_ALS, K_ALS, rg) %>%
    years()
})

# Expand timelines for each individual and combine
expanded_years = individuals %>%
  group_split(id) %>%
  map_df(~expand_individual_years(.x, current_year))

# Aggregate counts: number of people alive per age-year
age_year_population_ordered <- expanded_years %>%
  group_by(age, year) %>%
  summarise(population = sum(alive), .groups = "drop") %>%
  pivot_wider(
    names_from = year,
    values_from = population,
    values_fill = 0
  ) %>%
  reorder_year_columns()

## Functions to create deaths in population 
age_year_death_ordered <- expanded_years %>%
  group_by(age, year) %>%
  summarise(deaths = sum(death), .groups = "drop") %>%
  pivot_wider(
    names_from = year,
    values_from = deaths,
    values_fill = 0
  ) %>%
  reorder_year_columns()

# Aggregate ALS onset counts by age-year
age_year_als_onset_ordered <- expanded_years %>%
  group_by(age, year) %>%
  summarise(als_onsets = sum(als_onset), .groups = "drop") %>%
  pivot_wider(
    names_from = year,
    values_from = als_onsets,
    values_fill = 0
  ) %>%
  reorder_year_columns()

## STEP 3: calculate population metrics
# Mortality rate (= age-specific deaths in a year / total population for that age in a year)

age_specific_mortality_rate <- function(year_of_interest, deaths_df, pop_df) {
  # Ensure year is character (matches column names)
  year_col <- as.character(year_of_interest)
  
  # Check if the year exists in both dataframes
  if (!(year_col %in% colnames(deaths_df))) stop("Year not found in deaths dataframe")
  if (!(year_col %in% colnames(pop_df))) stop("Year not found in population dataframe")
  
  tibble::tibble(
    age = deaths_df$age,
    deaths = deaths_df[[year_col]],
    population = pop_df[[year_col]],
    mortality_rate = ifelse(pop_df[[year_col]] > 0, deaths_df[[year_col]] / pop_df[[year_col]], NA_real_)
  )
}

als_incidence_rate_period <- function(start_year, end_year, als_onset_df, pop_df) {
  # Ensure years are character to match column names
  years <- as.character(seq(as.numeric(start_year), as.numeric(end_year)))
  
  # Check if all years exist in both dataframes
  missing_als <- setdiff(years, colnames(als_onset_df))
  missing_pop <- setdiff(years, colnames(pop_df))
  if (length(missing_als) > 0) stop(paste("Years not found in ALS onset dataframe:", paste(missing_als, collapse = ", ")))
  if (length(missing_pop) > 0) stop(paste("Years not found in population dataframe:", paste(missing_pop, collapse = ", ")))
  
  # Calculate total cases and person-years for each age
  total_cases <- als_onset_df %>%
    dplyr::select(age, dplyr::all_of(years)) %>%
    dplyr::mutate(total_cases = rowSums(dplyr::across(-age)))
  
  total_personyears <- pop_df %>%
    dplyr::select(age, dplyr::all_of(years)) %>%
    dplyr::mutate(total_personyears = rowSums(dplyr::across(-age)))
  
  # Join and calculate incidence rate
  result <- total_cases %>%
    dplyr::select(age, total_cases) %>%
    dplyr::left_join(
      total_personyears %>% dplyr::select(age, total_personyears),
      by = "age"
    ) %>%
    dplyr::mutate(
      incidence_rate_per_100k = ifelse(
        total_personyears > 0,
        (total_cases / total_personyears) * 100000,
        NA_real_
      )
    )
  
  return(result)
}

## STEP 4: calculate lifetime risk following the current probability method
calculate_lifetime_risk <- function(als_onset_df, death_df, pop_df, reference_year = NULL) {
  if (is.null(reference_year)) {
    reference_year <- tail(colnames(pop_df), 1)
  }
  age <- pop_df$age
  pop <- pop_df[[reference_year]]
  als_onsets <- als_onset_df[[reference_year]]
  deaths <- death_df[[reference_year]]
  
  als_incidence_hazard <- ifelse(pop > 0, als_onsets / pop, 0)
  mortality_hazard <- ifelse(pop > 0, deaths / pop, 0)
  
  survival_to_age <- cumprod(1 - dplyr::lag(mortality_hazard, default = 0))
  prob_als_at_age <- survival_to_age * als_incidence_hazard
  
  tibble::tibble(
    age = age,
    pop = pop,
    als_onsets = als_onsets,
    deaths = deaths,
    als_incidence_hazard = als_incidence_hazard,
    mortality_hazard = mortality_hazard,
    survival_to_age = survival_to_age,
    prob_als_at_age = prob_als_at_age
  ) -> risk_table
  
  lifetime_risk <- sum(risk_table$prob_als_at_age, na.rm = TRUE)
  cat(sprintf("Estimated lifetime risk of ALS: %.2f%%\n", 100 * lifetime_risk))
  return(list(lifetime_risk = lifetime_risk, risk_table = risk_table))
}