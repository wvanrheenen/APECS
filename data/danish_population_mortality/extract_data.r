library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# Read the deaths data, forcing all columns to character to avoid type issues
deaths_data <- read_excel(
  "danish_deaths_1980_2021.xlsx",
   col_types = "text"
)

# Clean age column
deaths_data <- deaths_data %>%
  mutate(Age = str_trim(Age))

# Pivot longer, excluding 'Subtotal'
deaths_long <- deaths_data %>%
  pivot_longer(
    cols = -Age,
    names_to = "year",
    values_to = "deaths"
  ) %>%
  filter(year != "Subtotal") %>%
  mutate(
    year = as.integer(year),
    deaths = as.numeric(deaths),
    age = ifelse(str_detect(Age, "99"), "99+", str_extract(Age, "\\d+"))
  ) %>%
  select(age, year, deaths)

# Repeat similar steps for the population file
population_data <- read_excel(
  "danish_population_1980_2021.xlsx",
  col_types = "text"
) %>%
  mutate(Age = str_trim(Age))

population_long <- population_data %>%
  pivot_longer(
    cols = -Age,
    names_to = "year",
    values_to = "population"
  ) %>%
  filter(year != "Subtotal") %>%
  mutate(
    year = as.integer(year),
    population = as.numeric(population),
    age = ifelse(str_detect(Age, "99"), "99+", str_extract(Age, "\\d+"))
  ) %>%
  select(age, year, population)

# Join deaths and population
mortality_data_yearly <- left_join(deaths_long, population_long, by = c("age", "year")) %>%
  mutate(mortality_rate = deaths / population)

als_cases <- data.frame(
  start_year = c(1980, 1984, 1989, 1994, 1999, 2004, 2009, 2014, 2018),
  end_year   = c(1983, 1988, 1993, 1998, 2003, 2008, 2013, 2017, 2021),
  cases      = c(270, 321, 395, 665, 708, 913, 992, 852, 827)
) %>%
  rowwise() %>%
  mutate(years = list(seq(start_year, end_year))) %>%
  unnest(years) %>%
  select(year = years, cases)  

mortality_weighted <- mortality_data_yearly %>%
  left_join(als_cases, by = "year") %>%
  filter(!is.na(cases))  # Keep only years with ALS data

weighted_mortality <- mortality_weighted %>%
  group_by(age) %>%
  summarise(
    weighted_mortality_rate = sum(mortality_rate * cases) / sum(cases)
  ) %>%
  ungroup() %>%
  mutate(
    # Convert "99+" to 99 for numeric sorting
    age_numeric = ifelse(age == "99+", 99, as.numeric(age))
  ) %>%
  arrange(age_numeric) %>%
  select(-age_numeric)  # Remove helper column if not needed

# Write to tab-separated txt file
write.table(weighted_mortality,
            file = "/hpc/hers_en/pbeele/simPed/data/danish_population_mortality/danish_population_deaths.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

cat("File 'danish_population_deaths.txt' created successfully with weighted mortality rate column.\n")
