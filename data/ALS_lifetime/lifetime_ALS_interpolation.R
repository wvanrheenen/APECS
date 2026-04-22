library(dplyr)
library(ggplot2)
library(scales)


# 1. Load the raw digitized data
raw_data <- read.csv("/hpc/hers_en/pbeele/simPed/data/ALS_lifetime/lifetime_ALS_risk.csv",
                     header = FALSE, col.names = c("age_raw", "incidence_raw"))

# 2. For each raw age, find floor and ceiling integer ages
raw_data <- raw_data %>%
  mutate(
    floor_age = floor(age_raw),
    ceil_age = ceiling(age_raw),
    weight_floor = 1 - (age_raw - floor_age),
    weight_ceil = 1 - (ceil_age - age_raw)
  )

# 3. Create weighted contributions for floor and ceiling ages
floor_data <- raw_data %>%
  select(age = floor_age, incidence_raw, weight = weight_floor)

ceil_data <- raw_data %>%
  select(age = ceil_age, incidence_raw, weight = weight_ceil)

# 4. Combine and calculate weighted incidence sums and weights by age
combined_data <- bind_rows(floor_data, ceil_data) %>%
  group_by(age) %>%
  summarise(
    weighted_sum = sum(incidence_raw * weight, na.rm = TRUE),
    total_weight = sum(weight, na.rm = TRUE)
  ) %>%
  mutate(
    incidence_per_100k_PY = weighted_sum / total_weight
  ) %>%
  ungroup()

# 5. Create full age sequence and interpolate missing ages as before
full_ages <- data.frame(age = seq(min(combined_data$age), max(combined_data$age)))

complete_data <- full_ages %>%
  left_join(combined_data %>% select(age, incidence_per_100k_PY), by = "age") %>%
  mutate(incidence_per_100k_PY = approx(age[!is.na(incidence_per_100k_PY)],
                                       incidence_per_100k_PY[!is.na(incidence_per_100k_PY)],
                                       xout = age)$y)


# 6. Fit quadratic model to ages 80–88 (captures natural curvature)
trend_data <- complete_data %>% filter(age >= 77 & age <= 88)
quad_model <- lm(incidence_per_100k_PY ~ poly(age, 2, raw = TRUE), data = trend_data)

# 7. Extrapolate quadratic model to ages 89–100
old_extension <- data.frame(age = 89:100) %>%
  mutate(
    incidence_per_100k_PY = predict(quad_model, newdata = data.frame(age)),
    # Ensure non-negative values and smooth decline to zero at 100
    incidence_per_100k_PY = pmax(incidence_per_100k_PY, 0),
    incidence_per_100k_PY = ifelse(age == 100, 0, incidence_per_100k_PY)
  )

# 8. Young age extension (18-33) remains linear
young_anchor_age <- 36

young_avg_incidence <- complete_data %>%
  filter(age >= 35 & age <= 37) %>%
  summarise(avg_incidence = mean(incidence_per_100k_PY, na.rm = TRUE)) %>%
  pull(avg_incidence)

# Calculate slope based on anchor point
slope <- young_avg_incidence / (young_anchor_age - 18)

young_extension <- data.frame(
  age = 18:33,
  incidence_per_100k_PY = slope * (18:33 - 18)
)

# 9. Combine all data
complete_data <- bind_rows(
  young_extension,
  complete_data %>% filter(age > 33 & age < 89),
  old_extension
) %>%
  arrange(age)

# 10. Force zero incidence at age 100
complete_data$incidence_per_100k_PY[complete_data$age == 100] <- 0

# 11. Save cleaned data for simulation use
write.csv(complete_data,
          "/hpc/hers_en/pbeele/simPed/data/ALS_lifetime/lifetime_ALS_risk_processed.csv",
          row.names = FALSE)

# 12. Plot to check
plot(complete_data$age, complete_data$incidence_per_100k_PY, type = "l", col = "blue",
     xlab = "Age", ylab = "ALS incidence rate per 100,000 PY",
     main = "Processed ALS Incidence Rates by Age",
     xlim = c(18, 100), ylim = c(0, max(complete_data$incidence_per_100k_PY)))
points(raw_data$age_raw, raw_data$incidence_raw, col = "red", pch = 19, cex = 0.5)
legend("topright", legend = c("Processed (interpolated)", "Raw digitized"),
       col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 19))


## Next step to evaluate if it matches Levison's age aggregated incidence rates ## 
processed_data <- read.csv("/hpc/hers_en/pbeele/simPed/data/ALS_lifetime/lifetime_ALS_risk_processed.csv")

# Define bins and labels as per Table 1
age_bins <- c(18, 40, 50, 60, 70, 80, 200)  # 200 to include all ages >=80
age_labels <- c("18-39", "40-49", "50-59", "60-69", "70-79", "80+")

processed_data <- processed_data %>%
  mutate(age_group = cut(age, breaks = age_bins, labels = age_labels, right = FALSE, include.lowest = TRUE))

agg_incidence <- processed_data %>%
  group_by(age_group) %>%
  summarise(mean_incidence = mean(incidence_per_100k_PY, na.rm = TRUE)) %>%
  ungroup()
print(agg_incidence)
