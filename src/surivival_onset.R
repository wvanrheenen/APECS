# ALS: median 2.5 years, sdlog = 0.6
meanlog_als <- log(2.5)
sdlog_als <- 0.6

# FTD: median 6.6 years, sdlog = 0.25
meanlog_ftd <- log(6.6)
sdlog_ftd <- 0.25

# Plot ALS survival distribution
png("survival_als.png", width=800, height=600)
x_als <- seq(0, 15, by = 0.1)
y_als <- dlnorm(x_als, meanlog = meanlog_als, sdlog = sdlog_als)
plot(x_als, y_als, type = "l", lwd = 2, col = "red",
     xlab = "Years after ALS onset", ylab = "Density",
     main = "ALS Post-Onset Survival (Log-normal)")
abline(v = exp(meanlog_als), col = "red", lty = 2) # median
dev.off()

# Plot FTD survival distribution
png("survival_ftd.png", width=800, height=600)
x_ftd <- seq(0, 20, by = 0.1)
y_ftd <- dlnorm(x_ftd, meanlog = meanlog_ftd, sdlog = sdlog_ftd)
plot(x_ftd, y_ftd, type = "l", lwd = 2, col = "blue",
     xlab = "Years after FTD onset", ylab = "Density",
     main = "FTD Post-Onset Survival (Log-normal)")
abline(v = exp(meanlog_ftd), col = "blue", lty = 2) # median
dev.off()

# Quantiles for IQR calculation
z25 <- qnorm(0.25)  # -0.6744898
z75 <- qnorm(0.75)  # 0.6744898

# ALS calculations
mean_als <- exp(meanlog_als + 0.5 * sdlog_als^2)
q1_als <- exp(meanlog_als + sdlog_als * z25)
q3_als <- exp(meanlog_als + sdlog_als * z75)
iqr_als <- q3_als - q1_als

# FTD calculations
mean_ftd <- exp(meanlog_ftd + 0.5 * sdlog_ftd^2)
q1_ftd <- exp(meanlog_ftd + sdlog_ftd * z25)
q3_ftd <- exp(meanlog_ftd + sdlog_ftd * z75)
iqr_ftd <- q3_ftd - q1_ftd

# Print results
cat("ALS post-onset survival (log-normal):\n")
cat(sprintf("  Median: %.2f years\n", exp(meanlog_als)))
cat(sprintf("  Mean:   %.2f years\n", mean_als))
cat(sprintf("  Q1:     %.2f years\n", q1_als))
cat(sprintf("  Q3:     %.2f years\n", q3_als))
cat(sprintf("  IQR:    %.2f years\n\n", iqr_als))

cat("FTD post-onset survival (log-normal):\n")
cat(sprintf("  Median: %.2f years\n", exp(meanlog_ftd)))
cat(sprintf("  Mean:   %.2f years\n", mean_ftd))
cat(sprintf("  Q1:     %.2f years\n", q1_ftd))
cat(sprintf("  Q3:     %.2f years\n", q3_ftd))
cat(sprintf("  IQR:    %.2f years\n", iqr_ftd))
