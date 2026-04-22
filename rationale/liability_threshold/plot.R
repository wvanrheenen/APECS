# Load required library
library(wesanderson)

# Open 6x6 inch PDF device 
pdf("liability_threshold.pdf", width = 5, height = 5)

# Parameters
mean_liab <- 0           
sd_liab <- 1             
threshold <- 2           

# Create data
x <- seq(mean_liab - 4*sd_liab, mean_liab + 4*sd_liab, length.out = 1000)
y <- dnorm(x, mean = mean_liab, sd = sd_liab)

# Get Darjeeling1 colors
fill_color <- wes_palette("Darjeeling1")[1]  # First Darjeeling1 color (#FF4136)

plot(x, y, type = "l", lwd = 2, col = "black",
     ylab = "Population proportion", 
     xlab = "Phenotypic liability",
     main = "Liability Threshold Model for Polygenic Inheritance",
     cex.main = 1.0, cex.sub = 1, cex.lab = 0.8, cex.axis = 0.8,
     mgp = c(2, 0.7, 0), mar = c(2, 2, 2, 2))  # ← Perfect tight margins

mtext("Phenotypic liability = Genetic (G) + Non-genetic (E) value", 
      side = 3, line = -1.2, cex = 0.8)

threshold_height <- dnorm(threshold, mean = mean_liab, sd = sd_liab)
lines(x = rep(threshold, 2), y = c(0, threshold_height), col = "black", lwd = 1.5, lty = 2)

x_fill <- seq(threshold, max(x), length.out = 500)
y_fill <- dnorm(x_fill, mean = mean_liab, sd = sd_liab)
polygon(c(x_fill, rev(x_fill)), c(y_fill, rep(0, length(y_fill))), 
        col = adjustcolor(fill_color, alpha.f = 0.4), border = NA)

text(threshold + 0.5, threshold_height / 1.7, "Affected (= K)", 
     col = fill_color, adj = c(0, 0), cex = 0.7)
     
points(threshold, threshold_height, pch = 19, col = "black", cex = 1)
text(threshold + 0.15, threshold_height * 1.05, "lt = 2", 
     col = "black", cex = 0.7, adj = c(0, 0), font = 3)

dev.off()
