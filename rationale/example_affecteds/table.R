library(gridExtra)
library(grid)

# Create a data frame with your data
df <- data.frame(
  Inheritance = c(rep("Monogenic", 3), rep("Polygenic", 3)),
  Affected_1st_degr = c(0,1,0,0,0,1),
  Affected_2nd_degr = c(1,0,1,0,1,0),
  Affected_3rd_degr = c(0,0,2,0,0,1)
)

colnames(df) <- c("Inheritance",
                  "Affected \nrelatives \n1st degree",
                  "Affected \nrelatives \n2nd degree",
                  "Affected \nrelatives \n3rd degree")

# Create table grob
tg <- tableGrob(df, rows = NULL)

# Increase row heights (multiply original heights by factor)
tg$heights[2:7] <- tg$heights[2:7] * 1.5

# Write to PDF
pdf("table_output.pdf", width=5, height=5)  # square-ish size
grid.draw(tg)
dev.off()
