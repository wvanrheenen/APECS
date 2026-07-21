# Load kinship2 package
library(kinship2)
library(dplyr)
library(pedtools)

# Plot helper function

## Function to create plot of pedigree
plot_pedigree <- function(core_ped, output_dir = ".", labels, colors) {
 
  # Prepare data for pedigree creation
  ped_plot <- ped(
    id = core_ped$id,
    fid = core_ped$dadid,
    mid = core_ped$momid,
    sex = core_ped$sex,
    isConnected = TRUE
  )
  label_vector = rep("", nrow(core_ped))
  # Convert relatives to character, replace empty strings with NA to suppress text
  text_inside <- as.character(core_ped$relatives)
  text_inside[text_inside == ""] <- NA
  names(text_inside) <- as.character(core_ped$id)

  # Save plot
  pdf(file = file.path(output_dir, "degree_relatives_EHSG.pdf"), width = 23.9/2.54, height = 11.4/2.54)  # Slightly wider plot
  
  plot(
    ped_plot,
    title = "Degree of relatives to index patient",
    fill = colors,
    cex = 2,                          # Increase overall text and symbol size
    symbolsize = 1,                   # Increase symbol (box/circle) size only    labs = rep("", nrow(core_ped)),
    labs = rep("", nrow(core_ped)),
    textAnnot = list(
        inside = text_inside,
        cex = 30  # increase text size *inside* symbols (works like base R text())
    ),
    showPedId = FALSE,
    packed = FALSE  # Key change: fixed spacing per generation
  )
  dev.off()
  
  cat("Plot saved \n") 
}

# Define individuals: IDs 1-16

pedigree_data <- data.frame(
    id = 1:9,
    dadid = c(0,0, 0,1,1,0, 3,3,5),
    momid = c(0,0, 0,2,2,0, 4,4,6),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2, 1,2,1,2, 1,2,1),
    affected = rep(0,9),
    relatives = c(2,2, 1,1,2,"", "",1,3)
)


print(pedigree_data)

## Labels for clarity
labels <- rep("", nrow(pedigree_data)) 

## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Phenotype_index = ifelse(id == 7, 1, 0),  # index = black
    Phenotype_ALS = ifelse(id == 99, 1, 0),  # ALS = red
    Phenotype_FTD = ifelse(id == 99, 1, 0),   # FTD = blue
    Phenotype_Dem = ifelse(id == 99, 1, 0),   # Dem = green
    Individual_Name = labels
  )

# Assign colors
colors <- case_when(
  pedigree_data_expanded$Phenotype_index == 1 ~ "black",
  pedigree_data_expanded$Phenotype_ALS == 1 ~ "red",
  pedigree_data_expanded$Phenotype_FTD == 1 ~ "blue",
  pedigree_data_expanded$Phenotype_Dem == 1 ~ "green",
  TRUE ~ "white"
)

label_vector <- setNames(labels, as.character(pedigree_data_expanded$id))

plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors)
