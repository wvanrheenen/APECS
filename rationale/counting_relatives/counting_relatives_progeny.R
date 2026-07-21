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
  pdf(file = file.path(output_dir, "degree_relatives_progeny.pdf"), width = 3, height = 2.5)  # Slightly wider plot
  
  plot(
    ped_plot,
    title = "Verwantschap familieleden",
    fill = colors,
    cex = 0.8,
    labs = rep("", nrow(core_ped)),
    textAnnot = list(inside = text_inside),
    showPedId = FALSE,
    packed = FALSE  # Key change: fixed spacing per generation
  )
  dev.off()
  
  cat("Plot saved \n") 
}

# Define individuals: IDs 1-16

pedigree_data <- data.frame(
    id = 1:16,
    dadid = c(0,0, 0,1,1,0, 0,3,3,0,5,0, 7,7,9,11),
    momid = c(0,0, 0,2,2,0, 0,4,4,0,6,0, 8,8,10,12),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2, 1,2,1,2, 1,2,1,2,1,2, 1,2,1,2),
    affected = rep(0,16),
    relatives = c(3,3, 2,2,3,"", 1,1,2,"",4,"", "",1,3,5)
)


print(pedigree_data)

## Labels for clarity
labels <- rep("", nrow(pedigree_data)) 

## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Phenotype_index = ifelse(id == 13, 1, 0),  # index = black
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
