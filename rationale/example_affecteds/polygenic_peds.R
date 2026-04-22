# Load kinship2 package
library(kinship2)
library(dplyr)
library(pedtools)

# Plot helper function

## Function to create plot of pedigree
plot_pedigree <- function(core_ped, output_dir = ".", labels, colors, carriers, title) {
 
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
  print(carriers)
  # Save plot
  pdf(file = file.path(output_dir, title), width = 4, height = 4)

  plot(
    ped_plot,
    title = "Polygenic pedigrees",
    fill = colors,
    cex = 1,
    labs = rep("", nrow(core_ped)),
    showPedId=FALSE
  )
  dev.off()
  
  cat("Plot saved \n") 
}

# Define individuals: IDs 1-16

pedigree_data <- data.frame(
    # Assign father and mother IDs, 0 for founders (no parents)
    id = 1:16,
    # Fathers of each individual (0 if founder)
    dadid = c(0,0,0,0,0,0, 1,1,3,3,5,5, 9,9,11,11),
    # Mothers of each individual (0 if founder)
    momid = c(0,0,0,0,0,0, 2,2,4,4,6,6, 8,8,10,10),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2,1,2,1,2, 1,2,1,2,1,2, 1,2,1,2),
    relatives = c(2,2,2,2,"", "",2,1,1,2,"","","",1,3,3)
)

## Labels for clarity
labels <- rep("", nrow(pedigree_data)) 
carriers = c(1,2)

## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Phenotype_index = ifelse(id == 99, 1, 0),  # index = black
    Phenotype_monogenic = ifelse(id == 99, 1, 0),  # monogenic = red
    Phenotype_polygenic = ifelse(id %in% c(10, 13, 16), 1, 0),   # polygenic = blue
    Individual_Name = labels
  )

# Assign colors
colors <- case_when(
  pedigree_data_expanded$Phenotype_index == 1 ~ "black",
  pedigree_data_expanded$Phenotype_monogenic == 1 ~ "red",
  pedigree_data_expanded$Phenotype_polygenic == 1 ~ "blue",
  TRUE ~ "white"
)

label_vector <- setNames(labels, as.character(pedigree_data_expanded$id))

title = "polygenic_1.pdf"
plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, title=title)

## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Phenotype_index = ifelse(id == 99, 1, 0),  # index = black
    Phenotype_monogenic = ifelse(id == 99, 1, 0),  # monogenic = red
    Phenotype_polygenic = ifelse(id %in% c(13,14), 1, 0),   # polygenic = blue
    Individual_Name = labels
  )


# Assign colors
colors <- case_when(
  pedigree_data_expanded$Phenotype_index == 1 ~ "black",
  pedigree_data_expanded$Phenotype_monogenic == 1 ~ "red",
  pedigree_data_expanded$Phenotype_polygenic == 1 ~ "blue",
  TRUE ~ "white"
)

label_vector <- setNames(labels, as.character(pedigree_data_expanded$id))

title = "polygenic_2.pdf"
plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, title=title)


## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Phenotype_index = ifelse(id == 99, 1, 0),  # index = black
    Phenotype_monogenic = ifelse(id == 99, 1, 0),  # monogenic = red
    Phenotype_polygenic = ifelse(id %in% c(1,13), 1, 0),   # polygenic = blue
    Individual_Name = labels
  )


# Assign colors
colors <- case_when(
  pedigree_data_expanded$Phenotype_index == 1 ~ "black",
  pedigree_data_expanded$Phenotype_monogenic == 1 ~ "red",
  pedigree_data_expanded$Phenotype_polygenic == 1 ~ "blue",
  TRUE ~ "white"
)

label_vector <- setNames(labels, as.character(pedigree_data_expanded$id))

title = "polygenic_3.pdf"
plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, title=title)