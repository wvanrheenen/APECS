# Load kinship2 package
library(kinship2)
library(dplyr)
library(pedtools)

# Plot helper function

## Function to create plot of pedigree
plot_pedigree <- function(core_ped, output_dir = ".", labels, colors, carriers, plottitle, title) {
 
  # Prepare data for pedigree creation
  ped_plot <- ped(
    id = core_ped$id,
    fid = core_ped$dadid,
    mid = core_ped$momid,
    sex = core_ped$sex,
    isConnected = TRUE
  )

  # Save plot
  png(file = file.path(output_dir, plottitle), width = 800, height = 400)

  plot(
    ped_plot,
    title = title,
    fill = colors,
    cex = 1.5,
    labs = rep("", nrow(core_ped)),
    textAnnot = list(bottom = list(core_ped$labels, offset = 0.5, cex = 0.8))
  )
  dev.off()
  
  cat("Plot saved \n") 
}

# Define individuals: IDs 1-8

pedigree_data <- data.frame(
    # Assign father and mother IDs, 0 for founders (no parents)
    id = 1:13,
    # Fathers of each individual (0 if founder)
    dadid = c(0,0, 0,1,1,0, 3,5,5,5,5,5,5),
    # Mothers of each individual (0 if founder)
    momid = c(0,0, 0,2,2,0, 4,6,6,6,6,6,6),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2, 1,2,1,2, 1,2,1,2,1,2,1),
    labels = c("Granddad", "Grandmom", "Dad", "Mom", "Uncle", "Aunt", "Patient\n\n", "Cousin1\n\n", "Cousin2\n\n", "Cousin3\n\n", "Cousin4\n\n", "Cousin5\n\n", "Cousin6\n\n")
)

# ## Create coloring
# pedigree_data_expanded <- pedigree_data %>%
#   mutate(
#     Healthy = ifelse(labels %in% c("Granddad", "Grandmom", "Mom", "Dad", "Uncle", "Aunt", "Cousin1\n\n"), 1, 0),  
#     Index = ifelse(labels %in% c("Patient\n\n"), 1, 0),  
#     Individual_Name = labels
#   )

# # Assign colors
# colors <- case_when(
#   pedigree_data_expanded$Healthy == 1 ~ "white",
#   pedigree_data_expanded$Index == 1 ~ "black",
#   TRUE ~ "white"
# )

# print(pedigree_data_expanded)
# print(pedigree_data_expanded$labels)

# plottitle = "example_ped.png"
# title = "Family Pedigree"
# plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, plottitle = plottitle, title=title)



# ## Create coloring
# pedigree_data_expanded <- pedigree_data %>%
#   mutate(
#     Healthy = ifelse(labels %in% c("Grandmom", "Mom", "Dad", "Uncle", "Aunt", "Cousin1\n\n"), 1, 0),  
#     Index = ifelse(labels %in% c("Patient\n\n"), 1, 0),  
#     ALS = ifelse(labels %in% c("Granddad"), 1, 0),  
#     Individual_Name = labels
#   )

# # Assign colors
# colors <- case_when(
#   pedigree_data_expanded$Healthy == 1 ~ "white",
#   pedigree_data_expanded$Index == 1 ~ "black",
#   pedigree_data_expanded$ALS == 1 ~ "red",
#   TRUE ~ "white"
# )

# print(pedigree_data_expanded)
# print(pedigree_data_expanded$labels)

# plottitle = "example_ped.png"
# title = "Family Pedigree"
# plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, plottitle = plottitle, title=title)


# ## Create coloring
# pedigree_data_expanded <- pedigree_data %>%
#   mutate(
#     Healthy = ifelse(labels %in% c("Grandmom", "Mom", "Dad", "Aunt", "Cousin1\n\n"), 1, 0),  
#     Index = ifelse(labels %in% c("Patient\n\n"), 1, 0),  
#     ALS = ifelse(labels %in% c("Granddad"), 1, 0),  
#     FTD = ifelse(labels %in% c("Uncle"), 1, 0),  
#     Individual_Name = labels
#   )

# # Assign colors
# colors <- case_when(
#   pedigree_data_expanded$Healthy == 1 ~ "white",
#   pedigree_data_expanded$Index == 1 ~ "black",
#   pedigree_data_expanded$ALS == 1 ~ "red",
#   pedigree_data_expanded$FTD == 1 ~ "pink",
#   TRUE ~ "white"
# )

# print(pedigree_data_expanded)
# print(pedigree_data_expanded$labels)

# plottitle = "example_ped.png"
# title = "Family Pedigree"
# plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, plottitle = plottitle, title=title)



## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Healthy = ifelse(labels %in% c("Grandmom", "Mom", "Dad", "Uncle", "Aunt", "Cousin1\n\n", "Cousin2\n\n", "Cousin3\n\n", "Cousin4\n\n", "Cousin5\n\n"), 1, 0),  
    Index = ifelse(labels %in% c("Patient\n\n"), 1, 0),  
    ALS = ifelse(labels %in% c("Granddad", "Cousin6\n\n"), 1, 0),  
    Individual_Name = labels
  )

# Assign colors
colors <- case_when(
  pedigree_data_expanded$Healthy == 1 ~ "white",
  pedigree_data_expanded$Index == 1 ~ "black",
  pedigree_data_expanded$ALS == 1 ~ "black",
  TRUE ~ "white"
)

print(pedigree_data_expanded)
print(pedigree_data_expanded$labels)

plottitle = "example_ped.png"
title = "Family Pedigree"
plot_pedigree(pedigree_data_expanded, labels = labels, colors = colors, carriers=carriers, plottitle = plottitle, title=title)
