# Load libraries
library(kinship2)
library(dplyr)
library(pedtools)
library(wesanderson)  # Added for Wes Anderson colors

# Plot helper function

## Function to create plot of pedigree
plot_pedigree <- function(core_ped, output_dir = ".", plottitle, title) {
  ped_plot <- ped(
    id = core_ped$id, fid = core_ped$dadid, mid = core_ped$momid,
    sex = core_ped$sex, isConnected = TRUE
  )
  
  wes_cols <- wes_palette("Darjeeling1", 4, type = "discrete")
  colors <- ifelse(core_ped$Core_ped == 1, "white",
            ifelse(core_ped$Inlaws == 1, wes_cols[1],
            ifelse(core_ped$Ext_branches == 1, wes_cols[2], "white")))
  
  # 3x3 white-background PDF for Illustrator
  pdf(file.path(output_dir, plottitle), width = 3, height = 3, bg = "white")
  
  # Force white background + tight margins
  par(mar = c(1, 1, 2, 1), bg = "white")  
  
  plot(ped_plot,
       title = title,
       fill = colors,
       cex = 0.8,
       symbolsize = 0.9,
       branch = 0.4,
       packed = TRUE,
       align = c(1, 1.5),
       width = 6,
       labs = rep("", nrow(core_ped)),
       textAnnot = list(bottom = list(core_ped$Individual_Name, offset = 0.3, cex = 0.5))
  )
  dev.off()

  
  cat("Compact plot saved:", plottitle, "\n")
}



# Define individuals: IDs 1-8

pedigree_data <- data.frame(
    # Assign father and mother IDs, 0 for founders (no parents)
    id = 1:8,
    # Fathers of each individual (0 if founder)
    dadid = c(0,0, 0,1,1,0, 3,5),
    # Mothers of each individual (0 if founder)
    momid = c(0,0, 0,2,2,0, 4,6),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2, 1,2,1,2, 1,2),
    labels = c("P0", "C0", "P0_1", "C0_1", "C0_2", "P0_2", "C0_1_1", "C0_2_1")
)

## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Core_ped = 1,  # All core pedigree
    Inlaws = 0,
    Ext_branches = 0,
    Individual_Name = labels
  )

print(pedigree_data_expanded)
plot_pedigree(pedigree_data_expanded, plottitle = "core_ped.pdf", title = "Step 1 - Core Ped")


# Define individuals: IDs 1-12

pedigree_data <- data.frame(
    # Assign father and mother IDs, 0 for founders (no parents)
    id = 1:12,
    # Fathers of each individual (0 if founder)
    dadid = c(0,0,0,0,0,0, 1,3,3,5, 7,9),
    # Mothers of each individual (0 if founder)
    momid = c(0,0,0,0,0,0, 2,4,4,6, 8,10),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2,1,2,1,2, 1,2,1,2, 1,2),
    labels = c("P0_1_p", "P0_1_m", "P0", "C0", "P0_2_p", "P0_2_m", "P0_1", "C0_1", "C0_2", "P0_2", "C0_1_1", "C0_2_1")
)

## Create coloring
pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Core_ped = ifelse(labels %in% c("P0", "C0", "P0_1", "C0_1", "C0_2", "P0_2", "C0_1_1", "C0_2_1"), 1, 0),  
    Inlaws = ifelse(labels %in% c("P0_1_p", "P0_1_m", "P0_2_p", "P0_2_m"), 1, 0),  
    Ext_branches = ifelse(id == 99, 1, 0), 
    Individual_Name = labels
  )

pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Core_ped = labels %in% c("P0", "C0", "P0_1", "C0_1", "C0_2", "P0_2", "C0_1_1", "C0_2_1"),
    Inlaws = labels %in% c("P0_1_p", "P0_1_m", "P0_2_p", "P0_2_m"),
    Ext_branches = 0,
    Individual_Name = labels
  )

print(pedigree_data_expanded)
plot_pedigree(pedigree_data_expanded, plottitle = "inlaws.pdf", title = "Step 2 - In-laws")

# Define individuals: IDs 1-12

pedigree_data <- data.frame(
    # Assign father and mother IDs, 0 for founders (no parents)
    id = 1:15,
    # Fathers of each individual (0 if founder)
    dadid = c(0,0,0,0,0,0, 1,3,3,5,5,0, 7,9,11),
    # Mothers of each individual (0 if founder)
    momid = c(0,0,0,0,0,0, 2,4,4,6,6,0, 8,10,12),
    # Sex vector: 1 = male, 2 = female
    sex = c(1,2,1,2,1,2, 1,2,1,2,1,2, 1,2,1),
    labels = c("P0_1_p", "P0_1_m", "P0", "C0", "P0_2_p", "P0_2_m", "P0_1", "C0_1", "C0_2", "P0_2", "P0_2_m_1", "PP0_2_m_1", "C0_1_1", "C0_2_1", "P0_2_m_1_1")
)

pedigree_data_expanded <- pedigree_data %>%
  mutate(
    Core_ped = labels %in% c("P0", "C0", "P0_1", "C0_1", "C0_2", "P0_2", "C0_1_1", "C0_2_1"),
    Inlaws = labels %in% c("P0_1_p", "P0_1_m", "P0_2_p", "P0_2_m"),
    Ext_branches = labels %in% c("P0_2_m_1", "PP0_2_m_1", "P0_2_m_1_1"),
    Individual_Name = labels
  )

print(pedigree_data_expanded)
plot_pedigree(pedigree_data_expanded, plottitle = "ext_branches.pdf", title = "Step 3 - External Branches")