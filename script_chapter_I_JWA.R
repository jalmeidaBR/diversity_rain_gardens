#Supplementary Material III - STATISTICAL ANALYSIS 
#PAPER: "Limited and socially unequal contributions of rain gardens to native biodiversity in a tropical megacity"; AUTHOR: Jonathan Wilson de Almeida (almeidajonathan39@usp.br)

## Install and load necessary packages
library(tibble)
library(dplyr)
library(vegan)
library(tidyr)
library(ggpubr)
library(stats)
library(bbmle)
library(MuMIn)
library(visreg)
library(DHARMa)
library(lme4)
library(effects)
library(MASS)
library(ggplot2)
library(survival)
library(car)
library(spaMM)
library(GGally)
library(performance)
library(corrplot)
library(Hmisc)
library(ggeffects)
library(gridExtra)
library(merTools)
library(glmmTMB)
library(patchwork)
library(gamlss)
library(UpSetR)
library(reshape2)
library(ggcorrplot)


# Load the data (replace "inside_data.csv" and "outside_data.csv" with your actual file names)
inside_data <- read.csv("plant_diversity_RG.csv", sep=";")
outside_data <- read.csv("plant_diversity_ADJ.csv", sep=";")

###ID_RG is represented from 1 to 30, but in the script, the IDs mainly correspond to those previously used to select the sampled areas, based on the total number of rain gardens in São Paulo. The IDs from 1 to 30 represent the initial identification in ascending order.

# Add a column to indicate whether the data is from inside or outside the garden
inside_data <- inside_data %>% mutate(area_type = "inside")
outside_data <- outside_data %>% mutate(area_type = "outside")

# Combine the datasets into a single data frame
data <- bind_rows(inside_data, outside_data)

# Step 1: Consolidate species data within each area_type and ID_RG
consolidated_data <- data %>%
  group_by(ID_RG, area_type, powo_species_name) %>%
  summarise(
    total_abund = sum(abund, na.rm = TRUE),  # Sum abundances for the same species
    total_cover = sum(cover_sp, na.rm = TRUE),  # Sum coverages for the same species
    .groups = "drop"
  )

# Step 2: Consolidate data for each ID_RG (combine inside and outside)
data_per_point <- consolidated_data %>%
  group_by(ID_RG, powo_species_name) %>%
  summarise(
    total_abund_point = sum(total_abund, na.rm = TRUE),  # Sum abundances for each species
    total_cover_point = sum(total_cover, na.rm = TRUE),  # Sum coverages for each species
    .groups = "drop"
  ) %>%
  filter(total_abund_point > 0 | total_cover_point > 0)  # Exclude species with no valid data

# Step 3: Calculate diversity indices for each point (ID_RG) separately for abundance and coverage
diversity_per_point <- data_per_point %>%
  group_by(ID_RG) %>%
  summarise(
    # Diversity indices for abundance
    shannon_abund = if (sum(total_abund_point > 0) > 0) 
      diversity(total_abund_point[total_abund_point > 0], index = "shannon") 
    else NA,
    simpson_abund = if (sum(total_abund_point > 0) > 0) 
      diversity(total_abund_point[total_abund_point > 0], index = "simpson") 
    else NA,
    n_species_abund = sum(total_abund_point > 0),  # Number of species with abundance data
    
    # Diversity indices for cover
    shannon_cover = if (sum(total_cover_point > 0) > 0) 
      diversity((total_cover_point[total_cover_point > 0] / 100), index = "shannon") 
    else NA,
    simpson_cover = if (sum(total_cover_point > 0) > 0) 
      diversity((total_cover_point[total_cover_point > 0] / 100), index = "simpson") 
    else NA,
    n_species_cover = sum(total_cover_point > 0),  # Number of species with cover data
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Total number of species for weighting
    n_species = n_species_abund + n_species_cover,
    
    # Weights for abundance and cover
    weight_abund = ifelse(n_species > 0, n_species_abund / n_species, 0),
    weight_cover = ifelse(n_species > 0, n_species_cover / n_species, 0),
    
    # Combined diversity indices
    shannon_index = coalesce(weight_abund * shannon_abund, 0) + coalesce(weight_cover * shannon_cover, 0),
    simpson_index = coalesce(weight_abund * simpson_abund, 0) + coalesce(weight_cover * simpson_cover, 0)
  ) %>%
  select(ID_RG, shannon_index, simpson_index, n_species)  # Final output columns

# Step 4: Calculate diversity indices for each area_type (inside vs. outside)
# Abundance data
abundance_data <- consolidated_data %>%
  filter(total_abund > 0)  # Filter species with abundance data

abundance_diversity <- abundance_data %>%
  group_by(ID_RG, area_type) %>%
  summarise(
    shannon_index_abund = diversity(total_abund, index = "shannon"),
    simpson_index_abund = diversity(total_abund, index = "simpson"),
    n_abund = n_distinct(powo_species_name),  # Count number of species with abundance data
    .groups = "drop"
  )

# Cover data
cover_data <- consolidated_data %>%
  filter(total_cover > 0) %>%  # Filter species with cover data
  mutate(cover_norm = total_cover / 100)  # Normalize cover to proportions

cover_diversity <- cover_data %>%
  group_by(ID_RG, area_type) %>%
  summarise(
    shannon_index_cover = diversity(cover_norm, index = "shannon"),
    simpson_index_cover = diversity(cover_norm, index = "simpson"),
    n_cover = n_distinct(powo_species_name),  # Count number of species with cover data
    .groups = "drop"
  )

# Combine diversity indices for abundance and cover with weights at the area_type level
combined_diversity <- abundance_diversity %>%
  full_join(cover_diversity, by = c("ID_RG", "area_type")) %>%
  rowwise() %>%
  mutate(
    total_species = coalesce(n_abund, 0) + coalesce(n_cover, 0),  # Total number of species
    weight_abund = coalesce(n_abund, 0) / total_species,
    weight_cover = coalesce(n_cover, 0) / total_species,
    combined_shannon = coalesce(weight_abund, 0) * coalesce(shannon_index_abund, 0) +
      coalesce(weight_cover, 0) * coalesce(shannon_index_cover, 0),
    combined_simpson = coalesce(weight_abund, 0) * coalesce(simpson_index_abund, 0) +
      coalesce(weight_cover, 0) * coalesce(simpson_index_cover, 0)
  ) %>%
  select(ID_RG, area_type, combined_shannon, combined_simpson, total_species)

# Step 5: General diversity indices by area type (inside vs. outside)
general_abundance <- abundance_data %>%
  group_by(area_type) %>%
  summarise(
    shannon_abundance = diversity(total_abund, index = "shannon"),
    simpson_abundance = diversity(total_abund, index = "simpson"),
    total_species_abund = n_distinct(powo_species_name),
    .groups = "drop"
  )

general_cover <- cover_data %>%
  group_by(area_type) %>%
  summarise(
    shannon_cover = diversity(cover_norm, index = "shannon"),
    simpson_cover = diversity(cover_norm, index = "simpson"),
    total_species_cover = n_distinct(powo_species_name),
    .groups = "drop"
  )

general_diversity <- general_abundance %>%
  full_join(general_cover, by = "area_type") %>%
  rowwise() %>%
  mutate(
    total_species = coalesce(total_species_abund, 0) + coalesce(total_species_cover, 0),
    weight_abund = coalesce(total_species_abund, 0) / total_species,
    weight_cover = coalesce(total_species_cover, 0) / total_species,
    combined_shannon_general = coalesce(weight_abund, 0) * coalesce(shannon_abundance, 0) +
      coalesce(weight_cover, 0) * coalesce(shannon_cover, 0),
    combined_simpson_general = coalesce(weight_abund, 0) * coalesce(simpson_abundance, 0) +
      coalesce(weight_cover, 0) * coalesce(simpson_cover, 0)
  ) %>%
  select(area_type, combined_shannon_general, combined_simpson_general, total_species)

###########################################################################################
#ALPHA DIVERSITY PER STATUS

# Define native and exotic species subsets based on `final_status`
native_data <- data %>% filter(final_status %in% c("native", "native_BR"))
exotic_data <- data %>% filter(final_status %in% c("cultivated", "naturalized", "invasive"))

# Function to calculate diversity for a given subset of data
calculate_diversity <- function(subset_data) {
  # Step 1: Consolidate species data within each area_type and ID_RG
  consolidated_data <- subset_data %>%
    group_by(ID_RG, area_type, powo_species_name) %>%
    summarise(
      total_abund = sum(abund, na.rm = TRUE),  # Sum abundances for the same species
      total_cover = sum(cover_sp, na.rm = TRUE),  # Sum coverages for the same species
      .groups = "drop"
    )
  
  # Step 2: Calculate diversity indices for each ID_RG separately for abundance and coverage
  diversity_per_point <- consolidated_data %>%
    group_by(ID_RG) %>%
    summarise(
      # Diversity indices for abundance
      shannon_abund = if (sum(total_abund > 0) > 0) 
        diversity(total_abund[total_abund > 0], index = "shannon") 
      else NA,
      simpson_abund = if (sum(total_abund > 0) > 0) 
        diversity(total_abund[total_abund > 0], index = "simpson") 
      else NA,
      n_species_abund = sum(total_abund > 0),  # Number of species with abundance data
      
      # Diversity indices for cover
      shannon_cover = if (sum(total_cover > 0) > 0) 
        diversity((total_cover[total_cover > 0] / 100), index = "shannon") 
      else NA,
      simpson_cover = if (sum(total_cover > 0) > 0) 
        diversity((total_cover[total_cover > 0] / 100), index = "simpson") 
      else NA,
      n_species_cover = sum(total_cover > 0),  # Number of species with cover data
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      # Total number of species for weighting
      n_species = n_species_abund + n_species_cover,
      
      # Weights for abundance and cover
      weight_abund = ifelse(n_species > 0, n_species_abund / n_species, 0),
      weight_cover = ifelse(n_species > 0, n_species_cover / n_species, 0),
      
      # Combined diversity indices
      shannon_index = coalesce(weight_abund * shannon_abund, 0) + coalesce(weight_cover * shannon_cover, 0),
      simpson_index = coalesce(weight_abund * simpson_abund, 0) + coalesce(weight_cover * simpson_cover, 0)
    ) %>%
    select(ID_RG, shannon_index, simpson_index, n_species)
  
  # Step 3: Calculate diversity indices for each area_type (inside vs. outside)
  diversity_per_area_type <- consolidated_data %>%
    group_by(ID_RG, area_type) %>%
    summarise(
      # Diversity indices for abundance
      shannon_abund = if (sum(total_abund > 0) > 0) 
        diversity(total_abund[total_abund > 0], index = "shannon") 
      else NA,
      simpson_abund = if (sum(total_abund > 0) > 0) 
        diversity(total_abund[total_abund > 0], index = "simpson") 
      else NA,
      n_species_abund = sum(total_abund > 0),
      
      # Diversity indices for cover
      shannon_cover = if (sum(total_cover > 0) > 0) 
        diversity((total_cover[total_cover > 0] / 100), index = "shannon") 
      else NA,
      simpson_cover = if (sum(total_cover > 0) > 0) 
        diversity((total_cover[total_cover > 0] / 100), index = "simpson") 
      else NA,
      n_species_cover = sum(total_cover > 0),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      # Total number of species for weighting
      n_species = n_species_abund + n_species_cover,
      
      # Weights for abundance and cover
      weight_abund = ifelse(n_species > 0, n_species_abund / n_species, 0),
      weight_cover = ifelse(n_species > 0, n_species_cover / n_species, 0),
      
      
      # Combined diversity indices
      shannon_index = coalesce(weight_abund * shannon_abund, 0) + coalesce(weight_cover * shannon_cover, 0),
      simpson_index = coalesce(weight_abund * simpson_abund, 0) + coalesce(weight_cover * simpson_cover, 0)
    ) %>%
    select(ID_RG, area_type, shannon_index, simpson_index, n_species)
  
  # Step 4: Calculate general diversity indices for inside and outside areas
  general_diversity <- consolidated_data %>%
    group_by(area_type) %>%
    summarise(
      shannon_abund = if (sum(total_abund > 0) > 0) 
        diversity(total_abund[total_abund > 0], index = "shannon") 
      else NA,
      simpson_abund = if (sum(total_abund > 0) > 0) 
        diversity(total_abund[total_abund > 0], index = "simpson") 
      else NA,
      n_species_abund = sum(total_abund > 0),
      
      shannon_cover = if (sum(total_cover > 0) > 0) 
        diversity((total_cover[total_cover > 0] / 100), index = "shannon") 
      else NA,
      simpson_cover = if (sum(total_cover > 0) > 0) 
        diversity((total_cover[total_cover > 0] / 100), index = "simpson") 
      else NA,
      n_species_cover = sum(total_cover > 0),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      n_species = n_species_abund + n_species_cover,
      weight_abund = ifelse(n_species > 0, n_species_abund / n_species, 0),
      weight_cover = ifelse(n_species > 0, n_species_cover / n_species, 0),
      shannon_index = coalesce(weight_abund * shannon_abund, 0) + coalesce(weight_cover * shannon_cover, 0),
      simpson_index = coalesce(weight_abund * simpson_abund, 0) + coalesce(weight_cover * simpson_cover, 0)
    ) %>%
    select(area_type, shannon_index, simpson_index, n_species)
  
  return(list(
    diversity_per_point = diversity_per_point,
    diversity_per_area_type = diversity_per_area_type,
    general_diversity = general_diversity
  ))
}

# Calculate diversity for native and exotic species
native_diversity <- calculate_diversity(native_data)
exotic_diversity <- calculate_diversity(exotic_data)

#################################################################################################
#BETA DIVERSITY
# Prepare the data for Bray-Curtis calculations
bray_curtis_data <- data %>%
  group_by(ID_RG, area_type, powo_species_name) %>%
  summarise(
    total_abund = sum(abund, na.rm = TRUE),  # Sum abundances
    total_cover = sum(cover_sp, na.rm = TRUE),  # Sum coverages
    .groups = "drop"
  ) %>%
  mutate(total_cover_norm = total_cover / 100)  # Normalize coverage

# Create Matrices for Abundance and Coverage
bray_curtis_data <- bray_curtis_data %>%
  group_by(ID_RG, area_type, powo_species_name) %>%
  summarise(
    total_abund = sum(total_abund, na.rm = TRUE),
    total_cover_norm = sum(total_cover_norm, na.rm = TRUE),
    .groups = "drop"
  )


####################### Function to calculate Bray-Curtis diversity for inside vs outside of each point
calculate_beta_diversity <- function(data) {
  results <- list()
  unique_points <- unique(data$ID_RG)
  
  for (point in unique_points) {
    point_data <- data %>% filter(ID_RG == point)
    
    # Filtering inside and outside data
    inside_data <- point_data %>% filter(area_type == "inside")
    outside_data <- point_data %>% filter(area_type == "outside")
    
    # Creating abundance matrix
    inside_matrix <- inside_data %>%
      select(powo_species_name, total_abund) %>%
      spread(powo_species_name, total_abund, fill = 0)
    
    outside_matrix <- outside_data %>%
      select(powo_species_name, total_abund) %>%
      spread(powo_species_name, total_abund, fill = 0)
    
    # Guarantee consistant collums
    common_columns <- union(colnames(inside_matrix), colnames(outside_matrix))
    
    inside_matrix <- inside_matrix %>%
      add_column(!!!setNames(rep(0, length(setdiff(common_columns, colnames(inside_matrix)))), 
                             setdiff(common_columns, colnames(inside_matrix)))) %>%
      select(all_of(common_columns))
    
    outside_matrix <- outside_matrix %>%
      add_column(!!!setNames(rep(0, length(setdiff(common_columns, colnames(outside_matrix)))), 
                             setdiff(common_columns, colnames(outside_matrix)))) %>%
      select(all_of(common_columns))
    
    # Combine matrix
    abund_matrix <- rbind(inside_matrix, outside_matrix)
    rownames(abund_matrix) <- c("inside", "outside")
    
    # Calculate Bray-Curtis for abundance
    bray_abund <- vegdist(as.matrix(abund_matrix), method = "bray")
    
    # Repeat the process to cover
    inside_matrix <- inside_data %>%
      select(powo_species_name, total_cover_norm) %>%
      spread(powo_species_name, total_cover_norm, fill = 0)
    
    outside_matrix <- outside_data %>%
      select(powo_species_name, total_cover_norm) %>%
      spread(powo_species_name, total_cover_norm, fill = 0)
    
    inside_matrix <- inside_matrix %>%
      add_column(!!!setNames(rep(0, length(setdiff(common_columns, colnames(inside_matrix)))), 
                             setdiff(common_columns, colnames(inside_matrix)))) %>%
      select(all_of(common_columns))
    
    outside_matrix <- outside_matrix %>%
      add_column(!!!setNames(rep(0, length(setdiff(common_columns, colnames(outside_matrix)))), 
                             setdiff(common_columns, colnames(outside_matrix)))) %>%
      select(all_of(common_columns))
    
    cover_matrix <- rbind(inside_matrix, outside_matrix)
    rownames(cover_matrix) <- c("inside", "outside")
    
    bray_cover <- vegdist(as.matrix(cover_matrix), method = "bray")
    
    # Combine Bray-Curtis indices
    n_species_abund <- ncol(abund_matrix)
    n_species_cover <- ncol(cover_matrix)
    total_species <- n_species_abund + n_species_cover
    weight_abund <- n_species_abund / total_species
    weight_cover <- n_species_cover / total_species
    bray_combined <- weight_abund * as.numeric(bray_abund) + weight_cover * as.numeric(bray_cover)
    
    results[[as.character(point)]] <- list(
      bray_abund = as.numeric(bray_abund),
      bray_cover = as.numeric(bray_cover),
      bray_combined = bray_combined
    )
  }
  
  return(results)
}

# Calculate beta diversity for inside vs outside
beta_diversity_results <- calculate_beta_diversity(bray_curtis_data)

#################################################################################################
#BETA DIVERSITY STATUS

# Filtering native species
data_native <- data %>%
  filter(final_status %in% c("native", "native_BR"))

#JACCARD INDEX
presence_absence_matrix <- data_native %>%
  mutate(unique_id = paste(ID_RG, area_type, sep = "_")) %>%  # Combine ID_RG and area_type
  group_by(unique_id, powo_species_name) %>%
  summarise(presence = 1, .groups = "drop") %>%  # presence=1
  pivot_wider(names_from = powo_species_name, values_from = presence, values_fill = 0) %>%  # Convert to matrix
  column_to_rownames(var = "unique_id")  # Define unique_id as rownames

# Identify species name in the presence absence matriz (column)
species_names <- colnames(presence_absence_matrix)

# Create one row to make sure that the point 51 inside will be considered even without native data inside (only absence information)
zero_row <- data.frame(matrix(0, nrow = 1, ncol = length(species_names)))
colnames(zero_row) <- species_names
rownames(zero_row) <- "51_inside"

# Add a row
presence_absence_matrix <- rbind(presence_absence_matrix, zero_row)

# Create a function to calculate Jaccard for inside vs. outside
calculate_jaccard_inside_outside <- function(matrix) {
  # Identify the sampled points
  unique_points <- unique(sub("_.*", "", rownames(matrix)))
  
  # results in a data frame
  results <- data.frame(ID_RG = character(), jaccard_dissimilarity = numeric(), stringsAsFactors = FALSE)
  
  for (point in unique_points) {
    # Specify rows corresponding to inside and outside per point
    inside_row <- paste(point, "inside", sep = "_")
    outside_row <- paste(point, "outside", sep = "_")
    
    # Create one submatrix to add zeros to absent communities
    submatrix <- matrix[c(inside_row, outside_row), , drop = FALSE]
    
    if (!inside_row %in% rownames(matrix)) {
      
      submatrix <- rbind(submatrix, inside = rep(0, ncol(matrix)))
      rownames(submatrix)[nrow(submatrix)] <- inside_row
    }
    
    if (!outside_row %in% rownames(matrix)) {
      
      submatrix <- rbind(submatrix, outside = rep(0, ncol(matrix)))
      rownames(submatrix)[nrow(submatrix)] <- outside_row
    }
    
    # Calculate Jaccard index 
    jaccard_dist <- vegdist(as.matrix(submatrix), method = "jaccard")
    
    # Save the results in data frame
    jaccard_matrix <- as.matrix(jaccard_dist)
    jaccard_value <- jaccard_matrix[1, 2]
    results <- rbind(results, data.frame(ID_RG = point, jaccard_dissimilarity = jaccard_value))
  }
  
  return(results)
}

# Calculte jaccard to native inside vs. outside
jaccard_results <- calculate_jaccard_inside_outside(presence_absence_matrix)

#############################################################################################
#BETA DIVERSITY - JACCARD INDEX EXOTIC
presence_absence_matrix_exo <- data_exotic %>%
  mutate(unique_id = paste(ID_RG, area_type, sep = "_")) %>%  
  group_by(unique_id, powo_species_name) %>%
  summarise(presence = 1, .groups = "drop") %>%  
  pivot_wider(names_from = powo_species_name, values_from = presence, values_fill = 0) %>%  
  column_to_rownames(var = "unique_id") 

# Create a function
calculate_jaccard_inside_outside <- function(matrix) {
  
  unique_points <- unique(sub("_.*", "", rownames(matrix)))
  
  results <- data.frame(ID_RG = character(), jaccard_dissimilarity = numeric(), stringsAsFactors = FALSE)
  
  for (point in unique_points) {
    
    inside_row <- paste(point, "inside", sep = "_")
    outside_row <- paste(point, "outside", sep = "_")
    
    submatrix <- matrix[c(inside_row, outside_row), , drop = FALSE]
    
    if (!inside_row %in% rownames(matrix)) {
      
      submatrix <- rbind(submatrix, inside = rep(0, ncol(matrix)))
      rownames(submatrix)[nrow(submatrix)] <- inside_row
    }
    
    if (!outside_row %in% rownames(matrix)) {
      
      submatrix <- rbind(submatrix, outside = rep(0, ncol(matrix)))
      rownames(submatrix)[nrow(submatrix)] <- outside_row
    }
    
    jaccard_dist <- vegdist(as.matrix(submatrix), method = "jaccard")
    
    jaccard_matrix <- as.matrix(jaccard_dist)
    
    jaccard_value <- jaccard_matrix[1, 2]
    
    results <- rbind(results, data.frame(ID_RG = point, jaccard_dissimilarity = jaccard_value))
  }
  
  return(results)
}

# Calculate index
jaccard_results_exo <- calculate_jaccard_inside_outside(presence_absence_matrix_exo)
##########################################################################################################
#CALCULATE RICHNESS

####LOADING THE DATA AGAIN
dados_RG <- read.csv("plant_diversity_RG.csv", sep=";")
dados_ADJ <- read.csv("plant_diversity_ADJ.csv", sep=";")
# Add a column to indicate whether the data is from inside or outside the garden
inside_data <- inside_data %>% mutate(area_type = "inside")
outside_data <- outside_data %>% mutate(area_type = "outside")
# Combine the datasets into a single data frame
data <- bind_rows(inside_data, outside_data)


#Counting species richness in gardens
riqueza_RG <- dados_RG %>%
  group_by(ID_RG) %>%
  summarise(riqueza_jardim = n_distinct(powo_species_name))

riqueza_invasive_RG <- dados_RG %>%
  filter(final_status == "invasive") %>%
  group_by(ID_RG) %>%
  summarise(riqueza_invasive_RG = n_distinct(powo_species_name))

#Counting species richness in adjacent areas
riqueza_ADJ <- dados_ADJ %>%
  group_by(ID_RG) %>%
  summarise(riqueza_adj = n_distinct(powo_species_name))

riqueza_invasive_ADJ <- dados_ADJ %>%
  filter(final_status == "invasive") %>%
  group_by(ID_RG) %>%
  summarise(riqueza_invasive_ADJ = n_distinct(powo_species_name))

#Merge data
riqueza_final <- full_join(riqueza_RG, riqueza_ADJ, by = "ID_RG")
riqueza_final <- full_join(riqueza_final, riqueza_invasive_ADJ, by = "ID_RG")
riqueza_final <- full_join(riqueza_final, riqueza_invasive_RG, by = "ID_RG")
riqueza_final[is.na(riqueza_final)] <- 0

#Calculate richness per each sampled area
riqueza_buffer <- data %>%
  group_by(ID_RG) %>%
  summarise(riqueza_buffer = n_distinct(powo_species_name))

#Calculate richness per status and sampled area
riqueza_exo_buffer <- data %>%
  filter(status == "exotic") %>%
  group_by(ID_RG) %>%
  summarise(riqueza_exo_buffer = n_distinct(powo_species_name))

riqueza_nat_buffer <- data %>%
  filter(status == "native") %>%
  group_by(ID_RG) %>%
  summarise(riqueza_nat_buffer = n_distinct(powo_species_name))

#Calculate native richness in gardens 
riqueza_native_RG <- dados_RG %>%
  filter(status == "native") %>%
  group_by(ID_RG) %>%
  summarise(native_RG = n_distinct(powo_species_name), .groups = "drop")

#Calculate exotic richness in gardens
riqueza_exotic_RG <- dados_RG %>%
  filter(status == "exotic") %>%
  group_by(ID_RG) %>%
  summarise(exotic_RG = n_distinct(powo_species_name), .groups = "drop")

#Calculate native richness in adjacent areas
riqueza_native_ADJ <- dados_ADJ %>%
  filter(status == "native") %>%
  group_by(ID_RG) %>%
  summarise(native_ADJ = n_distinct(powo_species_name), .groups = "drop")

#Calculate exotic richness in adjacent areas
riqueza_exotic_ADJ <- dados_ADJ %>%
  filter(status == "exotic") %>%
  group_by(ID_RG) %>%
  summarise(exotic_ADJ = n_distinct(powo_species_name), .groups = "drop")

#Merge data
riqueza_type <- full_join(riqueza_exotic_ADJ, riqueza_native_RG, by = c("ID_RG")) %>%
  full_join(riqueza_exotic_RG, by = c("ID_RG")) %>%
  full_join(riqueza_native_ADJ, by = c("ID_RG")) %>%
  replace_na(list(native_RG = 0, exotic_RG = 0, native_ADJ = 0, exotic_ADJ = 0))  # Substitui NAs por 0

#ID_RG as character to merge
riqueza_type <- riqueza_type %>%
  mutate(ID_RG = as.character(ID_RG))
riqueza_final <- riqueza_final %>%
  mutate(ID_RG = as.character(ID_RG))
riqueza_nat_buffer <- riqueza_nat_buffer %>%
  mutate(ID_RG = as.character(ID_RG))
riqueza_exo_buffer <- riqueza_exo_buffer %>%
  mutate(ID_RG = as.character(ID_RG))
riqueza_buffer <- riqueza_buffer %>%
  mutate(ID_RG = as.character(ID_RG))

#MERGE DATA
riqueza_final<- full_join(riqueza_final, riqueza_type, by = "ID_RG")
riqueza_final<- full_join(riqueza_final, riqueza_nat_buffer, by = "ID_RG")
riqueza_final<- full_join(riqueza_final, riqueza_exo_buffer, by = "ID_RG")
riqueza_final<- full_join(riqueza_final, riqueza_buffer, by = "ID_RG")

###########################################################################################################
#STUDENT T-TEST 

#RICHNESS - Comparing the richness of the urban vegetation with and without rain gardens contribution
buffer <- t.test(riqueza_final$riqueza_adj, 
                 riqueza_final$riqueza_buffer, 
                 paired = TRUE)
nat <- t.test(riqueza_final$native_ADJ, 
              riqueza_final$riqueza_nat_buffer, 
              paired = TRUE)
exo <- t.test(riqueza_final$exotic_ADJ, 
              riqueza_final$riqueza_exo_buffer,
              paired = TRUE)

#DIVERSITY INDEX - Comparing the indeces inside and outside rain gardens
indices <- read.csv("teste_indices.csv", sep = ";")

# Function to aplly t-test and save specifics results as data.frame
run_paired_ttest <- function(var1, var2, indices) {
  resultado <- t.test(indices[[var1]], indices[[var2]], paired = TRUE)
  resumo <- data.frame(
    Comparação = paste(var1, "vs", var2),
    t_value = round(resultado$statistic, 3),
    p_value = round(resultado$p.value, 4),
    CI_lower = round(resultado$conf.int[1], 3),
    CI_upper = round(resultado$conf.int[2], 3)
  )
  return(resumo)
}

# List of variables pairwise
pares <- list(
  c("general_richness_RG", "general_richness_ADJ"),
  c("general_simpson_RG", "general_simpson_ADJ"),
  c("richness_native_RG", "richness_native_ADJ"),
  c("simpson_native_RG", "simpson_native_ADJ"),
  c("richness_exotic_RG", "richness_exotic_ADJ"),
  c("simpson_exotic_RG", "simpson_exotic_ADJ"),
  c("bray_curtis_native", "bray_curtis_exotic")
)

# Applying t-test to the pairwise list
resultados <- map_dfr(pares, ~run_paired_ttest(.x[1], .x[2], indices))

# Results
print(resultados)

###########################################################################################################
#NMDS AND PERMANOVA
# Load the data (replace with your actual file names)
inside_data <- read.csv("plant_diversity_RG.csv", sep=";")
outside_data <- read.csv("plant_diversity_ADJ.csv", sep=";")

# Add a column to indicate whether the data is from inside or outside the garden
inside_data <- inside_data %>% mutate(area_type = "inside")
outside_data <- outside_data %>% mutate(area_type = "outside")

# Combine the datasets into a single data frame
data <- bind_rows(inside_data, outside_data)

matrix_NMDS <- data %>%
  mutate(unique_id = paste(ID_RG, area_type, sep = "_")) %>%  # CombineID_RG and area_type
  group_by(unique_id, powo_species_name) %>%
  summarise(presence = 1, .groups = "drop") %>%  # presence=1
  pivot_wider(names_from = powo_species_name, values_from = presence, values_fill = 0) %>%  # Convert to matrix
  column_to_rownames(var = "unique_id")

# Calculate NMDS by data matrix
nmds <- metaMDS(matrix_NMDS, distance = "bray", k = 2, trymax = 100)

# Diferenciate inside e outside by colors
colors <- ifelse(grepl("inside", rownames(matrix_NMDS)), "blue", "red")

# Create labels to associate to each garden ID
labels <- sub("_inside|_outside", "", rownames(matrix_NMDS))

# Plot NMDS
ordiplot(nmds, type = "n", main = "NMDS - Jardins e Áreas Adjacentes")
points(nmds, col = colors, pch = 16)  
text(nmds, labels = labels, pos = 4, cex = 0.8)  
legend("topright", legend = c("Inside", "Outside"), col = c("blue", "red"), pch = 16)

#Trying to separate groups
ordiellipse(nmds, groups, col = c("blue", "red"), label = TRUE)

#CALCULATE PERMANOVA - testing the significance of the group differences
permanova_result <- adonis2(matrix_NMDS ~ groups, data = data.frame(groups = groups), method = "bray")
print(permanova_result)

###NMDS PER STATUS
native_inside <- inside_data %>% filter(status %in% c("native"))
native_outside <- outside_data %>% filter(status %in% c("native"))
exotic_inside <- inside_data %>% filter(status %in% c("exotic"))
exotic_outside <- outside_data %>% filter(final_status %in% c("exotic"))

# Add a column to indicate whether the data is from inside or outside the garden
native_inside <- native_inside %>% mutate(area_type = "inside")
native_outside <- native_outside %>% mutate(area_type = "outside")
exotic_inside <- exotic_inside %>% mutate(area_type = "inside")
exotic_outside <- exotic_outside %>% mutate(area_type = "outside")

# Combine the datasets into a single data frame
data_native <- bind_rows(native_inside, native_outside)
data_exotic <- bind_rows(exotic_inside, exotic_outside)

#NMDS- NATIVE SPECIES
matrix_NMDS_native <- data_native %>%
  mutate(unique_id = paste(ID_RG, area_type, sep = "_")) %>%  # Combine ID_RG and area_type
  group_by(unique_id, powo_species_name) %>%
  summarise(presence = 1, .groups = "drop") %>%  #Presence=1
  pivot_wider(names_from = powo_species_name, values_from = presence, values_fill = 0) %>%  # Convert to matrix
  column_to_rownames(var = "unique_id")

# Calculate NMDS 
nmds_native <- metaMDS(matrix_NMDS_native, distance = "bray", k = 2, trymax = 100)

################PLOT NMDS TEST
# Differenciate inside e outside by colors
colors <- ifelse(grepl("inside", rownames(matrix_NMDS_native)), "blue", "red")

# Create labels to each garden ID
labels <- sub("_inside|_outside", "", rownames(matrix_NMDS_native))

# Plot NMDS
ordiplot(nmds_native, type = "n", main = "NMDS nativas - Jardins e Áreas Adjacentes")
points(nmds_native, col = colors, pch = 16)  
text(nmds_native, labels = labels, pos = 4, cex = 0.8)  
legend("topright", legend = c("Inside", "Outside"), col = c("blue", "red"), pch = 16)

#Trying to separate groups 
ordiellipse(nmds_native, groups, col = c("blue", "red"), label = TRUE)

#PERMANOVA
permanova_result_native <- adonis2(matrix_NMDS_native ~ groups, data = data.frame(groups = groups), method = "bray")

###NMDS EXOTIC
matrix_NMDS_exotic <- data_exotic %>%
  mutate(unique_id = paste(ID_RG, area_type, sep = "_")) %>%  
  group_by(unique_id, powo_species_name) %>%
  summarise(presence = 1, .groups = "drop") %>%  
  pivot_wider(names_from = powo_species_name, values_from = presence, values_fill = 0) %>% 
  column_to_rownames(var = "unique_id")

# Calculate NMDS by data matrix
nmds_exotic <- metaMDS(matrix_NMDS_exotic, distance = "bray", k = 2, trymax = 100)

# Differenciate inside e outside by colors
colors <- ifelse(grepl("inside", rownames(matrix_NMDS_exotic)), "blue", "red")

# Labels to each garden ID
labels <- sub("_inside|_outside", "", rownames(matrix_NMDS_exotic))

################PLOT NMDS TEST
# Plota NMDS
ordiplot(nmds_exotic, type = "n", main = "NMDS exoticas - Jardins e Áreas Adjacentes")
points(nmds_exotic, col = colors, pch = 16)  # Adiciona pontos coloridos
text(nmds_exotic, labels = labels, pos = 4, cex = 0.8)  # Adiciona rótulos com ID_RG
legend("topright", legend = c("Inside", "Outside"), col = c("blue", "red"), pch = 16)

# Criate groups vector(inside/outside)
groups <- ifelse(grepl("inside", rownames(matrix_NMDS_exotic)), "Inside", "Outside")

#Trying visualize separated groups
ordiellipse(nmds_exotic, groups, col = c("blue", "red"), label = TRUE)

#PERMANOVA
permanova_result_exotic <- adonis2(matrix_NMDS_exotic ~ groups, data = data.frame(groups = groups), method = "bray")

####################################################################################################################################################################################################################
#MODELS SELECTION
#Exploratory analysis: Choosing the explanatory variables
##The distributions families were tested for all response variables using the package fitdistrplus, alternatively betareg and gamlss. We tested NORMAL, BETA, GAMMA, ZAGA, BEINF, NBI, POISSON

models_data <- read.csv("models_data.csv", sep=";")

#Testing models without correlated combinations

##MODELS ALPHA DIVERSITY
# Explanation variables list
variables <- c("area_m2","IDHM_DIST", "IDHM_R_DIST", "BHM_m3", "VHM_m3", "VEG_PROP", 
               "LOW_COVER_ARB", "MEDIUM_COVER_ARB", "COVER_HERB_SHRUB", 
               "B_PROP", "WATER_PROP", "ROAD_DIST_m", "ATL_FOREST_DIST_m", 
               "NEAR_ATL_FOREST_AREA_m2", "PARK_DIST_m", "APA_DIST_m", 
               "APA_AREA_m2", "WT_DIST_m", "WT_PER", "WT_AREA_km2", "manag_level_ADJ", "manag_level_RG")

# Response variables list
response_variables <- c("combined_simpson_inside", "combined_simpson_outside", 
                        "combined_shannon_inside", "combined_shannon_outside")

# Create a correlation matrix
cor_matrix <- cor(models_data[, variables], use = "pairwise.complete.obs")

# Indentify high correlations
tem_correlacao_alta <- function(combo, cor_matrix, limiar = 0.6) {
  if (length(combo) == 1) return(FALSE)  # allows "correlations" with one variable 
  
  sub_matrix <- cor_matrix[combo, combo]
  sub_matrix[lower.tri(sub_matrix, diag = TRUE)] <- NA  # Ignore diagonal and duplicated values
  
  return(any(abs(sub_matrix) >= limiar, na.rm = TRUE))
}

# Filtered combinations
combinations_filtradas <- list()

# 🔹 Generate combinations of 1, 2 e 3 explanation variables
for (k in 1:3) {
  comb_k <- combn(variables, k, simplify = FALSE)  
  
  # Filtering the correlations by the level (high correlated exclude)
  comb_k_filtradas <- comb_k[!sapply(comb_k, tem_correlacao_alta, cor_matrix = cor_matrix)]
  
  # final combinations
  combinations_filtradas <- c(combinations_filtradas, comb_k_filtradas)
}

# data.frame to save the results 
results_dredge <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_adj = numeric(),
  stringsAsFactors = FALSE
)

# Function to run the models with selected combinations of explanation variables
for (response_var in response_variables) {
  for (combo in combinations_filtradas) {
    
    # Model formula
    formula <- as.formula(paste(response_var, "~", paste(combo, collapse = " + ")))
    
    # Adjust the model
    model <- lm(formula, data = models_data)
    
    # AIC e R² adjusted
    aic_value <- AIC(model)
    r2_adj_value <- summary(model)$adj.r.squared
    
    # Save results
    results_dredge <- results_dredge %>%
      add_row(
        response_variable = response_var,
        model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  # Remove unnecessary spaces
        variables = paste(combo, collapse = ", "),
        AIC = aic_value,
        R2_adj = r2_adj_value
      )
  }
}

# best models
print(head(results_dredge))

######Calculating null model
for (response_var in response_variables) {
  
  # Model formula
  formula <- as.formula(paste(response_var, "~", 1))
  
  # Adjust the model
  model <- lm(formula, data = models_data)
  
  # AIC e R² adjusted
  aic_value <- AIC(model)
  r2_adj_value <- summary(model)$adj.r.squared
  
  # Save results
  results_dredge <- results_dredge %>%
    add_row(
      response_variable = response_var,
      model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  # Remove unnecessary spaces
      variables = NA,
      AIC = aic_value,
      R2_adj = r2_adj_value
    )
}

############Selection by AIC
# Selection models with ΔAIC < 2
best_models_dredge <- results_dredge %>%
  group_by(response_variable) %>%  # Grouping by response
  mutate(
    delta_AIC = AIC - min(AIC),  # Calculate ΔAIC in each group
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))  # Calculate w
  ) %>%
  filter(delta_AIC < 2) %>%  # Filtering models with ΔAIC < 2 per response variable
  arrange(response_variable, AIC) %>%  # Order by variable name and lower AIC
  ungroup()

# Counting number of models by response
best_models_summary <- best_models_dredge %>%
  group_by(response_variable) %>%
  summarise(num_models = n(), .groups = "drop")

# Selected models
print(best_models_dredge)

##########################################################################################################
##MODELS BETA DIVERSITY

#Data frame to save results
results_beta_dredge <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_pseudo_adj = numeric(),
  stringsAsFactors = FALSE
)

# Function to run the models with selected combinations of explanation variables

for (combo in combinations_filtradas) {
  
  # Model formula
  formula <- as.formula(paste("bray_combined ~", paste(combo, collapse = " + ")))
  
  # Adjust model
  model <- betareg(formula, data = models_data)
  
  # AIC e R² adjusted
  aic_value <- AIC(model)
  r2_pseudo <- PseudoR2(model, which = "McFaddenAdj")
  
  # Salving results
  results_beta_dredge <- rbind(results_beta_dredge, data.frame(
    response_variable = response_var,
    model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),
    variables = paste(combo, collapse = ", "),
    AIC = aic_value,
    R2_pseudo_adj = r2_pseudo
  ))
}

###NULL MODEL
# Calculate null model
modelo_nulo <- betareg(bray_combined ~ 1, data = models_data)

# Verify the adjustment
if (!is.null(modelo_nulo)) {
  
  # Null model AIC
  aic_nulo <- AIC(modelo_nulo)
  
  # R2
  r2_pseudo_nulo <- 0 
  
  # Create onde data.frame with null model results
  df_nulo <- data.frame(
    response_variable = response_var,  
    model = "Modelo Nulo",            
    variables = "Intercepto Apenas",   
    AIC = aic_nulo,                     
    R2_pseudo_adj = r2_pseudo_nulo    
  )
  
  # Add null model to all results 
  results_beta_dredge <- rbind(results_beta_dredge, df_nulo)
}

############MODELOS SELECTION BY AIC
# Select models with ΔAIC < 2
# Calculate ΔAIC
results_beta_dredge <- results_beta_dredge %>%
  mutate(
    delta_AIC = AIC - min(AIC),  # Calculate ΔAIC in each group
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))  # Calculate W
  )

# Select ΔAIC < 2
best_beta_dredge <- results_beta_dredge %>%
  filter(delta_AIC < 2) %>%
  arrange(AIC)  # Ordenar pelo menor AIC

# See selected models
print(best_beta_dredge)

##########################################################################################
#MODELS BY STATUS
#NORMAL DISTRIBUTION
# List of response variables
response_variables <- c("shannon_inside_exotic", "shannon_outside_exotic", 
                        "shannon_outside_native", "simpson_outside_exotic")


# Data frame to store the results
results_dredge_status <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_adj = numeric(),
  stringsAsFactors = FALSE
)

# Iterate over the response variables and the filtered combinations of explanatory variables
for (response_var in response_variables) {
  for (combo in combinations_filtradas) {
    
    # Model formula
    formula <- as.formula(paste(response_var, "~", paste(combo, collapse = " + ")))
    
    # Fit the model
    model <- lm(formula, data = models_data)
    
    # Evaluate the model (AIC and adjusted R²)
    aic_value <- AIC(model)
    r2_adj_value <- summary(model)$adj.r.squared
    
    # Save the results
    results_dredge_status <- results_dredge_status %>%
      add_row(
        response_variable = response_var,
        model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  # Removes multiple spaces
        variables = paste(combo, collapse = ", "),
        AIC = aic_value,
        R2_adj = r2_adj_value
      )
  }
}

######NULL MODELS
for (response_var in response_variables) {
  
  # Model formula
  formula <- as.formula(paste(response_var, "~", 1))
  
  # Fit the model
  model <- lm(formula, data = models_data)
  
  # Evaluate the model (AIC and adjusted R²)
  aic_value <- AIC(model)
  r2_adj_value <- summary(model)$adj.r.squared
  
  # Save the results
  results_dredge_status <- results_dredge_status %>%
    add_row(
      response_variable = response_var,
      model = "Null Model",
      variables = "Intercept Only",
      AIC = aic_value,
      R2_adj = r2_adj_value
    )
}

############MODEL DIAGNOSTICS AND ADJUSTMENTS
# Select models with ΔAIC < 2
best_models_dredge_status <- results_dredge_status %>%
  group_by(response_variable) %>%  # Group by response variable
  mutate(
    delta_AIC = AIC - min(AIC),  # Compute ΔAIC within each group
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))  # Compute Akaike weights
  ) %>%
  filter(delta_AIC < 2) %>%  # Filter models with ΔAIC < 2 within each response variable
  arrange(response_variable, AIC) %>%  # Order by response variable and lowest AIC
  ungroup()

# Count number of models by response variable
best_models_summary <- best_models_dredge_status %>%
  group_by(response_variable) %>%
  summarise(num_models = n(), .groups = "drop")


##########################################################################################################
##MODELS BETA DISTRIBUTION

# Response Variables list
response_variables <- c("simpson_outside_native", "simpson_inside_exotic")

# Data frame to save results
results_beta_dist <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_pseudo_adj = numeric(),  
  stringsAsFactors = FALSE
)

# Looping to run models with filtered combinations and response variables
for (response_var in response_variables) {
  for (combo in combinations_filtradas) {
    
    # Model formula
    formula <- as.formula(paste(response_var, "~", paste(combo, collapse = " + ")))
    
    # Adjust model
    model <- betareg(formula, data = models_data)
    
    # Calculate AIC e R² adjusted
    aic_value <- AIC(model)
    r2_pseudo <- PseudoR2(model, which = "McFaddenAdj")
    
    # Save results
    results_beta_dist <- rbind(results_beta_dist, data.frame(
      response_variable = response_var,
      model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  
      variables = paste(combo, collapse = ", "),
      AIC = aic_value,
      R2_pseudo_adj = r2_pseudo
    ))
  }
}

######NULL MODELS
for (response_var in response_variables) {
  
  # Model formula
  formula <- as.formula(paste(response_var, "~", 1))
  
  # Adjust model
  model <- betareg(formula, data = models_data)
  
  #Calculate AIC e R² adjusted
  aic_value <- AIC(model)
  r2_adj_value <- summary(model)$adj.r.squared
  
  # Save results
  results_beta_dist <- results_beta_dist %>%
    add_row(
      response_variable = response_var,
      model = "Modelo Nulo",
      variables = "Intercepto Apenas",
      AIC = aic_value,
      R2_adj = r2_adj_value
    )
}

############SELECT BY AIC
# Select by ΔAIC < 2
# Calculate ΔAIC
best_beta_dist <- results_beta_dist %>%
  group_by(response_variable) %>%  
  mutate(
    delta_AIC = AIC - min(AIC),  
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)) 
  ) %>%
  filter(delta_AIC < 2) %>%  
  arrange(response_variable, AIC) %>% 
  ungroup()


# Count number of models by response variable
best_models_summary <- best_beta_dist %>%
  group_by(response_variable) %>%
  summarise(num_models = n(), .groups = "drop")

##########################################################################################
####MODELS ZAGA DISTRIBUTION
# Create data frame to save results
r_shn_native_in <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_pseudo_adj = numeric(),  
  stringsAsFactors = FALSE
)

# Run models with filtered combinations
for (combo in combinations_filtradas) {
  
  # Model formula
  formula <- as.formula(paste("shannon_inside_native ~", paste(combo, collapse = " + ")))
  
  # Adjust ZAGA model
  model <- tryCatch(
    gamlss(formula, family = ZAGA, data = models_data, control = gamlss.control(n.cyc = 100, trace = FALSE)),
    error = function(e) NULL
  )
  
  # Calculate models metrics
  aic_value <- AIC(model)
  r2_pseudo <- 1 - (logLik(model) / logLik(update(model, . ~ 1)))
  
  # Save results
  r_shn_native_in <- rbind(r_shn_native_in, data.frame(
    response_variable = "shannon_inside_native",
    model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  
    variables = paste(combo, collapse = ", "),
    AIC = aic_value,
    R2_pseudo_adj = r2_pseudo
  ))
}

###NULL MODEL
modelo_nulo <- tryCatch(
  gamlss(shannon_inside_native ~ 1, 
         family = ZAGA, 
         data = models_data, 
         control = gamlss.control(n.cyc = 100, trace = FALSE)),
  error = function(e) NULL
)
# Verify adjustment
if (!is.null(modelo_nulo)) {
  
  # Calculate AIC
  aic_nulo <- AIC(modelo_nulo)
  
  # Define R²
  r2_pseudo_nulo <- 0 
  
  # Save results
  df_nulo <- data.frame(
    response_variable = response_var,  # Variável resposta
    model = "Modelo Nulo",             # Nome do modelo
    variables = "Intercepto Apenas",   # Sem variáveis explicativas
    AIC = aic_nulo,                     # AIC do modelo nulo
    R2_pseudo_adj = r2_pseudo_nulo      # Pseudo R² McFadden Ajustado
  )
  
  # Add null models to general models list
  r_shn_native_in <- rbind(r_shn_native_in, df_nulo)
}


############SELECT MODELS BY AIC
# SeLct models ΔAIC < 2
# Calculate ΔAIC
r_shn_native_in <- r_shn_native_in %>%
  mutate(
    delta_AIC = AIC - min(AIC), 
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)) 
  )

# Select ΔAIC < 2
best_shn_native_in<- r_shn_native_in %>%
  filter(delta_AIC < 2) %>%
  arrange(AIC)  

write.csv(best_shn_native_in, "dredge_shn_in_nat.csv", row.names = FALSE)

###############################################################################################
####MODELS DISTRIBUTION BEINF

# Response variables list
response_variables <- c("simpson_inside_native", "jaccard_exotic","jaccard_native")

# Data frame to save results
res_BEINF <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_pseudo_adj = numeric(),  
  stringsAsFactors = FALSE
)

# Run models with response variables and filtered combinations
for (response_var in response_variables) {
  for (combo in combinations_filtradas) {
    
    # Model formula
    formula <- as.formula(paste(response_var, "~", paste(combo, collapse = " + ")))
    
    # Adjustment ZAGA models
    model <- tryCatch(
      gamlss(formula, family = BEINF, data = models_data, control = gamlss.control(n.cyc = 100, trace = FALSE)),
      error = function(e) NULL
    )
    
    # Calculate metrics
    aic_value <- AIC(model)
    r2_pseudo <- 1 - (logLik(model) / logLik(update(model, . ~ 1)))
    
    # Save results
    res_BEINF <- rbind(res_BEINF, data.frame(
      response_variable = response_var,
      model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  # Remove múltiplos espaços
      variables = paste(combo, collapse = ", "),
      AIC = aic_value,
      R2_pseudo_adj = r2_pseudo
    ))
  }
}
###NULL MODELS -- RUN WITH VARIABLE SEPARATELY
# Create null model
modelo_nulo <- tryCatch(
  gamlss(simpson_inside_native ~ 1, #change response variable for each variable 
         family = BEINF, 
         data = models_data, 
         control = gamlss.control(n.cyc = 100, trace = FALSE)),
  error = function(e) NULL
)

# Verify adjustment
if (!is.null(modelo_nulo)) {
  
  # Calculate AIC
  aic_nulo <- AIC(modelo_nulo)
  
  # R² 
  r2_pseudo_nulo <- 0 
  
  # Data.frame to save
  df_nulo <- rbind(df_nulo,data.frame(
    response_variable = "simpson_inside_native",  # Response variable name
    model = "Modelo Nulo",             
    variables = "Intercepto Apenas",   
    AIC = aic_nulo,                    
    R2_pseudo_adj = r2_pseudo_nulo     
  ))
}

# Add null model to the general model list
res_BEINF <- rbind(res_BEINF, df_nulo)

############SELECT MODELS BY AIC
# Select ΔAIC < 2
# Calculate ΔAIC
best_BEINF <- res_BEINF %>%
  group_by(response_variable) %>%  
  mutate(
    delta_AIC = AIC - min(AIC),  
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))  
  ) %>%
  filter(delta_AIC < 2) %>%  
  arrange(response_variable, AIC) %>%  
  ungroup()

# Count selected models by response variable
best_models_summary <- best_BEINF %>%
  group_by(response_variable) %>%
  summarise(num_models = n(), .groups = "drop")


#########################################################################################################
#MODELS RICHNESS

#RUN MODELS USING DREDGE FUNCTION

exp_data_dredge <- models_data %>% dplyr::select(area_m2, IDHM_DIST,IDHM_R_DIST, BHM_m3, VHM_m3,        VEG_PROP, LOW_COVER_ARB, MEDIUM_COVER_ARB,COVER_HERB_SHRUB, B_PROP, WATER_PROP, ROAD_DIST_m, ATL_FOREST_DIST_m,NEAR_ATL_FOREST_AREA_m2, PARK_DIST_m, APA_DIST_m, APA_AREA_m2, WT_DIST_m, WT_PER, WT_AREA_km2, manag_level_ADJ, manag_level_RG)


#Create a data.frame with explanation variables
data_pred_na <- na.omit(exp_data_dredge)

# Correlation analysis
cor.table <- cor(data_pred_na)
p.table <- rcorr(as.matrix(data_pred_na,type="pearson"))$P
corrplot::corrplot(cor(data_pred_na), method = 'color',addgrid.col = 'black',type = 'upper', diag = F, order = 'original',tl.col = 'black',
                   tl.srt = 90, tl.cex = 0.55, p.mat = p.table, sig.level = 0.05, insig = 'pch', pch.cex = 1.25,
                   addCoef.col='black', number.cex = 0.75, number.font = 2, number.digits = 2)


# Create correlation matrix to use on dredge
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.6) {
  if (j >= i) return(NA)  
  
  # Select column
  x <- pull(data, i)  
  y <- pull(data, j)
  
  # Correlation test
  ct <- cor.test(x, y, conf.level = conf.level)
  
  # TRUE for high correlations
  return(ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff)
}

# Function as vector to use in outer()
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Correlation matrix
smat <- outer(1:ncol(data_pred_na), 1:ncol(data_pred_na), vCorrelated, data = data_pred_na)


nm <- colnames(data_pred_na)

dimnames(smat) <- list(nm, nm)

##MODELS RICHNESS 
#RUN MODELS WITH FAMILY = POISSON SEPARATELY FOR :riqueza_jardim, riqueza_adj, native_ADJ, exotic_RG, exotic_ADJ - CHANGE VARIABLE INDIVIDUALY 
# FULL MODEL
options(na.action = "na.omit")
full_model <- glm (riqueza_jardim ~ #replace response variable (see line 1421)
                     COVER_HERB_SHRUB +
                     ROAD_DIST_m +
                     WATER_PROP +
                     WT_DIST_m + 
                     WT_AREA_km2 + 
                     WT_PER +
                     ATL_FOREST_DIST_m +
                     VHM_m3 +
                     NEAR_ATL_FOREST_AREA_m2 +
                     APA_AREA_m2 +
                     area_m2+
                     VEG_PROP+
                     IDHM_R_DIST+
                     B_PROP+
                     PARK_DIST_m+
                     BHM_m3+
                     MEDIUM_COVER_ARB+
                     LOW_COVER_ARB+
                     APA_DIST_m+
                     IDHM_DIST,
                   family = poisson, data = models_data) 

#Summary full model
summary(full_model)

#Select using dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,3), subset = smat)

candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# Dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

####MODEL RICHNESS OF NATIVES SPECIES IN THE RAIN GARDENS -> I DIDN'T USE DREDGE CAUSE  NEGATIVE BINOMIAL DISTRIBUTION RUN WITH ERROR USING GAMLSS
# Dataframe to save results
res_orn <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_pseudo_adj = numeric(),  
  stringsAsFactors = FALSE
)

#Run models with response variable and filtered combinations
for (combo in combinations_filtradas) {
  
  # Model formula
  formula <- as.formula(paste("native_RG ~", paste(combo, collapse = " + ")))
  
  # Adjust NBI
  model <- tryCatch(
    gamlss(formula, family = NBI, data = models_data))
  
  # Calculate metrics
  aic_value <- AIC(model)
  r2_pseudo <- 1 - (logLik(model) / logLik(update(model, . ~ 1)))
  
  # Save results
  res_orn <- rbind(res_orn, data.frame(
    response_variable = "native_RG",
    model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  # Remove múltiplos espaços
    variables = paste(combo, collapse = ", "),
    AIC = aic_value,
    R2_pseudo_adj = r2_pseudo
  ))
}

###NULL MODEL
modelo_nulo <- tryCatch(
  gamlss(native_RG ~ 1, 
         family = NBI, 
         data = models_data) 
)

# VerifY adjustment
if (!is.null(modelo_nulo)) {
  
  # Calculate AIC
  aic_nulo <- AIC(modelo_nulo)
  
  # R²
  r2_pseudo_nulo <- 0 
  
  # save null model in dataframe
    df_nulo <- rbind(df_nulo,data.frame(
    response_variable = "native_RG",  
    model = "Modelo Nulo",            
    variables = "Intercepto Apenas",   
    AIC = aic_nulo,                    
    R2_pseudo_adj = r2_pseudo_nulo     
  ))
}

# Add model to general models list
res_orn<- rbind(res_orn, df_nulo)

############SELECT MODELS BY AIC
# Select ΔAIC < 2
# Calculate ΔAIC
best_orn <- res_orn %>%
  group_by(response_variable) %>%  
  mutate(
    delta_AIC = AIC - min(AIC),  
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))  
  ) %>%
  filter(delta_AIC < 2) %>%  
  arrange(response_variable, AIC) %>% 
  ungroup()

# Count models selected by each response 
best_models_summary <- best_orn %>%
  group_by(response_variable) %>%
  summarise(num_models = n(), .groups = "drop")


##########################################################################################################
#TESTING CORRELATION AMONG RESPONSE VARIABLES

# Define groups to test
grupos <- list(
  exotic_outside = c("simpson_outside_exotic", "shannon_outside_exotic", "exotic_ADJ"),
  exotic_inside = c("shannon_inside_exotic", "simpson_inside_exotic", "exotic_RG"),
  native_inside = c("simpson_inside_native", "shannon_inside_native", "native_RG"),
  native_outside = c("simpson_outside_native", "shannon_outside_native", "native_ADJ"),
  community_inside = c("combined_simpson_inside", "combined_shannon_inside", "riqueza_jardim"),
  community_outside = c("combined_simpson_outside", "combined_shannon_outside", "riqueza_adj")
)

# List to save correlations
correlation_tables <- list()

# Directory to save graphs
dir.create("heatmaps", showWarnings = FALSE)

# Run correlations by group
for (nome_grupo in names(grupos)) {
  vars <- grupos[[nome_grupo]]
  
  # Create correlation matrix
  cor_matrix <- cor(models_data[, vars], use = "pairwise.complete.obs", method = "spearman")
  
  # Convert matrix to table
  cor_table <- melt(cor_matrix, varnames = c("Variável 1", "Variável 2"), value.name = "Correlação")
  
  # Remove correlation between equal variables (diagonal)
  cor_table <- cor_table %>% filter(`Variável 1` != `Variável 2`)
  
  # Add group name to the table
  cor_table$Grupo <- nome_grupo
  
  # Save table in the list
  correlation_tables[[nome_grupo]] <- cor_table
  
  # heatmap highlighting correlations > 0.6
  heatmap <- ggcorrplot(cor_matrix, method = "circle", type = "lower",
                        lab = TRUE, outline.color = "white", 
                        colors = c("blue", "white", "pink"),
                        insig = "blank",
                        legend.title = "Correlação") +
    ggtitle(paste("Mapa de Calor -", nome_grupo)) +
    theme_bw() +  # Define fundo branco
    theme(panel.grid = element_blank())  # Remove as grades de fundo
  
  # Salve heatmap
  ggsave(filename = paste0("heatmaps/", nome_grupo, "_heatmap.jpg"), plot = heatmap, width = 6, height = 5)
}

# Combine all tables in one data.frame
correlation_df <- bind_rows(correlation_tables)

# Save table
#write.xlsx(correlation_df, "cor_var_respostas.xlsx", row.names = FALSE)

####################################################################################################################################################################################################################
#SELECTION MODELS
#FINAL MODELS - SELECTED EXPLANATION VARIABLES (MOST FREQUENT IN PREVIOUS ANALYSIS)

#List of variables
exp_data_dredge <- models_data %>% dplyr::select(B_PROP, ROAD_DIST_m, APA_DIST_m, WT_DIST_m, area_m2,
                                                   COVER_HERB_SHRUB, IDHM_R_DIST, manag_level_ADJ,manag_level_RG)


#List with explanation variables
data_pred_na <- na.omit(exp_data_dredge)

#Correlation analysis
cor.table <- cor(data_pred_na)
p.table <- rcorr(as.matrix(data_pred_na,type="pearson"))$P
corrplot::corrplot(cor(data_pred_na), method = 'color',addgrid.col = 'black',type = 'upper', diag = F, order = 'original',tl.col = 'black',
                   tl.srt = 90, tl.cex = 0.55, p.mat = p.table, sig.level = 0.05, insig = 'pch', pch.cex = 1.25,
                   addCoef.col='black', number.cex = 0.75, number.font = 2, number.digits = 2)

#Correlation matrix

# Create a function to test correlation
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.6) {
  if (j >= i) return(NA)  
  
  # Select columns
  x <- pull(data, i)  
  y <- pull(data, j)
  
  # Correlation test
  ct <- cor.test(x, y, conf.level = conf.level)
  
  # TRUE if is not high correlated
  return(ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff)
}

# function as vector
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Logical correlation matrix
smat <- outer(1:ncol(data_pred_na), 1:ncol(data_pred_na), vCorrelated, data = data_pred_na)


nm <- colnames(data_pred_na)

dimnames(smat) <- list(nm, nm)

####MODELS NORMAL DISTRIBUTION
#FULL MODEL
options(na.action = "na.omit")

full_model <- lm(simpson_outside_exotic ~ B_PROP+
                   ROAD_DIST_m+
                   APA_DIST_m+
                   WT_DIST_m+
                   area_m2+
                   COVER_HERB_SHRUB+
                   IDHM_R_DIST+
                   manag_level_ADJ,
                 data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models by delta AIC<2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- lm(combined_simpson_inside ~ B_PROP+
                   ROAD_DIST_m+
                   APA_DIST_m+
                   WT_DIST_m+
                   area_m2+
                   COVER_HERB_SHRUB+
                   IDHM_R_DIST+
                   manag_level_RG,
                 data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- lm(combined_simpson_outside ~ B_PROP+
                   ROAD_DIST_m+
                   APA_DIST_m+
                   WT_DIST_m+
                   area_m2+
                   COVER_HERB_SHRUB+
                   IDHM_R_DIST+
                   manag_level_ADJ,
                 data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- glm(exotic_RG ~ B_PROP+
                    ROAD_DIST_m+
                    APA_DIST_m+
                    WT_DIST_m+
                    area_m2+
                    COVER_HERB_SHRUB+
                    IDHM_R_DIST+
                    manag_level_RG, 
                  family = poisson,
                  data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- glm(native_ADJ ~ B_PROP+
                    ROAD_DIST_m+
                    APA_DIST_m+
                    WT_DIST_m+
                    area_m2+
                    COVER_HERB_SHRUB+
                    IDHM_R_DIST+
                    manag_level_ADJ, 
                  family = poisson,
                  data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- glm(riqueza_jardim ~ B_PROP+
                    ROAD_DIST_m+
                    APA_DIST_m+
                    WT_DIST_m+
                    area_m2+
                    COVER_HERB_SHRUB+
                    IDHM_R_DIST+
                    manag_level_RG, 
                  family = poisson,
                  data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- gamlss(simpson_inside_native ~ B_PROP+
                       ROAD_DIST_m+
                       APA_DIST_m+
                       WT_DIST_m+
                       area_m2+
                       COVER_HERB_SHRUB+
                       IDHM_R_DIST+
                       manag_level_RG, 
                     family = BEINF, data = models_data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE))

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- gamlss(simpson_outside_native ~ B_PROP+
                       ROAD_DIST_m+
                       APA_DIST_m+
                       WT_DIST_m+
                       area_m2+
                       COVER_HERB_SHRUB+
                       IDHM_R_DIST+
                       manag_level_ADJ, 
                     family = BE, data = models_data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE))

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- gamlss(simpson_inside_exotic ~ B_PROP+
                       ROAD_DIST_m+
                       APA_DIST_m+
                       WT_DIST_m+
                       area_m2+
                       COVER_HERB_SHRUB+
                       IDHM_R_DIST+
                       manag_level_RG, 
                     family = BE, data = models_data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE))

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- gamlss(bray_combined ~ B_PROP+
                       ROAD_DIST_m+
                       APA_DIST_m+
                       WT_DIST_m+
                       area_m2+
                       COVER_HERB_SHRUB+
                       IDHM_R_DIST+
                       manag_level_RG + manag_level_ADJ, 
                     family = BE, data = models_data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE))

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- glm(riqueza_adj ~ B_PROP+
                    ROAD_DIST_m+
                    APA_DIST_m+
                    WT_DIST_m+
                    area_m2+
                    COVER_HERB_SHRUB+
                    IDHM_R_DIST+
                    manag_level_ADJ, 
                  family = poisson,
                  data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- glm(exotic_ADJ ~ B_PROP+
                    ROAD_DIST_m+
                    APA_DIST_m+
                    WT_DIST_m+
                    area_m2+
                    COVER_HERB_SHRUB+
                    IDHM_R_DIST+
                    manag_level_ADJ, 
                  family = poisson,
                  data = models_data)

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- gamlss(jaccard_native~ B_PROP+
                       ROAD_DIST_m+
                       APA_DIST_m+
                       WT_DIST_m+
                       area_m2+
                       COVER_HERB_SHRUB+
                       IDHM_R_DIST+
                       manag_level_ADJ+
                       manag_level_RG, 
                     family = BEINF, data = models_data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE))

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
options(na.action = "na.omit")

full_model <- gamlss(native_RG~ B_PROP+
                       ROAD_DIST_m+
                       APA_DIST_m+
                       WT_DIST_m+
                       area_m2+
                       COVER_HERB_SHRUB+
                       IDHM_R_DIST+
                       manag_level_RG, 
                     family = NBI, data = models_data,
                     control = gamlss.control(n.cyc = 100, trace = FALSE))

#Summary full model
summary(full_model)

#Selection by dredge

options(na.action = "na.pass")
candidate_models <- dredge(full_model, m.lim = c(0,2), subset = smat)

#Select models delta AIC <2
candidate_models_sub <- subset(candidate_models,delta<2)
test_get_models <- get.models(candidate_models_sub, subset = TRUE)
str(candidate_models_sub)

# dataframe with best models (ΔAIC < 2)
best_models_df <- candidate_models_sub %>% as.data.frame()

############################################################################################
######MODELOS COM VARIÁVEIS SELECIONADAS

# List of explanation variables
variables <- c("B_PROP", "ROAD_DIST_m", "APA_DIST_m", "WT_DIST_m", "area_m2",
               "COVER_HERB_SHRUB", "IDHM_R_DIST", "manag_level_ADJ","manag_level_RG")


# Correlation matrix
cor_matrix <- cor(models_data[, variables], use = "pairwise.complete.obs")

#Verify high correlations
tem_correlacao_alta <- function(combo, cor_matrix, limiar = 0.6) {
  if (length(combo) == 1) return(FALSE)  # Allows combination with one variable
  
  sub_matrix <- cor_matrix[combo, combo]
  sub_matrix[lower.tri(sub_matrix, diag = TRUE)] <- NA  
  
  return(any(abs(sub_matrix) >= limiar, na.rm = TRUE))
}

# Filtered combinations list
combinations_filtradas_S <- list()

#Generate combination with 1 or 2 variables
for (k in 1:2) {
  comb_k <- combn(variables, k, simplify = FALSE) 
  
  # Filtering combinations high correlated
  comb_k_filtradas <- comb_k[!sapply(comb_k, tem_correlacao_alta, cor_matrix = cor_matrix)]
  
  # Final combinations
  combinations_filtradas_S <- c(combinations_filtradas_S, comb_k_filtradas)
}


# Dataframe to save results
r_JAC_exo <- data.frame(
  response_variable = character(),
  model = character(),
  variables = character(),
  AIC = numeric(),
  R2_pseudo_adj = numeric(),
  stringsAsFactors = FALSE
)

#Run models with filtered combinations
for (combo in combinations_filtradas_S) {
  
  # Model formula
  formula <- as.formula(paste("jaccard_exotic ~", paste(combo, collapse = " + ")))
  
  # Adjustment
  model <- tryCatch(
    gamlss(formula, family = BEINF, data = models_data, control = gamlss.control(n.cyc = 100, trace = FALSE)),
    error = function(e) NULL
  )
  
  # Calculate metrics
  aic_value <- AIC(model)
  r2_pseudo <- 1 - (logLik(model) / logLik(update(model, . ~ 1)))
  
  # Save results
  r_JAC_exo <- rbind(r_JAC_exo, data.frame(
    response_variable = "jaccard_exotic",
    model = gsub("\\s+", " ", paste(deparse(formula), collapse = "")),  
    variables = paste(combo, collapse = ", "),
    AIC = aic_value,
    R2_pseudo_adj = r2_pseudo
  ))
}

###NULL MODEL
modelo_nulo <- tryCatch(
  gamlss(jaccard_exotic ~ 1, 
         family = BEINF, 
         data = models_data, 
         control = gamlss.control(n.cyc = 100, trace = FALSE)),
  error = function(e) NULL
)
# Verify adjustment
if (!is.null(modelo_nulo)) {
  
  # AIC
  aic_nulo <- AIC(modelo_nulo)
  
  # R²
  r2_pseudo_nulo <- 0 
  
  # Dataframe results null model
  df_nulo <- data.frame(
    response_variable = "jaccard_exotic",  
    model = "Modelo Nulo",             
    variables = "Intercepto Apenas",   
    AIC = aic_nulo,                    
    R2_pseudo_adj = r2_pseudo_nulo     
  )
  
  # Add null model to general models list
  r_JAC_exo <- rbind(r_JAC_exo, df_nulo)
}


############SELECT BY AIC
# Select ΔAIC < 2
# Calculte ΔAIC
r_JAC_exo <- r_JAC_exo %>%
  mutate(
    delta_AIC = AIC - min(AIC),  
    W = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))  
  )

# Select models ΔAIC < 2
best_JAC_exo<- r_JAC_exo %>%
  filter(delta_AIC < 2) %>%
  arrange(AIC) 

##########################################################################################################
###PLOTS MODELOS PAPER
############EXOTIC INSIDE MODELS
######1
exotic_RG_1<-glm(exotic_RG~COVER_HERB_SHRUB + IDHM_R_DIST, family=poisson, data=models_data)
summary(exotic_RG_1)

exotic_RG_2<-glm(exotic_RG~COVER_HERB_SHRUB + ROAD_DIST_m, family=poisson, data=models_data)
summary(exotic_RG_2)

### Graphic
# Predictions
pred_exotic_RG_1 <- ggpredict(exotic_RG_1, terms = c("COVER_HERB_SHRUB[all]"))

pred_exotic_RG_1_IDH  <- ggpredict(exotic_RG_1, terms = c("IDHM_R_DIST[all]"))

pred_exotic_RG_2_ROAD  <- ggpredict(exotic_RG_2, terms = c("ROAD_DIST_m[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_exotic_RG_1$conf.low, pred_exotic_RG_1_IDH $conf.low, pred_exotic_RG_2_ROAD$conf.low,  na.rm = TRUE)
y_max_riq_jardim <- max(pred_exotic_RG_1$conf.high, pred_exotic_RG_1_IDH $conf.high, pred_exotic_RG_2_ROAD$conf.high, na.rm = TRUE)
my_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 24),  
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# First term
plot_exotic_RG_1 <- ggplot(pred_exotic_RG_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Cover herbaceous and shrub", y = "Exotic richness inside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_exotic_RG_1)

# Second term
plot_exotic_RG_1_IDH <- ggplot(pred_exotic_RG_1_IDH, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "IDHM-R", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_exotic_RG_1_IDH )

#Third term
plot_exotic_RG_2_ROAD <- ggplot(pred_exotic_RG_2_ROAD, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Road distance (m)", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_exotic_RG_2_ROAD)

# Combined graphic
plot_exotic_RG_1_final<- (plot_exotic_RG_1 + plot_exotic_RG_1_IDH + plot_exotic_RG_2_ROAD) + 
  plot_layout(ncol = 3) 

print(plot_exotic_RG_1_final)

############NATIVE OUTSIDE MODELS
######2
native_ADJ_2<-glm(native_ADJ~area_m2 + manag_level_ADJ, family=poisson, data=models_data)
native_ADJ_3<-glm(native_ADJ~area_m2 + B_PROP, family=poisson, data=models_data)
native_ADJ_4<-glm(native_ADJ~area_m2 + IDHM_R_DIST, family=poisson, data=models_data)

### Graphic
# Predictions
pred_native_ADJ_2 <- ggpredict(native_ADJ_2, terms = c("manag_level_ADJ[all]"))

pred_native_ADJ_2_area  <- ggpredict(native_ADJ_2, terms = c("area_m2[all]"))

pred_native_ADJ_3 <- ggpredict(native_ADJ_3, terms = c("B_PROP[all]"))

pred_native_ADJ_4 <- ggpredict(native_ADJ_4, terms = c("IDHM_R_DIST[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_native_ADJ_2$conf.low, pred_native_ADJ_2_area$conf.low, pred_native_ADJ_4$conf.low, pred_native_ADJ_3$conf.low, na.rm = TRUE)
y_max_riq_jardim <- max(pred_native_ADJ_2$conf.high, pred_native_ADJ_2_area$conf.high, pred_native_ADJ_4$conf.high, pred_native_ADJ_3$conf.high, na.rm = TRUE)

# First term
plot_native_ADJ_2_area <- ggplot(pred_native_ADJ_2_area, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Rain garden area (m²)", y = "Native richness outside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_native_ADJ_2_area )

# Second term
plot_native_ADJ_2 <- ggplot(pred_native_ADJ_2, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Management level", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_native_ADJ_2)


# Third term
plot_native_ADJ_3 <- ggplot(pred_native_ADJ_3, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Building proportion", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_native_ADJ_3)

# Fourth term
plot_native_ADJ_4 <- ggplot(pred_native_ADJ_4, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "IDHM-R", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)+
  theme(plot.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "cm"))

plot(plot_native_ADJ_4)

# Combined graphic
plot_native_ADJ_final<- (plot_native_ADJ_2_area + plot_native_ADJ_2 +plot_native_ADJ_3 +plot_native_ADJ_4) + 
  plot_layout(ncol = 4) 

print(plot_native_ADJ_final)

############NATIVE INSIDE MODELS
######1
native_RG_1<-gamlss(native_RG~ROAD_DIST_m, family=NBI, data=models_data, control = gamlss.control(n.cyc = 100, trace = FALSE))
summary(native_RG_1)

native_RG_3<-gamlss(native_RG~COVER_HERB_SHRUB, family=NBI, data=models_data, control = gamlss.control(n.cyc = 100, trace = FALSE))

### Graphic
# Predictions
pred_native_RG_1 <- ggpredict(native_RG_1, terms = c("ROAD_DIST_m[all]"))
pred_native_RG_3 <- ggpredict(native_RG_3, terms = c("COVER_HERB_SHRUB[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_native_RG_1$conf.low, pred_native_RG_3$conf.low, na.rm = TRUE)
y_max_riq_jardim <- max(pred_native_RG_1$conf.high, pred_native_RG_3$conf.high, na.rm = TRUE)

# First term
plot_native_RG_1 <- ggplot(pred_native_RG_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Road distance (m)", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_native_RG_1)

# second term
plot_native_RG_3 <- ggplot(pred_native_RG_3, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Cover herbaceous and shrub", y = "Native richness inside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_native_RG_3)


# match graphic
plot_native_RG_final<- (plot_native_RG_3 + plot_native_RG_1) + 
  plot_layout(ncol = 2) 

print(plot_native_RG_final)

####GRÁFICOS ÍNDICE DE SIMPSON
#SIMPSON INSIDE
combined_simpson_inside_2 <-lm(combined_simpson_inside~COVER_HERB_SHRUB + ROAD_DIST_m, data=models_data)

### Graphic
# Predictions
pred_simp_RG_1 <- ggpredict(combined_simpson_inside_2, terms = c("ROAD_DIST_m[all]"))
pred_simp_RG_2 <- ggpredict(combined_simpson_inside_2, terms = c("COVER_HERB_SHRUB[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_simp_RG_1$conf.low, pred_simp_RG_2$conf.low, na.rm = TRUE)
y_max_riq_jardim <- max(pred_simp_RG_1$conf.high, pred_simp_RG_2$conf.high, na.rm = TRUE)

# First term
plot_simp_RG_1 <- ggplot(pred_simp_RG_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Road distance (m)", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_RG_1)

# second term
plot_simp_RG_2 <- ggplot(pred_simp_RG_2, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Cover herbaceous and shrub", y = "Simpson index inside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_RG_2)


# match graphic
plot_simp_RG_final<- (plot_simp_RG_2 + plot_simp_RG_1) + 
  plot_layout(ncol = 2) 

print(plot_simp_RG_final)

#SIMPSON OUTSIDE
combined_simpson_outside <- lm(combined_simpson_outside~COVER_HERB_SHRUB + WT_DIST_m, data=models_data)

### Graphic
# Predictions
pred_simp_ADJ_1 <- ggpredict(combined_simpson_outside, terms = c("WT_DIST_m[all]"))
pred_simp_ADJ_2 <- ggpredict(combined_simpson_outside, terms = c("COVER_HERB_SHRUB[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_simp_ADJ_1$conf.low, pred_simp_ADJ_2$conf.low, na.rm = TRUE)
y_max_riq_jardim <- max(pred_simp_ADJ_1$conf.high, pred_simp_ADJ_2$conf.high, na.rm = TRUE)

# First term
plot_simp_ADJ_1 <- ggplot(pred_simp_ADJ_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Watercourse distance (m)", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_ADJ_1)

# second term
plot_simp_ADJ_2 <- ggplot(pred_simp_ADJ_2, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Cover herbaceous and shrub", y = "Simpson index outside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_ADJ_2)


# match graphic
plot_simp_ADJ_final<- (plot_simp_ADJ_2 + plot_simp_ADJ_1) + 
  plot_layout(ncol = 2) 

print(plot_simp_ADJ_final)


#SIMPSON_OUTSIDE_EXOTIC
simpson_outside_exotic_1 <- lm(simpson_outside_exotic~B_PROP + COVER_HERB_SHRUB, data=models_data)

### Graphic
# Predictions
pred_simp_EXO_1 <- ggpredict(simpson_outside_exotic_1, terms = c("B_PROP[all]"))
pred_simp_EXO_2 <- ggpredict(simpson_outside_exotic_1, terms = c("COVER_HERB_SHRUB[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_simp_EXO_1$conf.low, pred_simp_EXO_2$conf.low, na.rm = TRUE)
y_max_riq_jardim <- max(pred_simp_EXO_1$conf.high, pred_simp_EXO_2$conf.high, na.rm = TRUE)

# First term
plot_simp_EXO_1 <- ggplot(pred_simp_EXO_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Building proportion", y = " ") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_EXO_1)

# second term
plot_simp_EXO_2 <- ggplot(pred_simp_EXO_2, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Cover herbaceous and shrub", y = "Exotic simpson index outside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_EXO_2)


# match graphic
plot_simp_EXO_final<- (plot_simp_EXO_2 + plot_simp_EXO_1) + 
  plot_layout(ncol = 2) 

print(plot_simp_EXO_final)


#SIMPSON INSIDE NATIVE
simpson_inside_native_1 <- gamlss(simpson_inside_native~COVER_HERB_SHRUB, family=BEINF, data=models_data, control = gamlss.control(n.cyc = 100, trace = FALSE))

### Graphic
# Predictions
pred_simp_NAT_1 <- ggpredict(simpson_inside_native_1, terms = c("COVER_HERB_SHRUB[all]"))

#FUNCTION TO ADJUST Y-AXIS BY CONFIDENCE INTERVAL
y_min_riq_jardim <- min(pred_simp_NAT_1$conf.low, na.rm = TRUE)
y_max_riq_jardim <- max(pred_simp_NAT_1$conf.high, na.rm = TRUE)

# First term
plot_simp_NAT_1 <- ggplot(pred_simp_NAT_1, aes(x = x, y = predicted)) +
  geom_line(color = "black", size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "gray") +  
  labs(x = "Cover herbaceous and shrub", y = "Native simpson index inside rain gardens") +
  my_theme +
  ylim(y_min_riq_jardim, y_max_riq_jardim)

plot(plot_simp_NAT_1)

# match graphic
plot_simp_NAT_final<- (plot_simp_NAT_1) + 
  plot_layout(ncol = 2) 

print(plot_simp_NAT_final)


######BETA DIVERSITY GRAPHS
dados <- read.csv("beta_diversity.csv", sep=";")
# Calculate the median and quartiles for each group
summary_data <- dados %>%
  group_by(group) %>%
  summarise(
    mediana = median(bray_combined, na.rm = TRUE),  # Median
    Q1 = quantile(bray_combined, 0.25, na.rm = TRUE),  # 1º quartil
    Q3 = quantile(bray_combined, 0.75, na.rm = TRUE)   # 3º quartil
  ) 

#   GRAPH
beta_diversity_range_0 <- ggplot(summary_data, aes(x = mediana, y = group)) +
  # HORIZONTAL LINES
  geom_errorbarh(aes(xmin = Q1, xmax = Q3), height = 0.2, color = "black", size = 1.2) +
  # VERTICAL BARS IN THE EXTREMITIES
  geom_segment(aes(x = Q1, xend = Q1, y = as.numeric(factor(group)) - 0.1, 
                   yend = as.numeric(factor(group)) + 0.1), color = "black", size = 1.2) +
  geom_segment(aes(x = Q3, xend = Q3, y = as.numeric(factor(group)) - 0.1, 
                   yend = as.numeric(factor(group)) + 0.1), color = "black", size = 1.2) +
  # CENTRAL POINT
  geom_point(size = 4, color = "black") +  
  # AXIS ADJUSTMENT
  labs(x = "Bray-Curtis Dissimilarity", y = " ") +
  scale_x_continuous(limits = c(0, 1)) +  # Define o eixo x de 0 a 1
  theme_minimal() +
  my_theme


# GRAPH
beta_diversity_adjusted <- ggplot(summary_data, aes(x = mediana, y = group)) +
  # HORIZONTAL LINES
  geom_errorbarh(aes(xmin = Q1, xmax = Q3), height = 0.2, color = "black", size = 1.2) +
  # VERTICAL BARS
  geom_segment(aes(x = Q1, xend = Q1, y = as.numeric(factor(group)) - 0.1, 
                   yend = as.numeric(factor(group)) + 0.1), color = "black", size = 1.2) +
  geom_segment(aes(x = Q3, xend = Q3, y = as.numeric(factor(group)) - 0.1, 
                   yend = as.numeric(factor(group)) + 0.1), color = "black", size = 1.2) +
  # CENTRAL POINT
  geom_point(size = 4, color = "black") +  
  # AXIS
  labs(x = "Bray-Curtis Dissimilarity", y = " ", title = " ") +
  theme_minimal() +
  my_theme

ggsave("beta_diversity_range_0.png", 
       plot = beta_diversity_range_0, 
       width = 10, height = 7, dpi = 300)

ggsave("beta_diversity_adjusted.png", 
       plot = beta_diversity_adjusted, 
       width = 10, height = 7, dpi = 300)

######## OVERALL RESULTS GRAPH
input <- c (
  Ornamental= 120,
  Ruderal= 86,
  Native=94,
  Exotic=112,
  Invasive=11,
  Rain.Garden=128,
  Adjacent.area=152,
  "Ornamental&Exotic&Rain.Garden"=20,
  "Ornamental&Native&Rain.Garden"=11,
  "Ruderal&Exotic&Rain.Garden"=7,
  "Ruderal&Native&Rain.Garden"=16,
  "Ornamental&Exotic&Adjacent.area"=39,
  "Ornamental&Native&Adjacent.area"=17,
  "Ruderal&Exotic&Adjacent.area"=5,
  "Ruderal&Native&Adjacent.area"=17,
  "Ornamental&Exotic&Rain.Garden&Adjacent.area"=21,
  "Ornamental&Native&&Rain.Garden&Adjacent.area"=12,
  "Ruderal&Exotic&&Rain.Garden&Adjacent.area"=20,
  "Ruderal&Native&&Rain.Garden&Adjacent.area"=21,
  "Ornamental&Native"=40,
  "Ornamental&Exotic"=80,
  "Ruderal&Native"=54,
  "Ruderal&Exotic"=32,
  "Ornamental&Invasive"= 8,
  "Ruderal&Invasive"= 3,
  "Ornamental&Invasive&Rain.Garden"= 4,
  "Ornamental&Invasive&Adjacent.area"= 5,
  "Ruderal&Invasive&Rain.Garden"= 3,
  "Ruderal&Invasive&Adjacent.area"= 2
)

# Plot
png("upset_plot.png", width = 1900, height = 1700, res = 300) # Define size and resolution
upset(fromExpression(input), 
      nsets = 7, 
      sets = c("Ruderal","Ornamental","Invasive","Exotic","Native","Adjacent.area", "Rain.Garden"),
      keep.order = TRUE,
      order.by = "freq", 
      sets.bar.color = NA,
      sets.x.label = NULL,
      mainbar.y.label = "Number of species",
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1, 
      number.angles = 0,  # Remove lables
)

dev.off()

