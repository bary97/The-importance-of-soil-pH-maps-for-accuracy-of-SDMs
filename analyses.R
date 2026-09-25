library(terra)
library(data.table)
library(dismo)
library(raster)
library(sf)
library(RCzechia)

# STUDY AREA SETUP
# Define the target region for this run. 
# Change to "Europe", "Austria", or "Czech.Republic"
target_region <- "Europe"



# 1. LOAD DATA AND R FUNCTIONS
# Load data table containing species occurrences (column 'species') and study area (column 'study_area'), decimalLongitude and decimalLatitude
all_data <- fread("./data/data_species_SDM_20260812.csv") 

# Load raster stack of environmental variables
env_predictors <- rast("./data/rasters/climate_only_Europe.tif") # variables from CHELSA (bio1, bio4, bio12, bio15) at spatial resolution of ~1 km
# Other options
# LUCAS <- rast("./data/rasters/climate_LUCAS_Europe.tif") # variables from CHELSA and soil pH from LUCAS at spatial resolution of 500 m
# SoilGrids <- rast("./data/rasters/climate_SoilGrids_Europe.tif") # variables from CHELSA and soil pH from SoilGrids at spatial resolution of 250 m

# RASTER STACK CROPPING (Only if target_region is NOT Europe)
if (target_region != "Europe") {
  cat("\nTarget region is", target_region, "- cropping environmental rasters...\n")
  
  # Load the shapefile for the specific country
  if (target_region == "Austria") {
    country_border <- read_sf("./data/shp/austria_administrative_boundaries_national_polygon.shp")
  } else if (target_region == "Czech.Republic") {
    country_border <- republika() # Adjust path if needed
  }
  
  # Crop and mask the rasters exactly to the country borders
  env_predictors <- terra::crop(env_predictors, country_border) 
  env_predictors <- terra::mask(env_predictors, country_border)
}

# Load custom functions for filtering and cross-validation
source("envi_filtering_occ_bg.R")
source("cross_val.R")



# 2. RANDOMISED pH GENERATION (Optional)
# NOTE: Run this block only if you are running the randomised pH model.
# Using "spatial shuffling", this step takes the real pH values, 
# and randomly shuffles them across all valid pixels.
LUCAS_values <- values(env_predictors[[5]])

# Identify indices of valid pixels (excluding NA values such as oceans or areas outside the study region)
valid_indices <- which(!is.na(LUCAS_values))
# Create a copy of the original values to store the randomized data
shuffled_values <- LUCAS_values

# Randomly shuffle only the values within the valid pixels
set.seed(42) 
shuffled_values[valid_indices] <- sample(LUCAS_values[valid_indices])

# Assign the shuffled values back to the pH raster layer
LUCAS_pH_random <- setValues(env_predictors[[5]], shuffled_values)
# Replace the original pH layer with the newly randomized pH layer in the stack
env_predictors[[5]] <- LUCAS_pH_random 



# 3. PREPARE DATA FOR ENVIRONMENTAL FILTERING
# Extract only columns with spatial coordinates
coords <- all_data[, c("decimalLongitude", "decimalLatitude")] 

# Scale CHELSA raster stack (standardize to mean = 0, sd = 1)
raster_mean <- global(rast("./data/rasters/climate_only_Europe.tif"), fun = "mean", na.rm = TRUE) 
raster_sd <- global(rast("./data/rasters/climate_only_Europe.tif"), fun = "sd", na.rm = TRUE) 
CHELSA_scaled <- (rast("./data/rasters/climate_only_Europe.tif") - raster_mean[, 1]) / raster_sd[, 1] 

# Extract scaled climate values to species occurrences (ID = F prevents adding an extra ID column)
occ_scale <- terra::extract(CHELSA_scaled, coords, ID = F) 
# Add extracted climate values to the original data table
occ_scale <- cbind(as.data.frame(all_data), occ_scale) 

# Generate candidate background points across the study area (Europe)
set.seed(123) 
bg_candidates <- terra::spatSample(
  CHELSA_scaled,  # For study areas other than Europe, use the corresponding regional raster stack
  size = 300000,
  method = "random",
  xy = TRUE,
  values = FALSE,
  na.rm = TRUE
)
# Extract scaled environmental values for background candidates
bg_env <- terra::extract(CHELSA_scaled, bg_candidates)
bg_candidates_df <- cbind(bg_candidates, bg_env)

# Define the CHELSA bioclimatic predictors used for filtering
predictor_columns <- c("bio1", "bio4", "bio12", "bio15")

# Filter occurrences only for Europe (matching the candidate background points)
occ_scale_study_area <- occ_scale[occ_scale$study_area == target_region, ]

# Create a list of unique species for the study area
species_list <- unique(occ_scale_study_area$species)



# 4. ENVIRONMENTAL FILTERING

# Create an empty list to store the filtered subsets
env_filtered_list <- list()

# Loop through each species directly
for (current_sp in species_list) {
  
  # Subset occurrence data for the current species
  sub_data <- occ_scale_study_area[occ_scale_study_area$species == current_sp, ]
  
  # Run the function for ALL species to correctly filter background points
  # (This ensures background points do not overlap with any known occurrences)
  filtered_data <- enviFilter(
    occ_data = sub_data,
    pred_cols = predictor_columns,
    bg_candidates = bg_candidates_df, 
    n_bg = 50000,     # Adjust the number of background points based on the study area        
    step = 0.1,
    seed = 700
  )
  
  # If the species has fewer than 100 records, occurrences are NOT filtered.
  # The background points, however, remain properly filtered by the function.
  if (nrow(sub_data) < 100) {
    filtered_data$occ_filtered <- sub_data  # Overwrite the filtered occurrences with the raw data (sub_data)
  }
  
  # Save the final results for the current species into the list
  env_filtered_list[[current_sp]] <- filtered_data
}



# 5. MAXENT CROSS-VALIDATION, FINAL MODEL FITTING AND MODEL PREDICTIONS
env_stack <- raster::stack(env_predictors)

# Create empty lists to store both cross-validation results and final models
maxent_cv_results <- list()
final_maxent_models <- list()

# Loop through each species
for (current_sp in species_list) {
  
  cat("\nProcessing species:", current_sp, "\n")
  
  # Extract the environmentally filtered occurrences AND background points for the current species
  sp_occ_data <- env_filtered_list[[current_sp]]$occ_filtered
  sp_bg_data  <- env_filtered_list[[current_sp]]$bg_filtered
  
  # Extract only the spatial coordinates for MaxEnt
  occ_coords <- sp_occ_data[, c("decimalLongitude", "decimalLatitude")]
  bg_coords  <- sp_bg_data[, c("x", "y")] # Ensure column names match your background dataset
  
  # ---------------------------------------------------------
  # A) CROSS-VALIDATION (For performance evaluation)
  # ---------------------------------------------------------
  set.seed(123) 
  folds <- dismo::kfold(occ_coords, k = 5) 
  p_data <- cbind(occ_coords, fold = folds)
  
  cv_result <- cv.maxent(
    p = p_data,
    env = env_stack,
    bg = bg_coords, # Using the species-specific filtered background points
    removeDuplicates = FALSE, 
    maxent.args = c(
      "addsamplestobackground=true",
      "autofeature=true",
      "betamultiplier=1"
    ),
    print.progress = TRUE
  )
  
  maxent_cv_results[[current_sp]] <- cv_result
  

  # B) FINAL MODEL (For spatial prediction)
  # Fit the final model using 100% of the occurrence data
  cat("Fitting final model on all data...\n")
  
  final_model <- dismo::maxent(
    x = env_stack,
    p = occ_coords,
    a = bg_coords, # Using the species-specific filtered background points
    removeDuplicates = FALSE,
    args = c(
      "addsamplestobackground=true",
      "autofeature=true",
      "betamultiplier=1"
    )
  )
  
  # Save the final model object for future predictions
  final_maxent_models[[current_sp]] <- final_model
}



# 6. SPATIAL PREDICTIONS
cat("\n--- STARTING SPATIAL PREDICTIONS ---\n")

# Create a directory to save the output prediction maps if it doesn't exist
if (!dir.exists("./predictions")) {
  dir.create("./predictions")
}

# Loop through each species to generate and save predictions
for (current_sp in species_list) {
  
  cat("Generating spatial prediction for:", current_sp, "...\n")
  
  # 1. Extract the fitted final MaxEnt model for the current species
  current_model <- final_maxent_models[[current_sp]]
  
  # Generate the prediction map using the environmental raster stack
  # This calculates habitat suitability (typically values from 0 to 1) across the study area
  pred_map <- raster::predict(
    object = current_model, 
    x = env_stack
  )
 
  # Define the filename for saving the raster to disk
  output_filename <- paste0("./predictions/pred_", current_sp, ".tif")
  
  # Save the prediction as a GeoTIFF file
  raster::writeRaster(
    x = pred_map, 
    filename = output_filename, 
    format = "GTiff", 
    overwrite = TRUE
  )
}
