################################################################################
#
# Load in and process the environmental data to ready it for prediction
#
################################################################################

# Load necessary libraries
library(terra)
library(openxlsx)
library(sf)
library(tidyverse)
library(stringr)
source("R Files/Helper Functions.R")

# Read in the survey data
survey_data_combined_cropped <- read_rds("Data/processed_survey_data.rds")

# Load in and name the locally-saved environmental covariate rasters
covars.name <- read.xlsx(paste0("W:/OffshoreWind/Boundaries_Covariates/",
                                "Envt_Covariates/OSW_Gridded_metadata-April-11-2026.xlsx"),
                         sheet = "Metadata_seacuke_model")
covars.name <- covars.name %>% # Change directory from drive to local
  mutate(Directory = str_replace(Directory,
                                 fixed("W:/OffshoreWind/Boundaries_Covariates"),
                                 "Data"), 
         Filename = str_replace(Filename, "_0([1-9])", "_\\1"))
covars.names <- covars.name$Name            
covars.name$FilenameAll <- paste0(covars.name$Filename, covars.name$extension)
covars.name <- covars.name %>% dplyr::select("FilenameAll", "Directory")
ras_list <- paste(covars.name$Directory, covars.name$FilenameAll, sep="/")
r <- lapply(ras_list, terra::rast)
for (w in 1:length(r)) {
  names(r[[w]]) <- covars.names[w]
}

# Grab RV survey domain shape files, subset to study area, and format
survey.strata <- st_read(paste0("R:/Science/Population Ecology Division/Shared/Ha",
                                "ddock 4X5Y/2023/Shapefiles/MaritimesRegionEcosyst",
                                "emAssessmentStrata(2014-)NAD83.shp"))
survey.strata.cropped <- survey.strata %>% 
  filter(StrataID%in%unique(survey_data_combined_cropped$StrataID)) %>%
  filter(!TYPE == 1)

# Calculate bottom temperature range raster (using 2000-2024 data)
temp_rasters_2000on <- do.call(c, r[88:387])
temp_range <- app(temp_rasters_2000on, fun = max, na.rm = TRUE) - 
  app(temp_rasters_2000on, fun = min, na.rm = TRUE)

# Calculate yearly average bottom temperature rasters (for 1993-2024)
btm_temp_rasters <- r[4:387]
year1 <- as.numeric(substr(names(btm_temp_rasters[[1]]),9,12))
year2 <- as.numeric(substr(names(btm_temp_rasters[[length(btm_temp_rasters)]]),9,12))
avg_temp_rasters <- c()
for(i in 1:length(year1:year2))
{
  temp_rasters <- do.call(c, btm_temp_rasters[(1+12*(i-1)):(12*i)])
  avg_temp <- app(temp_rasters, fun = mean, na.rm = TRUE)
  names(avg_temp) <- paste0("mean",(year1+i-1))
  avg_temp_rasters <- c(avg_temp_rasters, avg_temp)
}

# Process the depth raster (crop, project, aggregate)
processed_depth <- r[[1]] %>%
  mask(survey.strata.cropped) %>%
  crop(survey.strata.cropped) %>%
  terra::project("EPSG:32620") %>%
  aggregate(fact = 10) # Roughly 5 km resolution

# Process all other rasters, using depth as the reference
rasters_to_process <- c(r[2:3], avg_temp_rasters, temp_range)
processed_others <- list()
# Loop through each raster
for (i in seq_along(rasters_to_process)) {
  
  # Crop and mask first to reduce size
  rr <- rasters_to_process[[i]]
  rr_cropped <- rr %>%
    mask(survey.strata.cropped) %>%
    crop(survey.strata.cropped)
  
  # Project directly onto processed_depth grid and write to disk
  rr_projected <- terra::project(rr_cropped,
                                 processed_depth,
                                 method = "bilinear",
                                 filename = paste0("Data/Temp/processed_raster_", 
                                                   i, ".tif"), 
                                 overwrite = TRUE)
  
  # Add processed raster to list
  processed_others[[i]] <- rr_projected
}

# Setup a dataframe for each year
years <- unique(survey_data_combined_cropped$year)
out_list <- vector("list", length(years))
idx <- 1
for (yr in years) {
  # Grab the average bottom temperature raster for the year
  ras_idx <- which(names(do.call(c, processed_others)) == paste0("mean", yr))
  btm_temp_raster <- processed_others[[ras_idx]]
  
  # Stack rasters and convert to a dataframe
  tmp_stack <- c(processed_depth, do.call(c, processed_others[c(1:2,length(processed_others))]),
                 btm_temp_raster,processed_others[length(processed_others)])
  names(tmp_stack) <- c("Depth","ShearVelocity","Phi","BtmTempRange","BtmTemp")
  df <- as.data.frame(tmp_stack, xy = TRUE, na.rm = TRUE)
  
  # Add year and store in list
  df$year <- yr
  out_list[[idx]] <- df
  idx <- idx + 1
}

# Combine data from all years and further process
raster_df <- bind_rows(out_list) %>%
  mutate(
    x = x / 1000,
    y = y / 1000,
    f_year  = factor(year),
    Depth = -Depth, # Positive depths
    LogDepth = log(Depth),
    LogShearVelocity = log(ShearVelocity),
    snowcrab = TRUE
  ) %>%
  filter(!if_any(c("Depth", "LogDepth", "Phi", "ShearVelocity", 
                   "LogShearVelocity", "BtmTemp", "BtmTempRange"), is.na))

# Scale covariates using means and SDs from survey data
vars_to_scale <- c("Depth", "LogDepth", "Phi", "ShearVelocity", 
                   "LogShearVelocity", "BtmTemp", "BtmTempRange")
scaling_stats <- survey_data_combined_cropped %>%
  summarise(across(all_of(vars_to_scale),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE))))
raster_df <- raster_df %>%
  mutate(across(all_of(vars_to_scale),
                ~ (.x - scaling_stats[[paste0(cur_column(), "_mean")]]) /
                  scaling_stats[[paste0(cur_column(), "_sd")]],
                .names = "{.col}Scaled"))

# Write to rds file
write_rds(raster_df, "Data/processed_environmental_data.rds")