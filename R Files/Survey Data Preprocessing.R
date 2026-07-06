################################################################################
#
# Load in and process the sea cucumber survey data (snow crab and RV survey)
#
################################################################################

# Load necessary libraries
library(terra)
library(sf)
library(tidyverse)
library(gtools)
library(openxlsx)
source("R Files/Helper Functions.R")

# Read in RV survey data and 2005-2024 snow crab survey data
RV_dat <- read.csv("W:/OffshoreWind/STSDMs/SeaCucumber/data/Upto2025/RVSurveyData/RVsurvey.2000to2025.seacukes.Feb2026.csv")
snowcrab_dat1 <- read.csv("Data/snowcrab_data_clean_Jan242022 (1).csv") %>% 
  filter(year > 2004)
snowcrab_dat2 <- read.csv("W:/OffshoreWind/STSDMs/SeaCucumber/data/Upto2025/SnowCrabSurveyData/SnowCrabSurvey_combined_cukes_2021to2025.csv") %>% 
  filter(year < 2025)

# Grab month for snow crab survey datasets
snowcrab_dat1$month <- as.numeric(substr(snowcrab_dat1$TRIP,4,5))
snowcrab_dat2$month <- as.numeric(substr(snowcrab_dat2$TRIP,4,5))

# Combine the two snow crab survey datasets, removing any tows with NA weights 
# due to a net mensuration equipment error
snowcrab_dat2 <- snowcrab_dat2 |>
  dplyr::rename(
    LONGITUDE = START_LONG,
    LATITUDE  = START_LAT
  )
missing_cols <- setdiff(names(snowcrab_dat1), names(snowcrab_dat2))
snowcrab_dat2[missing_cols] <- NA
snowcrab_dat2 <- snowcrab_dat2[, names(snowcrab_dat1)]
snowcrab_dat <- rbind(snowcrab_dat1, snowcrab_dat2) %>%
  filter(!(is.na(Code6600.wgt)))

# Remove RV survey tows with recorded numbers but 0 weights, and then set weights
# to 0 in RV survey if weights are NA
RV_dat <- RV_dat %>% filter(std.WGT > 0 | is.na(std.WGT)) %>% 
  mutate(std.WGT = ifelse(is.na(std.WGT), 0, std.WGT))

# Convert weight units for RV survey to snow crab survey units
RV_dat$std.WGT_km_per_tow <- RV_dat$std.WGT/(3.241)
RV_dat$std.WGT <- RV_dat$std.WGT_km_per_tow/(0.0124968)

# Combine the two surveys
snowcrab_dat$snowcrab <- TRUE
colnames(snowcrab_dat)[c(6,7,10)] <- c("mid.lat","mid.lon","std.WGT") 
survey_data_combined <- smartbind(RV_dat, snowcrab_dat)
survey_data_combined$snowcrab <- ifelse(is.na(survey_data_combined$snowcrab) == TRUE, 
                                        FALSE, TRUE)

# Remove tows that likely recorded sea cucumber weights in error... then convert 
# to sf
survey_data_combined <- survey_data_combined %>%
  filter(!MISSIONSET %in% c("TEM2008830_94", "NED2010027_221")) %>%
  st_as_sf(coords = c("mid.lon", "mid.lat"), crs = 4326)

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

# Add depth, phi, and shear velocity covariates to the survey tows... log-transform 
# shear velocity and depth
survey_data_combined$Depth <- -terra::extract(r[[which(sapply(r, function(x) names(x) == "depth"))]], 
                                              survey_data_combined)[,2]
survey_data_combined$Phi <- terra::extract(r[[which(sapply(r, function(x) names(x) == "grainsize"))]], 
                                           survey_data_combined)[,2]
survey_data_combined$ShearVelocity <- terra::extract(r[[which(sapply(r, function(x) names(x) == "shear_velocity"))]], 
                                                     survey_data_combined)[,2]
survey_data_combined$LogShearVelocity <- log(survey_data_combined$ShearVelocity)
survey_data_combined$LogDepth <- log(survey_data_combined$Depth)

# Calculate bottom temperature range raster and add as a covariate
temp_rasters <- do.call(c, r[88:387])
temp_range <- app(temp_rasters, fun = max, na.rm = TRUE) - 
  app(temp_rasters, fun = min, na.rm = TRUE)
survey_data_combined$BtmTempRange <- terra::extract(temp_range,survey_data_combined)[,2]

# Add bottom temperature and salinity covariates (matched to month and year)
survey_data_combined$BtmTemp <- NA
survey_data_combined$BtmSalinity <- NA
for (i in seq_along(r)) {
  raster_yr <- substr(names(r[[i]]), 9, 12)
  raster_month <- substr(names(r[[i]]), 14, 15)
  
  idx <- which(survey_data_combined$year == raster_yr &
                 survey_data_combined$month == raster_month)
  
  if (length(idx) > 0) {
    values <- terra::extract(r[[i]], survey_data_combined[idx, ])[, 2]
    
    if (i < 388) {
      survey_data_combined$BtmTemp[idx] <- values
    } else {
      survey_data_combined$BtmSalinity[idx] <- values
    }
  }
}

# Further process the survey data, and get it into UTM coordinates
survey_data_combined <- survey_data_combined %>%
  filter(Depth > 0,
    if_all(c(Phi, BtmTemp, BtmTempRange, BtmSalinity, ShearVelocity), ~ !is.na(.))) %>%
  st_transform(crs = 32620) %>%
  convert_to_km()
survey_data_combined[c("x", "y")] <- st_coordinates(survey_data_combined$geometry) # Store UTM coords

# Grab RV survey domain shape files, subset to study area, and get into UTM coords
survey.strata <- st_read(paste0("R:/Science/Population Ecology Division/Shared/Ha",
                                "ddock 4X5Y/2023/Shapefiles/MaritimesRegionEcosyst",
                                "emAssessmentStrata(2014-)NAD83.shp"))
survey.strata.cropped <- survey.strata %>% 
  filter(StrataID%in%440:483) # Scotian Shelf without BoF
survey.strata.cropped.UTM <- survey.strata.cropped %>% 
  st_transform(crs = 32620) %>% 
  convert_to_km()

# Crop to the study area, and scale the covariates
survey_data_combined_cropped <- survey_data_combined %>% 
  st_intersection(st_make_valid(survey.strata.cropped.UTM)) %>%
  mutate(across(c(LogDepth, Phi, ShearVelocity, LogShearVelocity, 
                  BtmTemp, BtmTempRange, BtmSalinity), ~ scale(.), 
                .names = "{.col}Scaled"))

# Write to rds file
write_rds(survey_data_combined_cropped, "Data/processed_survey_data.rds")