################################################################################
#
# Model sea cucumber using a spatiotemporal model in sdmTMB
#
################################################################################

# Load necessary libraries
library(rnaturalearth)
library(sf)
library(terra)
library(tidyverse)
library(sdmTMB)
library(scales)
library(viridis)
library(patchwork)
library(raster)
source("R Files/Helper Functions.R")

# Set the font size for the figures
font_size <- 16

###################################Modelling####################################

# Read in the survey data and create a nonsf copy, since sf can wreck havoc with
# certain functions
survey_dat_cropped <- read_rds("Data/processed_survey_data.rds")
survey_dat_cropped_nonsf <- as.data.frame(survey_dat_cropped)

# Check the correlation between the environmental covariates
cor(survey_dat_cropped_nonsf[c("LogDepth","BtmTemp","BtmTempRange","ShearVelocity")])

# Setup and plot mesh
model_mesh <- sdmTMB::make_mesh(data = survey_dat_cropped_nonsf, 
                                xy_cols = c("x","y"), cutoff = 13)
plot(model_mesh)

# Fit the model
fit_seacuke <- sdmTMB(
  std.WGT ~ 0 + factor(year) + 
    snowcrab +
    LogDepthScaled + I(LogDepthScaled^2) +
    ShearVelocityScaled + I(ShearVelocityScaled^2) +
    BtmTempScaled + I(BtmTempScaled^2) +
    BtmTempRangeScaled + I(BtmTempRangeScaled^2),
  family = delta_lognormal(type = "poisson-link"), 
  data = survey_dat_cropped_nonsf, 
  mesh = model_mesh, spatial = "on",
  time = "year", spatiotemporal = "rw",
  silent = F, control = sdmTMBcontrol(newton_loops = 2)
) %>% run_extra_optimization(newton_loops = 2)
sanity(fit_seacuke)

# Check residuals for positive catch component
qqnorm(residuals(fit_seacuke, model = 2))
qqline(residuals(fit_seacuke, model = 2))

# Generate conditional effects plots of the covariates
plot_effects(fit_seacuke, xvars = c("LogDepthScaled","ShearVelocityScaled", 
                                    "BtmTempScaled","BtmTempRangeScaled"),
             xlab = c("Depth (m)","Shear velocity (m/s)", 
                      "Bottom temperature (°C)", "Bottom temperature range (°C)"),
             ylim = lapply(vector("list", length = 4), function(x) c(0,300)), 
             ylab = expression("Predicted relative biomass density (kg/km"^2*")"), 
             nrow = 2, ncol = 2, no_ylab = c(2,4), yr = 2024, text_size = font_size)
ggsave(paste0(getwd(),"/Figs/env_covariate_effects.jpeg"), plot=last_plot(), 
       width=9, height=7, units="in")

######################Predictions/indices using model###########################

# Grab land data for figures
land <- ne_countries(scale = "large", returnclass = "sf", 
                     continent = "North America") %>%
  st_crop(c(xmin = -75, ymin = 40, xmax = -50, ymax = 50)) %>% 
  st_transform(crs = 32620)
# Convert units from m to km
land_km <- st_geometry(land)/1000
st_crs(land_km) <- st_crs(land)
st_geometry(land) <- land_km

# Load wind and sea cucumber-relevant polygons and convert to UTM coordinates
WEA_poly <- st_read("W:/OffshoreWind/Boundaries_Covariates/WEAs/WEA_Tier_1_polygons_June_2025/Designated_WEAs_250627/Designated_WEAs_250627.shp") %>%
  st_transform(crs = 32620) %>% convert_to_km()
sea_cuke_boundaries <- st_read("Y:/Projects/OFI/BEcoME/SeaCucumber/data/OperationalSeaCucumberBoundaries/All Fishing Areas with Reserves/Fishing_Areas_2018.shp") %>% 
  convert_to_km()

# Read in the environmental data
raster_df <- read_rds("Data/processed_environmental_data.rds")

###Plot predictions across space and time###

# Make predictions across time and space (with uncertainty)
predictions <- predict(fit_seacuke, raster_df, type = "response") %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE)
prediction_sims <- predict(fit_seacuke, raster_df, type = "response", nsim = 500)
predictions$se <- apply(prediction_sims, 1, sd)

# Plot predictions for a subset of years... overlay WEA polygons
ggplot(data = predictions[predictions$year%in%seq(2000,2024,by = 6),]) + 
  geom_raster(aes(x = x, y = y, fill = est)) +
  ggtitle(expression(paste("Predicted relative biomass density (kg/",km^{2},")"))) +
  geom_sf(data = land, fill = "darkgrey") + 
  geom_sf(data = WEA_poly, fill = NA, col = "white", lwd = 0.75) + 
  theme(text = element_text(size = font_size), 
        plot.title = element_text(size = font_size), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(0.4, "cm")) +
  coord_sf(xlim = c(125,950), ylim = c(4670, 5230)) + 
  facet_wrap(~year) +
  scale_fill_gradientn("", colors = viridis(100),
                       trans = trans_new("log10", transform = function(x) log(x+1), 
                                         inverse = function(x) exp(x)-1), 
                       breaks = c(0, 100, 10000, NA), limits = c(0, NA))
ggsave(paste0(getwd(),"/Figs/spatial_predictions.jpeg"), plot=last_plot(), 
       width=6, height=4.4, units="in")

###Illustrate the distribution of predictions within each WEA###

# Crop predictions to the WEAs
pred_in_WEA <- predictions %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  st_intersection(st_make_valid(WEA_poly))

# Generate a density plot of the predictions for each WEA, all years together
ggplot(predictions, aes(x = est)) +
  geom_density(alpha = .6, fill = "grey") + 
  geom_density(data = pred_in_WEA, 
               aes(fill = WEA), alpha = .6) + 
  scale_x_continuous(trans = "log10", 
                     breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                     labels = function(x) x) +
  labs(x = expression(paste("Predicted relative biomass density (kg/",km^{2},")")),
       y = "") +
  scale_fill_hue(h = c(150, 300)) +
  facet_wrap( ~ WEA) + 
  theme(legend.position = "none", text = element_text(size = font_size),
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(getwd(),"/Figs/WEA_distribution.jpeg"), plot=last_plot(), 
       width=7.5, height=5.25, units="in")

###Plot an index for each WEA###

# Crop raster data to the WEAs
raster_df_WEA <- raster_df %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  st_intersection(st_make_valid(WEA_poly))

# Loop over WEAs to generate an index for each
area_cell <- diff(sort(unique(raster_df$x)))[1] * diff(sort(unique(raster_df$y)))[1]
index_list <- list()
for (area in unique(raster_df_WEA$WEA)) {
  
  # Subset to WEA
  df_subset <- raster_df_WEA[raster_df_WEA$WEA == area, ]
  
  # Compute index and store in list
  idx <- get_index(predict(fit_seacuke, df_subset[,,drop = T], 
                           return_tmb_object = TRUE),
                   area = area_cell, bias_correct = TRUE)
  idx$WEA <- area
  index_list[[area]] <- idx
}
# Combine into a single dataframe
index_df <- do.call(rbind, index_list)

# Plot indices, faceting by WEA
ggplot(index_df, aes(x = year, y = est/1e6, 
                        ymin = lwr/1e6, ymax = upr/1e6)) +
  geom_ribbon(alpha = 0.2, col = NA) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, NA), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Biomass index (kt)") +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        legend.text = element_text(size = font_size),
        legend.position = "bottom", 
        panel.spacing.x = unit(0.6, "cm")) + 
  facet_wrap(~WEA, scales = "free_y")
ggsave(paste0(getwd(),"/Figs/index_WEA.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

###Plot an index for the whole study area###

# Compute an index for the whole study area
total_index_df <- get_index(predict(fit_seacuke, raster_df[,,drop = T], 
                                    return_tmb_object = TRUE),
                            area = area_cell, bias_correct = TRUE)

# Plot the index
ggplot(total_index_df, aes(x = year, y = est/1e6, 
                           ymin = lwr/1e6, ymax = upr/1e6)) +
  geom_ribbon(alpha = 0.2, col = NA) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, NA), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Biomass index (kt)") +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        legend.text = element_text(size = font_size),
        legend.position = "bottom")
ggsave(paste0(getwd(),"/Figs/index_whole_area.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

####Plot an index for the reserves####

# Crop raster data to the reserves
reserves <- sea_cuke_boundaries[sea_cuke_boundaries$Type=="Reserve",]
raster_df_reserves <- raster_df %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  st_intersection(st_make_valid(reserves))

# Compute an index for the reserve areas
reserve_index_df <- get_index(predict(fit_seacuke, raster_df_reserves[,,drop = T], 
                                    return_tmb_object = TRUE), area = area_cell, 
                            bias_correct = TRUE)

# Plot the index
ggplot(reserve_index_df, aes(x = year, y = est/1e6, 
                             ymin = lwr/1e6, ymax = upr/1e6)) +
  geom_ribbon(alpha = 0.2, col = NA) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, NA), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Biomass index (kt)") +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        legend.text = element_text(size = font_size),
        legend.position = "bottom")
ggsave(paste0(getwd(),"/Figs/index_reserve.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

####Plot the percentage of total biomass in each WEA###

# Generate plot
index_df %>% 
  left_join(total_index_df, by = "year") %>%
  mutate(prop_biomass = est.x/est.y) %>%
  ggplot(aes(x = year, y = prop_biomass*100)) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, NA), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  scale_y_continuous(labels = \(x) sprintf("%.1f", x)) +
  labs(x = "Year", y = "Percentage of total biomass in WEA") +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        legend.text = element_text(size = font_size),
        legend.position = "bottom", 
        panel.spacing.x = unit(1, "cm")) + 
  facet_wrap(~WEA, scales = "free_y")
ggsave(paste0(getwd(),"/Figs/total_biomass_percentage_WEA.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

####Plot the percentage of total biomass in the reserves###

# Generate plot
reserve_index_df %>% 
  left_join(total_index_df, by = "year") %>%
  mutate(prop_biomass = est.x/est.y) %>%
  ggplot(aes(x = year, y = prop_biomass*100)) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, 40), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  scale_y_continuous(labels = \(x) sprintf("%.1f", x)) +
  labs(x = "Year", y = "Percentage of total biomass within reserves") +
  theme_classic() + geom_hline(yintercept = 30, col = "red", lwd = 1.5, 
                               linetype = "dashed") +
  theme(text = element_text(size = font_size), 
        legend.text = element_text(size = font_size),
        legend.position = "bottom", 
        panel.spacing.x = unit(1, "cm"))
ggsave(paste0(getwd(),"/Figs/total_biomass_percentage_reserves.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

###Use logbook data to determine a core habitat cutoff###

# Import the sea cucumber logbook data and convert to an sf object in UTM coordinates
# with units of km
sc_logbook_data <- read_csv("Data/sc_logbook_data.csv")
sc_logbook_sf <- st_as_sf(sc_logbook_data, coords = c("lon", "lat"),
                          crs = 4326, remove = FALSE) %>% 
  st_transform(32620) %>% convert_to_km()

# Grab the fishing area boundaries associated with at least one logbook entry
hits <- st_intersects(sea_cuke_boundaries, sc_logbook_sf)
fished_boundaries <- sea_cuke_boundaries %>%
  mutate(not_in_logbook = lengths(hits) == 0) %>%
  filter(not_in_logbook == FALSE) %>%
  filter(Type == "Fishing Area")

# Convert the model predictions to annual rasters, clip to the fishing areas above, 
# and then spatio-temporally intersect with the logbook data
pred_rasters <- lapply(unique(predictions$year), function(i) {
  d <- predictions[predictions$year == i, ]
  xy <- st_coordinates(d)
  rast(cbind(xy, est = d$est), type = "xyz") %>% 
    crop(fished_boundaries) %>% mask(fished_boundaries)
})
names(pred_rasters) <- unique(predictions$year)[order(unique(predictions$year))]
pred <- rep(NA_real_, nrow(sc_logbook_sf))
for (yr in unique(sc_logbook_sf$yr)) {
  idx <- sc_logbook_sf$yr == yr
  pts <- vect(sc_logbook_sf[idx, ])
  pred[idx] <- terra::extract(pred_rasters[[as.character(yr)]],pts)[, 2]
}
sc_logbook_sf$predicted <- pred

# Generate a dataframe storing both the clipped prediction rasters above, and the 
# raster values therein associated with logbook entries... needed for next plot
pred_rasters_df <- bind_rows(
  lapply(names(pred_rasters), function(yr) {
    data.frame(year = as.integer(yr),
               value = values(pred_rasters[[yr]]),
               source = "Entire area")})) %>% 
  na.omit()
logbook_df <- sc_logbook_sf %>%
  st_drop_geometry() %>%
  transmute(year = yr,
            est = predicted,
            source = "Logbook locations") %>%
  filter(!is.na(est))
plot_df <- bind_rows(pred_rasters_df, logbook_df)

# Generate a plot focused on the fishing areas with at least one logbook entry: 
# distribution of all predictions vs distribution of the predictions associated 
# with logbook entries... fishery focuses on higher density areas
ggplot(data = plot_df[plot_df$year%in%2012:2022,], 
       aes(x = log(est), fill = source, col = source)) + 
  geom_density(alpha = 0.3) + facet_wrap(~year, scales = "free_y") +
  labs(x = expression(paste("Log(predicted relative biomass density (kg/",km^{2},"))")),
       y = "", fill = "", col = "") +
  theme(text = element_text(size = font_size), 
        strip.text = element_text(size = font_size), 
        legend.position = "bottom")
ggsave(paste0(getwd(),"/Figs/logbook_pred_distribution.jpeg"), plot=last_plot(), 
       width=10, height=7, units="in")

# Calculate the yearly % of those logbook entries above, where the predictions are 
# greater than the thresholds of exp(6), exp(6.5), exp(7), and exp(7.5)... plot %s 
# as a boxplot for each threshold... black dots for means
plot_df %>%
  filter(source == "Logbook locations") %>%
  group_by(year) %>%
  summarize(
    pct_gt_6   = 100 * mean(log(est) > 6,   na.rm = TRUE),
    pct_gt_6_5 = 100 * mean(log(est) > 6.5, na.rm = TRUE),
    pct_gt_7   = 100 * mean(log(est) > 7,   na.rm = TRUE),
    pct_gt_7_5 = 100 * mean(log(est) > 7.5, na.rm = TRUE)
  ) %>%
  pivot_longer(-year, names_to = "threshold", values_to = "percent") %>%
  ggplot(aes(threshold, percent)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, aes(col = year)) +
  labs(x = expression(paste(
    "Biomass density threshold (kg/", km^{2}, "))")),
    y = "% of logbook entries above threshold",
    col = "Year") +
  stat_summary(fun = mean,
               geom = "point",
               shape = 16,
               size = 3,
               colour = "black") +
  scale_colour_viridis_c() +
  scale_x_discrete(labels = paste0("exp(", seq(6,7.5, by = 0.5),")")) +
  theme(text = element_text(size = font_size))
ggsave(paste0(getwd(),"/Figs/threshold_boxplots.jpeg"), plot=last_plot(), 
       width=10, height=7, units="in")

# Grab core areas based on three biomass density thresholds (i.e., exp(6), exp(6.5), 
# and exp(7))
core_maps <- predictions %>%
  mutate(core_6 = est > exp(6),
         core_6.5 = est > exp(6.5),
         core_7 = est > exp(7)) %>%
  pivot_longer(cols = starts_with("core_"),
               names_to = "threshold",
               values_to = "core")

# Grab core areas based on three biomass density thresholds (i.e., exp(6), exp(6.5), 
# and exp(7))... a pixel is core across the entire time series if it's above the 
# threshold at least once
core_maps <- predictions %>%
  sf::st_drop_geometry() %>%
  {
    df <- .
    df %>%
      dplyr::left_join(
        df %>%
          dplyr::group_by(x, y) %>%
          dplyr::summarise(
            core_6   = any(est > exp(6)),
            core_6.5 = any(est > exp(6.5)),
            core_7   = any(est > exp(7)),
            .groups = "drop"),
        by = c("x", "y"),
        copy = TRUE)
  } %>%
  pivot_longer(cols = starts_with("core_"),
               names_to = "threshold",
               values_to = "core")

# For each threshold, spatially plot the core area across the time series, zoomed 
# in on the WEAs
ggplot() + 
  geom_raster(data = core_maps %>% filter(year %in% seq(2000,2024, by = 6)),
              aes(x = x, y = y, fill = core)) +
  geom_sf(data = land, fill = "darkgrey") +
  geom_sf(data = WEA_poly, fill = NA, col = "white", lwd = 0.75) +
  facet_grid(
    threshold ~ year, switch = "y",
    labeller = labeller(threshold = function(x) {
      paste0("exp(", sub("^core_", "", x), ")")})) +
  coord_sf(xlim = c(525,875), ylim = c(4730, 5230)) + 
  labs(y = expression(paste("Biomass density threshold (kg/",km^{2},")")), 
       x = "", fill = "") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(0.4, "cm"), 
        text = element_text(size = font_size),
        strip.text = element_text(size = font_size*0.8)) +
  scale_fill_manual(values = c(`TRUE` = "#FDBE34",  # core
                               `FALSE` = "darkgrey")) # non-core
ggsave(paste0(getwd(),"/Figs/core_habitat.jpeg"), plot=last_plot(), 
       width=5.8, height=4.8, units="in")

# For each threshold, spatially plot the core area, with the WEAs in white
ggplot() +
  geom_raster(data = core_maps %>% filter(year %in% seq(2000, 2024, by = 6)),
              aes(x = x, y = y, fill = core)) +
  geom_sf(data = land, fill = "darkgrey") +
  geom_sf(data = WEA_poly, fill = NA, col = "white", lwd = 0.75) +
  facet_wrap(~threshold,
    ncol = 3,
    labeller = labeller(threshold = function(x) {paste0("exp(", sub("^core_", "", x), ")")})) +
  coord_sf(xlim = c(125, 950), ylim = c(4670, 5230)) +
  labs(title = expression(paste("Core area by biomass density threshold (kg/", km^2, ")")),
       y = "", x = "", fill = "") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), legend.position = "none",
    legend.key.width = unit(2, "cm"), legend.key.height = unit(0.4, "cm"),
    text = element_text(size = font_size), 
    strip.text = element_text(size = font_size * 0.8)) +
  scale_fill_manual(values = c(`TRUE` = "#FDBE34",`FALSE` = "darkgrey"))
ggsave(paste0(getwd(), "/Figs/core_habitat_2.jpeg"),
       plot = last_plot(), width = 6.3, height = 2.5, units = "in")

# For each threshold, plot the total core area over time
core_maps %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  filter(core == TRUE) %>%
  group_by(year, threshold) %>%
  summarise(n_cells = sum(core, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(thresh_num = as.numeric(substr(threshold, 6, 8)),
         threshold = factor(paste0("exp(", thresh_num, ")"),
                            levels = paste0("exp(", sort(unique(thresh_num)), ")"))) %>%
  ggplot(aes(x = year, y = n_cells*area_cell/1000, col = threshold)) + geom_path() +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, NA), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Total core area (1000s of square km)", 
       col = expression(paste("Threshold (kg/",km^{2},")"))) +
  theme(text = element_text(size = font_size), 
        panel.spacing.x = unit(0.6, "cm"),
        legend.text = element_text(size = font_size))
ggsave(paste0(getwd(),"/Figs/core_area_ts.jpeg"), plot=last_plot(), 
       width=8.4, height=4.8, units="in")

# Move forward with a biomass density threshold of exp(7)
core_maps_7 <- core_maps %>% filter(threshold == "core_7")

# Plot predictions across space for the core area
ggplot() + 
  geom_raster(data = core_maps_7[core_maps_7$year%in%seq(2000,2024,by = 6),], 
              aes(x = x, y = y), fill = "grey", col = NA) + 
  geom_raster(data = core_maps_7[core_maps_7$year%in%seq(2000,2024,by = 6)&
                                   core_maps_7$core==TRUE,], 
              aes(x = x, y = y, fill = est)) +
  ggtitle(expression(paste("Predicted relative biomass density (kg/",km^{2},")"))) +
  geom_sf(data = land, fill = "darkgrey") + 
  geom_sf(data = WEA_poly, fill = NA, col = "white", lwd = 0.75) + 
  theme(text = element_text(size = font_size), 
        plot.title = element_text(size = font_size), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(0.4, "cm")) +
  coord_sf(xlim = c(125,950), ylim = c(4670, 5230)) + 
  facet_wrap(~year) +
  scale_fill_gradientn("", colors = viridis(100), na.value = "grey",
                       trans = trans_new("log10", transform = function(x) log(x+1), 
                                         inverse = function(x) exp(x)-1), 
                       breaks = c(25, 500, 10000, NA), limits = c(25, NA))
ggsave(paste0(getwd(),"/Figs/spatial_predictions_core.jpeg"), plot=last_plot(), 
       width=6, height=4.4, units="in")

####Plot the percentage of total core biomass in the reserves and in each WEA###

# Merge core data with the raster data, and filter this data to the WEAs and reserves
core_raster_df <- raster_df %>% left_join(core_maps_7)
core_raster_df_WEA <- core_raster_df %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  st_intersection(st_make_valid(WEA_poly))
core_raster_df_reserves <- core_raster_df %>%
  st_as_sf(coords= c("x","y"), crs = 32620, remove = FALSE) %>%
  st_intersection(st_make_valid(reserves))

# Loop over WEAs to generate an index for each (core area only)
index_list2 <- list()
for (area in unique(core_raster_df_WEA$WEA)) {
      # Subset to WEA
      df_subset <- core_raster_df_WEA %>% 
        filter(WEA == area & core == TRUE)
      
      # If no area in core, index is 0
      if(nrow(df_subset)==0)
      {
        idx <- expand_grid(year = unique(raster_df$year), est = 0, lwr = 0,
                           upr = 0, log_est = 0, se = 0, se_natural = 0, 
                           type = "index",
                           WEA = area)
      }else{
        # Compute index and fill in missing years
        idx <- get_index(predict(fit_seacuke, df_subset[,,drop = T], 
                                 return_tmb_object = TRUE),
                         area = area_cell, bias_correct = TRUE) %>%
          tidyr::complete(year = unique(raster_df$year),
                          fill = list(est = 0, lwr = 0, upr = 0, 
                                      log_est = 0, se = 0, 
                                      se_natural = 0, type = "index"))
        idx$WEA <- area
      }
      
      # Store in list
      index_list2[[area]] <- idx
}
WEA_core_index_df <- do.call(rbind, index_list2)

# Get an overall core area index
combined_core_index_df <- get_index(predict(fit_seacuke, 
                                            core_raster_df[core_raster_df$core==TRUE,,
                                                           drop = T], 
                                            return_tmb_object = TRUE), area = area_cell, 
                                    bias_correct = TRUE)

# Get an overall index for the core area within the reserves
reserve_core_index_df <- get_index(predict(fit_seacuke, 
                                           core_raster_df_reserves[core_raster_df_reserves$core==TRUE,,
                                                                   drop = T], , 
                                            return_tmb_object = TRUE), area = area_cell, 
                                    bias_correct = TRUE)

# Calculate the percentage of core biomass in each WEA, and then plot using facet_wrap
WEA_core_index_df %>% left_join(combined_core_index_df, by = c("year")) %>% 
  mutate(prop = est.x/est.y) %>%
  ggplot(aes(x = year, y = prop*100)) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2026), ylim = c(0, 16), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Percentage of core area biomass", fill = "", col = "") +
  theme(text = element_text(size = font_size), panel.spacing.x = unit(0.6, "cm")) + 
  facet_wrap(~WEA)
ggsave(paste0(getwd(),"/Figs/core_biomass_proportion_WEAs.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

# Calculate the percentage of core biomass in the reserves, and then plot
reserve_core_index_df %>% left_join(combined_core_index_df, by = c("year")) %>% 
  mutate(prop = est.x/est.y) %>%
  ggplot(aes(x = year, y = prop*100)) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2026), ylim = c(0, 40), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Percentage of core area biomass within reserves", 
       fill = "", col = "") +
  theme_classic() + geom_hline(yintercept = 30, col = "red", lwd = 2, 
                               type = "dashed") +
  theme(text = element_text(size = font_size), panel.spacing.x = unit(0.6, "cm"))
ggsave(paste0(getwd(),"/Figs/core_biomass_proportion_reserves.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")

####Plot an index for the total core area###

# Plot the index
ggplot(combined_core_index_df, aes(x = year, y = est/1e6, 
                                   ymin = lwr/1e6, ymax = upr/1e6)) +
  geom_ribbon(alpha = 0.2, col = NA) +
  geom_path(size = 1.2) +
  coord_cartesian(xlim = c(2000, 2024), ylim = c(0, NA), expand = FALSE) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 10)) +
  labs(x = "Year", y = "Biomass index (kt)") +
  theme_classic() +
  theme(text = element_text(size = font_size), 
        legend.text = element_text(size = font_size),
        legend.position = "bottom")
ggsave(paste0(getwd(),"/Figs/index_core_area.jpeg"), plot=last_plot(), 
       width=7.2, height=5.5, units="in")