################################################################################
# Script: Boundary-Constrained Anisotropic Domain Analysis
# Description: This script delineates atmospheric domains based on 
#              anisotropic cost distance. 
#              * Optimized for memory stability (Sequential Processing).
#
# Data Input:  Terrain Rasters (TIF), Study Area (GPKG), Wind Data (RDS/CSV)
# Output:      Vectorized Domains (GPKG)
################################################################################

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  terra,       # Spatiotemporal data analysis
  sf,          # Simple features for vector data
  dplyr,       # Data manipulation
  gdistance,   # Cost distance and connectivity analysis
  here,        # Robust file path management
  lubridate,   # Date handling
  utils        # For progress bar
)

# 2. Data Acquisition and Pre-processing ---------------------------------------
message("Loading spatial data and defining study area...")

# Load study area boundary
study_area <- st_read(here("01data/area.gpkg"), quiet = TRUE)
target_crs <- st_crs(study_area)

# Load topographic rasters
path_topo  <- here("01data/inside_150")
dem        <- rast(file.path(path_topo, "dem.tif"))
slope      <- rast(file.path(path_topo, "slope.tif"))
chm        <- rast(file.path(path_topo, "chm.tif"))

# [Optional] Downsampling for Test (메모리가 부족하면 주석 해제하여 사용)
# message("Downsampling rasters for memory efficiency...")
# dem <- terra::aggregate(dem, fact = 2, fun = mean)
# slope <- terra::aggregate(slope, fact = 2, fun = mean)
# chm <- terra::aggregate(chm, fact = 2, fun = mean)

# 3. Surface Resistance Modeling -----------------------------------------------
message("Generating aerodynamic resistance surface...")

# Calculate aerodynamic resistance surface based on CHM and Slope
# Formula: exp((CHM * 0.1) + (Slope / 45))
base_res_terra <- exp((chm * 0.1) + (slope / 45))

# Constrain resistance surface to study area (Masking)
base_res_terra <- terra::crop(base_res_terra, vect(study_area))
base_res_terra <- terra::mask(base_res_terra, vect(study_area))

# Convert to 'raster' object for compatibility with gdistance package
base_raster <- raster::raster(base_res_terra)

# 4. Station Filtering and Alignment -------------------------------------------
message("Processing raw wind data and calculating prevailing wind directions...")

# 4-1. Load Raw Data (Daily Time-series)
wind_raw <- readRDS(here("01data/wind.rds"))

# 4-2. Quality Control: Filter stations with >80% data availability
# Calculate total span (wind_raw$Date is already Date format)
total_days_possible <- as.numeric(max(wind_raw$Date, na.rm = TRUE) - min(wind_raw$Date, na.rm = TRUE))
if(total_days_possible <= 0) total_days_possible <- 1 

data_threshold <- total_days_possible * 0.7

valid_sites <- wind_raw %>%
  group_by(site) %>%
  summarise(n_records = n()) %>%
  filter(n_records >= data_threshold) %>%
  pull(site)

message(paste("Filtered stations: Retaining", length(valid_sites), 
              "sites with >80% data completeness."))

# 4-3. Calculate Wind Regime (Vector Averaging)
wind_regime <- wind_raw %>%
  filter(site %in% valid_sites) %>%
  mutate(
    # Convert Degrees to Radians
    wd_rad = wd * pi / 180,
    
    # Decompose into U (East-West) and V (North-South) components
    u = -ws * sin(wd_rad),
    v = -ws * cos(wd_rad)
  ) %>%
  group_by(site) %>%
  summarise(
    mean_u   = mean(u, na.rm = TRUE),
    mean_v   = mean(v, na.rm = TRUE),
    mean_ws  = mean(ws, na.rm = TRUE),
    resultant_ws = sqrt(mean(u, na.rm = TRUE)^2 + mean(v, na.rm = TRUE)^2),
    # Prevailing Direction (Azimuth 0-360)
    prevailing = (atan2(mean(u, na.rm=TRUE), mean(v, na.rm=TRUE)) * 180 / pi + 180) %% 360,
    .groups  = "drop"
  ) %>%
  mutate(
    site = as.character(site),
    constancy = resultant_ws / mean_ws
  )

# 4-4. Spatial Join with Station Geometry
stations <- st_read(here("01data/stations.gpkg"), quiet = TRUE) %>%
  mutate(site = as.character(site)) %>%
  inner_join(wind_regime, by = "site") %>%
  st_transform(target_crs) %>%
  st_intersection(study_area)

message(paste("Final Station Count aligned with Study Area:", nrow(stations)))

# 5. Static Transition Matrix Calculation --------------------------------------
message("Calculating static transition matrix and geocorrection...")

# Create transition layer (Inverse of resistance = Conductance)
tr <- gdistance::transition(1/base_raster, transitionFunction = mean, directions = 8)

# Apply geocorrection for anisotropic distance
tr_corr <- gdistance::geoCorrection(tr, type = "c")

# 6. Sequential Cost Distance Computation (Memory Safe) ------------------------
message(paste("Initiating sequential cost calculation for", nrow(stations), "stations..."))

# Storage for results
cost_layers_list <- list()
n_stations <- nrow(stations)

# Progress bar setup
pb <- txtProgressBar(min = 0, max = n_stations, style = 3)

for(i in 1:n_stations) {
  
  # Station geometry extraction
  st_node <- stations[i, ]
  coords  <- st_coordinates(st_node)
  
  # Calculate cumulative cost surface
  # accCost returns a RasterLayer (gdistance/raster)
  cost_raster_layer <- gdistance::accCost(tr_corr, coords)
  
  # Convert to terra::rast and store
  cost_layers_list[[i]] <- terra::rast(cost_raster_layer)
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # [Critical] Memory Management: Garbage Collection every 10 iterations
  if(i %% 10 == 0) gc()
}

close(pb)

# 7. Domain Integration and Vectorization --------------------------------------
message("\nIntegrating domains and generating polygons...")

# Stack all cost layers
cost_stack <- terra::rast(cost_layers_list)

# Identify the index of the layer with the minimum cost for each pixel
domain_id_raster <- terra::which.min(cost_stack)

# Convert raster domains to polygons
domain_sf <- terra::as.polygons(domain_id_raster, aggregate = TRUE) %>%
  st_as_sf() %>%
  st_set_crs(5179) %>%                      # Explicitly define projected CRS
  rename(id = 1) %>%
  mutate(id = as.integer(id)) %>%
  inner_join(                               # Attach station attributes
    stations %>% 
      mutate(id = row_number()) %>% 
      st_drop_geometry(), 
    by = "id"
  ) %>%
  st_transform(target_crs) %>%              # Ensure CRS consistency
  st_intersection(study_area)               # Final clipping

# 8. Output Export -------------------------------------------------------------
output_path <- here("02Output/anisotropic_domains_final.gpkg")
st_write(domain_sf, output_path, delete_dsn = TRUE, quiet = TRUE)
write.csv(wind_regime, "02Output/wind_regime.csv")

message("Analysis Complete. Output saved to: ", output_path)

# Cleanup
rm(tr, tr_corr, cost_layers_list, cost_stack)
gc()
