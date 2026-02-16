################################################################################
# STEP 02: Wide-Area Aerodynamic Friction Surface (AFS) Generation
# Description: Generates a continuous friction surface extended to the 'outside'
#              domain to ensure connectivity with external ecological nodes.
#              Global normalization is applied to topographic variables.
#
# Logic:
#   1. Load 'outside' data (Broad-scale topography).
#   2. Normalize variables (Roughness, Slope, TPI) based on global min/max.
#   3. Calculate AFS without masking to preserve external continuity.
#
# Input:       Broad-scale DEM, Slope, CHM, TPI (tif)
# Output:      Wide-Area AFS (tif)
################################################################################

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, sf, dplyr, here, logger)

# Set logger format
log_formatter(formatter_glue)
log_appender(appender_console)

# 2. Load Broad-Scale Topographic Data -----------------------------------------
log_info("Step 2-1: Loading broad-scale (outside) topographic data...")

path_outside <- here("01data/outside")

# Check file existence
req_files <- c("dem.tif", "slope.tif", "chm.tif", "tpi.tif")
if (!all(file.exists(file.path(path_outside, req_files)))) {
  stop("Error: Missing required topographic files in '01data/outside'.")
}

# Load Rasters
dem_out   <- rast(file.path(path_outside, "dem.tif"))
slope_out <- rast(file.path(path_outside, "slope.tif"))
chm_out   <- rast(file.path(path_outside, "chm.tif"))
tpi_out   <- rast(file.path(path_outside, "tpi.tif"))

# 3. Global Variable Normalization ---------------------------------------------
log_info("Step 2-2: Normalizing physical variables across the wide domain...")

# [Roughness] Canopy-based Surface Roughness
# Assumption: Roughness is proportional to 10% of canopy height
roughness_raw <- chm_out * 0.1

# Min-Max Normalization (0 to 1)
r_minmax <- minmax(roughness_raw) # Returns 2x1 matrix
r_min <- r_minmax[1, ]
r_max <- r_minmax[2, ]

roughness_std <- (roughness_raw - r_min) / (r_max - r_min)

# [Obstacle] Slope Normalization
# Standardize slope (degrees) to 0-1 range (assuming max slope 90 deg)
slope_std <- slope_out / 90 

# [Modulation] Topographic Position Index (TPI) Scaling
# Scale TPI to range [-1, 1] to represent valleys (-1) to ridges (+1)
tpi_minmax <- minmax(tpi_out)
tpi_min <- tpi_minmax[1, ]
tpi_max <- tpi_minmax[2, ]

tpi_scaled <- 2 * (tpi_out - tpi_min) / (tpi_max - tpi_min) - 1

# 4. Aerodynamic Friction Surface (AFS) Calculation ----------------------------
log_info("Step 2-3: Synthesizing Wide-Area AFS...")

# Coefficients (Empirically derived or Literature-based)
k1 <- 2.5  # Weight for Vegetation Roughness (Exponential drag)
k2 <- 1.5  # Weight for Topographic Obstacles (Slope)
k3 <- 0.4  # Weight for Ridge Acceleration (TPI bonus)

# AFS Equation: 
# Friction increases exponentially with roughness and slope.
# Friction decreases at ridges (High TPI) due to wind speed-up effects.
afs_outside <- exp(k1 * roughness_std + k2 * slope_std) * (1 - k3 * tpi_scaled)

# Note: No masking is applied here to allow connectivity analysis 
#       originating from outside the immediate study area.

# 5. Export Results ------------------------------------------------------------
output_path <- here("02Output/AFS.tif")

log_info("Step 2-4: Exporting result to {output_path}")

writeRaster(afs_outside, output_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))



