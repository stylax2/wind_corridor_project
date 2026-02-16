################################################################################
# STEP 03: Fan-Shape Network Connectivity Generation
# Description: Generates a potential connectivity network by linking nodes 
#              in a 'fan-shape' pattern (connecting strictly to opposite sectors).
#              Optimized using vectorized operations for performance.
#
# Input:       Candidate Nodes (GPKG)
# Output:      Pairwise Connection Lines (GPKG)
################################################################################

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  sf,          # Spatial vector data
  dplyr,       # Data manipulation
  tidyr,       # Data tidying
  here,        # File path management
  logger       # Logging
)

# Initialize Logger
log_formatter(formatter_glue)
log_appender(appender_console)

# Parameters
INPUT_NODE_PATH  <- here("01data/node_cadidate.gpkg")
OUTPUT_LINE_PATH <- here("02Output/pair_node_fan_shape.gpkg")
OUTPUT_NODE_PATH <- here("02Output/final_node_points.gpkg")

TARGET_RANGE   <- 2    # Range of sectors to connect (Opposite +/- 2)
SECTOR_COUNT   <- 16   # Total number of sectors (16-wind directions)
DEFAULT_WEIGHT <- 1.0  # Default edge weight

# 2. Data Loading and Preprocessing --------------------------------------------
log_info("Step 3-1: Loading nodes and calculating spatial sectors...")

if (!file.exists(INPUT_NODE_PATH)) stop("Error: Input node file not found.")

# Load Nodes
nodes <- st_read(INPUT_NODE_PATH, quiet = TRUE)

# Ensure ID column exists and is unique
if (!"id" %in% names(nodes)) {
  nodes$id <- 1:nrow(nodes)
}
nodes <- nodes %>% mutate(id = as.integer(id)) # Ensure integer ID for speed

# Calculate Centroid of the Node Cluster
# (Used as the reference point for sector definition)
center_point <- st_union(nodes) %>% st_centroid() %>% st_coordinates()
cx <- center_point[1]
cy <- center_point[2]

# Extract Coordinates
coords <- st_coordinates(nodes)
nodes$x <- coords[, 1]
nodes$y <- coords[, 2]

# 3. Sector Assignment (Vectorized) --------------------------------------------
# Calculate angle relative to center (Azimuth)
dx <- nodes$x - cx
dy <- nodes$y - cy

# Atan2 returns radians (-pi to pi). Convert to degrees.
deg_math <- atan2(dy, dx) * (180 / pi)

# Convert Math Angle (CCW from East) to Azimuth (CW from North)
# Formula: Azimuth = (450 - Math_Angle) %% 360
azimuth <- (450 - deg_math) %% 360

# Assign Sector (1 to 16)
# Sector 1 centers on North (0/360 deg). Each sector is 360/16 = 22.5 deg wide.
nodes$sector <- floor((azimuth + 11.25) / 22.5) + 1
nodes$sector[nodes$sector > 16] <- 1

log_info("Sector assignment complete. Total nodes: {nrow(nodes)}")

# 4. Generate Connection Pairs (Vectorized) ------------------------------------
log_info("Step 3-2: Generating fan-shape connections (Target Range: +/- {TARGET_RANGE})...")

# Create Sector-to-Sector Connection Rule Table
sector_rules <- lapply(1:SECTOR_COUNT, function(src) {
  # Define Opposite Sector Index
  opposite <- (src + (SECTOR_COUNT / 2) - 1) %% SECTOR_COUNT + 1
  
  # Define Target Range (Opposite +/- Range)
  targets <- (opposite - 1 + (-TARGET_RANGE:TARGET_RANGE)) %% SECTOR_COUNT + 1
  
  data.frame(src_sector = src, target_sector = targets)
}) %>% bind_rows()

# Prepare Node Dataframes for Joining (Drop geometry for speed)
nodes_df <- nodes %>% st_drop_geometry()

# Join to create potential links:
links_df <- nodes_df %>%
  select(id, sector, x, y) %>%
  rename(src_id = id, src_sec = sector, src_x = x, src_y = y) %>%
  
  # 1st Join: Bind Source Nodes with Sector Rules
  inner_join(sector_rules, 
             by = c("src_sec" = "src_sector"), 
             relationship = "many-to-many") %>%
  
  # 2nd Join: Bind Target Nodes based on 'target_sector'
  # Note: The 'snk_sec' column is absorbed into 'target_sector' during join
  inner_join(nodes_df %>% 
               select(id, sector, x, y) %>%
               rename(snk_id = id, snk_sec = sector, snk_x = x, snk_y = y),
             by = c("target_sector" = "snk_sec"), 
             relationship = "many-to-many") %>%
  
  # Rename 'target_sector' back to 'snk_sec' for clarity
  rename(snk_sec = target_sector)

# Filter Constraints
# 1. Remove Self-loops (src != snk)
# 2. Remove Duplicates in Undirected Graph (keep only src_id < snk_id)
final_links <- links_df %>%
  filter(src_id < snk_id) %>%
  mutate(
    pair_id = paste0("S", src_sec, "_", src_id, "_to_S", snk_sec, "_", snk_id),
    type = "Fan_Shape",
    weight = DEFAULT_WEIGHT
  )

log_info("Connection generation complete. Total links: {nrow(final_links)}")

# 5. Geometry Construction and Export ------------------------------------------
log_info("Step 3-3: Constructing geometries and saving output...")

if (nrow(final_links) == 0) stop("Error: No connections were generated.")

# Create Linestrings from Coordinates (Efficient Method)
# Create a matrix of coordinates: columns (x1, y1, x2, y2)
coord_matrix <- as.matrix(final_links[, c("src_x", "src_y", "snk_x", "snk_y")])

# Function to create LINESTRING sfc from coordinate matrix row
make_line <- function(coords) {
  st_linestring(matrix(coords, ncol = 2, byrow = TRUE))
}

# Convert to sf object
# Note: Using lapply on matrix rows is faster than sf headers for massive data
geometry_list <- lapply(1:nrow(coord_matrix), function(i) make_line(coord_matrix[i,]))
geometry_sfc  <- st_sfc(geometry_list, crs = st_crs(nodes))

final_sf <- st_sf(final_links, geometry = geometry_sfc)

# Export Line Strings
st_write(final_sf, OUTPUT_LINE_PATH, delete_dsn = TRUE, quiet = TRUE)
log_info("Saved Line Geometries: {OUTPUT_LINE_PATH}")

# Export Unique Nodes (Optional: for visualization verification)
# We already have 'nodes' sf object, just save it cleanly
st_write(nodes, OUTPUT_NODE_PATH, delete_dsn = TRUE, quiet = TRUE)
log_info("Saved Node Geometries: {OUTPUT_NODE_PATH}")

log_info("✅ STEP 03 Completed Successfully.")