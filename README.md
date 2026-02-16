# Wind Corridor Project

## Overview
This project identifies and analyzes wind corridors using spatial network analysis and least-cost path algorithms. The analysis integrates geographic data, anisotropic domains, and fan-shaped network connectivity to determine optimal wind flow pathways.

## Project Structure

```
wind_corridor_project/
├── Final_WindCorridor_Heatmap.tif    # Final output heatmap
├── LICENSE                            # Project license
├── README.md                          # This file
├── anisotropic_domains_final.gpkg    # Anisotropic domain geometries
├── area.gpkg                         # Study area boundary
├── node_caddidate.gpkg               # Candidate node locations
├── pair_node_fan_shape.gpkg          # Fan-shaped network pairs
├── stations.gpkg                     # Station locations
├── step01_Boundary-Constrained A...  # Step 1 script
├── step02_Wid...                     # Step 2 script
├── step03_Fan-Shape Network Conn...  # Step 3 script (R)
├── step04_Least Cost Corridor (LCC)...  # Step 4 script
├── step04_ststic.R                   # Step 4 statistical analysis
└── wind.rds                          # Wind data (R object)
```

## Methodology

The analysis follows a four-step workflow:

### Step 1: Boundary-Constrained Analysis
Initial boundary definition and constraint application for the study area.

### Step 2: Width Calculation
Corridor width calculation based on spatial parameters.

### Step 3: Fan-Shape Network Connectivity Generation
- **Script**: `step03_Fan-Shape Network Connectivity Generation.R`
- Creates fan-shaped network connections between nodes
- Generates `pair_node_fan_shape.gpkg` output

### Step 4: Least Cost Corridor (LCC) Analysis
- **Main Script**: Calculates least-cost paths for wind corridors
- **Statistical Analysis**: `step04_ststic.R` - Post-processing and statistics
- Produces final wind corridor heatmap

## Data Files

### Input Data
- `area.gpkg` - Study area polygon
- `stations.gpkg` - Meteorological or sampling station locations
- `node_caddidate.gpkg` - Potential corridor nodes
- `anisotropic_domains_final.gpkg` - Directional flow domains
- `wind.rds` - Wind data in R serialized format

### Output Data
- `pair_node_fan_shape.gpkg` - Network connectivity results
- `Final_WindCorridor_Heatmap.tif` - Final corridor heatmap (GeoTIFF)

## Requirements

### R Packages
```r
# Spatial analysis
library(sf)
library(terra)
library(stars)

# Network analysis
library(igraph)
library(gdistance)

# Data manipulation
library(dplyr)
library(tidyr)

# Visualization
library(ggplot2)
library(tmap)
```

### Python Dependencies (if applicable)
```python
# Add if Python scripts are used
import geopandas as gpd
import rasterio
import numpy as np
```

## Usage

### Running the Analysis

1. **Prepare input data**
   ```r
   # Load study area and stations
   area <- st_read("area.gpkg")
   stations <- st_read("stations.gpkg")
   wind_data <- readRDS("wind.rds")
   ```

2. **Execute analysis steps sequentially**
   ```bash
   # Run each step in order
   Rscript step01_Boundary-Constrained_A....R
   Rscript step02_Wid....R
   Rscript step03_Fan-Shape_Network_Connectivity_Generation.R
   Rscript step04_Least_Cost_Corridor_LCC....R
   Rscript step04_ststic.R
   ```

3. **View results**
   ```r
   # Load and visualize final heatmap
   library(terra)
   heatmap <- rast("Final_WindCorridor_Heatmap.tif")
   plot(heatmap)
   ```

## Key Concepts

### Anisotropic Domains
Directional flow regions where wind movement is influenced by terrain and obstacles, stored in `anisotropic_domains_final.gpkg`.

### Fan-Shape Network
A network topology that models wind dispersal patterns radiating from source points, accounting for directional preferences.

### Least Cost Corridor (LCC)
Optimal pathways for wind flow calculated using cost-distance algorithms, minimizing resistance across the landscape.

## Output Interpretation

The final heatmap (`Final_WindCorridor_Heatmap.tif`) represents:
- **High values**: Primary wind corridor zones
- **Low values**: Areas with restricted wind flow
- **Spatial resolution**: 100m
- **CRS**: EPSG 5179

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/analysis-improvement`)
3. Commit your changes (`git commit -am 'Add new analysis method'`)
4. Push to the branch (`git push origin feature/analysis-improvement`)
5. Create a Pull Request


