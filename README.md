Methodology
The analysis follows a four-step workflow:
Step 1: Boundary-Constrained Analysis
Initial boundary definition and constraint application for the study area.
Step 2: Width Calculation
Corridor width calculation based on spatial parameters.
Step 3: Fan-Shape Network Connectivity Generation

Script: step03_Fan-Shape Network Connectivity Generation.R
Creates fan-shaped network connections between nodes
Generates pair_node_fan_shape.gpkg output

Step 4: Least Cost Corridor (LCC) Analysis

Main Script: Calculates least-cost paths for wind corridors
Statistical Analysis: step04_ststic.R - Post-processing and statistics
Produces final wind corridor heatmap

Data Files
Input Data

area.gpkg - Study area polygon
stations.gpkg - Meteorological or sampling station locations
node_caddidate.gpkg - Potential corridor nodes
anisotropic_domains_final.gpkg - Directional flow domains
wind.rds - Wind data in R serialized format

Output Data

pair_node_fan_shape.gpkg - Network connectivity results
Final_WindCorridor_Heatmap.tif - Final corridor heatmap (GeoTIFF)
