"""
Script: Least Cost Corridor (LCC)
Dependencies: rasterio, geopandas, shapely, skimage, scipy, numpy, tqdm
"""

import os
import shutil
import tempfile
import logging
import glob
import numpy as np
import rasterio
import geopandas as gpd
from shapely.geometry import LineString
from skimage.graph import MCP_Geometric
from scipy.spatial import cKDTree
from tqdm import tqdm
from typing import Tuple, List, Optional, Dict

# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

class EcologicalCorridorAnalyzer:
    """
    A geospatial analyzer for computing Least Cost Paths (LCP) and 
    Probabilistic Connectivity Corridors over a friction surface.
    """

    def __init__(self, friction_path: str, temp_dir: Optional[str] = None):
        """
        Initialize the analyzer, load friction data, and setup temporary storage.
        
        Args:
            friction_path (str): Path to the input friction raster (GeoTIFF).
            temp_dir (str, optional): Directory for temporary processing files.
                                      Recommended to use a fast SSD drive.
        """
        self.friction_path = friction_path
        self.acs_cache: Dict[str, str] = {}  # Cache for Accumulated Cost Surface paths

        # Setup Temporary Directory
        if temp_dir:
            os.makedirs(temp_dir, exist_ok=True)
            self.temp_dir_obj = tempfile.TemporaryDirectory(dir=temp_dir)
        else:
            self.temp_dir_obj = tempfile.TemporaryDirectory()
        
        self.temp_dir = self.temp_dir_obj.name
        logging.info(f"Workspace initialized at: {self.temp_dir}")

        # Load and Preprocess Friction Surface
        self._load_and_process_friction()

    def _load_and_process_friction(self):
        """
        Loads the friction raster, creates a validity mask, builds a KDTree for snapping,
        and saves the data as a memory-mapped file for RAM efficiency.
        """
        try:
            with rasterio.open(self.friction_path) as src:
                self.profile = src.profile.copy()
                self.transform = src.transform
                self.shape = src.shape
                self.nodata = src.nodata
                self.crs = src.crs
                
                logging.info("Step 1: Loading and preprocessing friction surface...")
                data = src.read(1).astype('float32')

                # Create Validity Mask (Exclude NoData and NaN)
                if self.nodata is not None:
                    valid_mask = (data > 0) & (data != self.nodata) & (~np.isnan(data))
                else:
                    valid_mask = (data > 0) & (~np.isnan(data))

                # Build KDTree for Coordinate Snapping (Nearest Neighbor)
                self.valid_indices = np.argwhere(valid_mask)
                self.tree = cKDTree(self.valid_indices)

                # Soft Barrier Implementation:
                # Assign a very high cost to invalid areas instead of strict NoData
                # to prevent graph disconnection errors in MCP.
                max_val = np.nanmax(data[valid_mask])
                self.fill_value = max_val * 100.0 
                data[~valid_mask] = self.fill_value

                # Save as Memmap to handle large rasters without OOM errors
                self.friction_memmap_path = os.path.join(self.temp_dir, 'friction.dat')
                friction_mm = np.memmap(self.friction_memmap_path, dtype='float32', mode='w+', shape=self.shape)
                friction_mm[:] = data[:]
                friction_mm.flush()
                
                logging.info(f"Friction surface loaded. Dimensions: {self.shape}")

        except Exception as e:
            logging.error(f"Failed to load friction surface: {e}")
            raise

    def _coord_to_idx(self, x: float, y: float) -> Tuple[int, int]:
        """Converts map coordinates (x, y) to raster indices (row, col)."""
        row, col = ~self.transform * (x, y)
        return int(row), int(col)

    def snap_to_valid(self, r: int, c: int) -> Tuple[int, int]:
        """
        Snaps a target pixel index to the nearest valid friction pixel.
        This ensures start/end points do not fall on barriers.
        """
        need_snap = False
        
        # Check bounds
        if not (0 <= r < self.shape[0] and 0 <= c < self.shape[1]):
            need_snap = True
        else:
            # Check value (Soft barrier detection)
            friction = np.memmap(self.friction_memmap_path, dtype='float32', mode='r', shape=self.shape)
            val = friction[r, c]
            if val >= (self.fill_value * 0.9): 
                need_snap = True
            del friction 

        if need_snap:
            # Query KDTree for nearest valid pixel
            dist, idx = self.tree.query([r, c], k=1)
            return tuple(self.valid_indices[idx])
        
        return (r, c)

    def calculate_acs(self, node_key: str, start_idx: Tuple[int, int]) -> str:
        """
        Calculates or retrieves the Accumulated Cost Surface (ACS) from a source node.
        Returns the file path to the memory-mapped ACS result.
        """
        # Return cached path if available
        if node_key in self.acs_cache:
            if os.path.exists(self.acs_cache[node_key]):
                return self.acs_cache[node_key]

        # Compute ACS using Minimum Cost Path algorithm
        friction = np.memmap(self.friction_memmap_path, dtype='float32', mode='r', shape=self.shape)
        mcp = MCP_Geometric(friction, fully_connected=True)
        cumulative_costs, _ = mcp.find_costs(starts=[start_idx])

        # Save result to memmap
        path = os.path.join(self.temp_dir, f'acs_{node_key}.dat')
        acs_mm = np.memmap(path, dtype='float32', mode='w+', shape=self.shape)
        acs_mm[:] = cumulative_costs[:]
        acs_mm.flush()

        # Update cache and cleanup
        self.acs_cache[node_key] = path
        del acs_mm, cumulative_costs, friction, mcp
        
        return path

    def run_analysis(self, gpkg_path: str, output_folder: str, lcc_percent: float = 0.05):
        """
        Executes the main analysis pipeline:
        1. Iterates through node pairs.
        2. Calculates Least Cost Corridors.
        3. Aggregates results into a regional heatmap.
        
        Args:
            gpkg_path (str): Path to the input GPKG file containing node pairs (lines).
            output_folder (str): Destination folder for results.
            lcc_percent (float): Threshold for corridor width definition (e.g., 0.05 = top 5%).
        """
        os.makedirs(output_folder, exist_ok=True)
        
        logging.info("Step 2: Loading vector pair data...")
        gdf = gpd.read_file(gpkg_path)
        if gdf.crs != self.crs:
            logging.warning("CRS mismatch detected. Reprojecting vector data to match raster...")
            gdf = gdf.to_crs(self.crs)

        lcp_results = []
        logging.info(f"Starting pairwise analysis for {len(gdf)} links. Threshold: {lcc_percent*100}%")

        for idx, row in tqdm(gdf.iterrows(), total=len(gdf), desc="Processing Corridors"):
            temp_resources = [] 
            
            try:
                pair_id = str(row.get('pair_id', idx))
                energy_weight = float(row.get('weight', 1.0))
                
                # Extract coordinates from LineString geometry
                coords = list(row.geometry.coords)
                src_x, src_y = coords[0]
                snk_x, snk_y = coords[-1]

                # Convert to raster indices and snap to valid pixels
                src_idx = self.snap_to_valid(*self._coord_to_idx(src_x, src_y))
                snk_idx = self.snap_to_valid(*self._coord_to_idx(snk_x, snk_y))
                
                # Generate unique keys for caching
                src_key = f"node_{src_idx[0]}_{src_idx[1]}"
                snk_key = f"node_{snk_idx[0]}_{snk_idx[1]}"

                # Load Accumulated Cost Surfaces (ACS)
                # Cost from Source -> All Pixels
                path_a = self.calculate_acs(src_key, src_idx)
                # Cost from Sink -> All Pixels
                path_b = self.calculate_acs(snk_key, snk_idx)

                # Map ACS files to memory
                acs_a = np.memmap(path_a, dtype='float32', mode='r', shape=self.shape)
                acs_b = np.memmap(path_b, dtype='float32', mode='r', shape=self.shape)
                temp_resources.extend([acs_a, acs_b])

                # Determine Minimum Cost (LCP Cost)
                min_cost = acs_a[snk_idx]
                
                if np.isinf(min_cost) or min_cost == 0:
                    continue

                # Calculate Corridor Surface: C(x) = Cost(Src, x) + Cost(x, Snk)
                corridor_surface = acs_a + acs_b
                
                # Define Threshold for Corridor Width
                threshold_val = min_cost * (1 + lcc_percent)
                
                # Generate Weight Map (Inverted Cost)
                # Normalizes values so LCP = 1 and Threshold Edge = 0
                denom = max(threshold_val - min_cost, 1e-7)
                weight_map = np.where(
                    corridor_surface <= threshold_val,
                    (threshold_val - corridor_surface) / denom, 
                    0
                )
                final_map = np.maximum(0, weight_map) * energy_weight

                # Save individual corridor to temporary TIFF
                temp_tif_path = os.path.join(output_folder, f"TEMP_corridor_{idx}.tif")
                self._save_temp_raster(temp_tif_path, final_map)

                # Extract Least Cost Path (LCP) Geometry (Traceback)
                friction = np.memmap(self.friction_memmap_path, dtype='float32', mode='r', shape=self.shape)
                mcp_trace = MCP_Geometric(friction, fully_connected=True)
                mcp_trace.find_costs(starts=[src_idx], ends=[snk_idx])
                lcp_indices = mcp_trace.traceback(snk_idx)
                
                # Convert indices back to coordinates
                lcp_coords = [self.transform * (c, r) for r, c in lcp_indices]

                if len(lcp_coords) > 1:
                    lcp_results.append({
                        'pair_id': pair_id,
                        'weight': energy_weight,
                        'cost': float(min_cost),
                        'geometry': LineString(lcp_coords)
                    })
                
                del friction, mcp_trace

            except Exception as e:
                logging.error(f"Error processing pair ID {pair_id}: {e}")
            
            finally:
                # Explicit cleanup to release file handles
                for res in temp_resources:
                    del res

        # Final Aggregation
        logging.info("Step 3: Aggregating individual corridors...")
        self._aggregate_rasters(output_folder)
        
        # Save LCP Vectors
        if lcp_results:
            out_gpkg = os.path.join(output_folder, "Final_LCP_Pathways.gpkg")
            gpd.GeoDataFrame(lcp_results, crs=self.crs).to_file(out_gpkg, driver="GPKG")
            logging.info(f"LCP vector data saved to: {out_gpkg}")

    def _save_temp_raster(self, path: str, data: np.ndarray):
        """Helper function to save temporary raster files with compression."""
        out_meta = self.profile.copy()
        out_meta.update(dtype='float32', nodata=0, compress='lzw')
        with rasterio.open(path, 'w', **out_meta) as dst:
            dst.write(data.astype('float32'), 1)

    def _aggregate_rasters(self, folder: str):
        """Sums all temporary corridor rasters to create the final connectivity heatmap."""
        search_pattern = os.path.join(folder, "TEMP_corridor_*.tif")
        files = glob.glob(search_pattern)
        
        if not files:
            logging.warning("No temporary raster files found for aggregation.")
            return

        # Initialize accumulator array
        aggregated_landscape = np.zeros(self.shape, dtype='float32')

        for f in tqdm(files, desc="Aggregating Heatmap"):
            try:
                with rasterio.open(f) as src:
                    aggregated_landscape += src.read(1)
                os.remove(f) # Delete temp file immediately after processing to save space
            except Exception as e:
                logging.warning(f"Failed to process/delete temp file {f}: {e}")

        # Save Final Output
        final_path = os.path.join(folder, "Final_WindCorridor_Heatmap.tif")
        self._save_temp_raster(final_path, aggregated_landscape)
        logging.info(f"Final Heatmap saved successfully: {final_path}")

    def cleanup(self):
        """Clean up temporary directory and resources."""
        try:
            self.temp_dir_obj.cleanup()
            logging.info("Temporary workspace cleaned up.")
        except Exception as e:
            logging.warning(f"Error cleaning up temp dir: {e}")

# =========================================================
# Main Execution Block
# =========================================================
if __name__ == "__main__":
    # Input Configuration
    # Ensure these paths match your previous R script outputs
    FRICTION_FILE = "AFS.tif" 
    PAIR_GPKG     = "pair_node_fan_shape.gpkg" 
    OUTPUT_DIR    = "Result_lcc"
    
    # Analysis Parameters
    # LCC_THRESHOLD: Defines the breadth of the corridor. 
    # 0.2 means the corridor includes paths up to 20% more costly than the optimal LCP.
    LCC_THRESHOLD = 0.2 

    # Initialize and Run
    analyzer = EcologicalCorridorAnalyzer(FRICTION_FILE, temp_dir="./temp_workspace_python")
    try:
        analyzer.run_analysis(PAIR_GPKG, OUTPUT_DIR, lcc_percent=LCC_THRESHOLD)
    finally:
        analyzer.cleanup()
