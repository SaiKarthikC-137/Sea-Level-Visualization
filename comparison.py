import xarray as xr
import rasterio
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os # For checking file existence for caching

# --- Configuration ---
YEAR1 = 2023
YEAR2 = 2024

# IMPORTANT: Adjust this pattern to match your RAW data file names and paths.
# This is for the original NetCDF files for each year.
RAW_DATA_FILE_PATTERN = "gebco_{year}.nc" # Example: "gebco_2020.nc", "gebco_2024.nc"
# If your files are like "GEBCO_YYYY_sub_ice_topo.nc", use:
# RAW_DATA_FILE_PATTERN = "GEBCO_{year}_sub_ice_topo.nc"


# Output filename for the difference GeoTIFF (using the diverging colormap)
# The name reflects that it uses user's specific transform/flip method
OUTPUT_FILENAME_PATTERN = "ocean_elevation_{year2}_vs_{year1}.tif"

# Directory for intermediate cached files
CACHE_DIR = "intermediate_cache"
os.makedirs(CACHE_DIR, exist_ok=True) # Ensure cache directory exists

# Dask chunk size for initial loading
CHUNK_SIZE = {'lat': 500, 'lon': 500}

# Coarsening factor
COARSEN_FACTOR = 10

# --- Helper function to load and preprocess a single year's data with caching ---
def load_and_preprocess_year(year, raw_file_pattern, cache_dir, chunk_size, coarsen_factor):
    """
    Loads, filters, and downsamples ocean elevation data for a given year.
    Implements caching for the processed 'ocean_coarse' DataArray.
    Latitude order is preserved as is from the coarsened data.
    """
    intermediate_filepath = os.path.join(cache_dir, f"ocean_coarse_{year}_intermediate.nc")

    if os.path.exists(intermediate_filepath):
        print(f"Loading pre-computed 'ocean_coarse' data for year {year} from {intermediate_filepath}...")
        ocean_coarse = xr.open_dataarray(intermediate_filepath)
        print(f"'ocean_coarse' for year {year} loaded successfully from cache.")
    else:
        print(f"No pre-computed data for year {year} found at {intermediate_filepath}. Computing from scratch...")
        raw_filepath = raw_file_pattern.format(year=year)
        print(f"Loading raw data for year {year} from {raw_filepath}")
        try:
            ds = xr.open_dataset(raw_filepath, chunks=chunk_size)
        except FileNotFoundError:
            print(f"ERROR: Raw data file not found: {raw_filepath}")
            raise
        except Exception as e:
            print(f"ERROR: Could not open raw dataset {raw_filepath}: {e}")
            raise

        if 'elevation' not in ds.data_vars:
            raise ValueError(f"Variable 'elevation' not found in {raw_filepath}. Available: {list(ds.data_vars)}")

        elevation = ds['elevation']
        ocean_elevation = elevation.where(elevation < 0) # Filter for ocean
        ocean_coarse = ocean_elevation.coarsen(lat=coarsen_factor, lon=coarsen_factor, boundary="trim").mean()

        try:
            ocean_coarse.to_netcdf(intermediate_filepath)
            print(f"Saved 'ocean_coarse' data for year {year} to {intermediate_filepath}")
        except Exception as e:
            print(f"Error saving 'ocean_coarse' data for year {year} to cache: {e}")
        print(f"Initial processing for year {year} complete. Coarse shape: {ocean_coarse.shape}")

    # We do NOT force latitude sort here, to allow user's flip logic to operate on original order.
    return ocean_coarse.load()

# --- Main script ---
output_filename = OUTPUT_FILENAME_PATTERN.format(year1=YEAR1, year2=YEAR2)

print(f"--- Preparing data for Year 1: {YEAR1} ---")
try:
    ocean_coarse_y1 = load_and_preprocess_year(YEAR1, RAW_DATA_FILE_PATTERN, CACHE_DIR, CHUNK_SIZE, COARSEN_FACTOR)
except Exception as e:
    print(f"Failed to process data for {YEAR1}. Exiting. Error: {e}")
    exit(1)

print(f"\n--- Preparing data for Year 2: {YEAR2} ---")
try:
    ocean_coarse_y2 = load_and_preprocess_year(YEAR2, RAW_DATA_FILE_PATTERN, CACHE_DIR, CHUNK_SIZE, COARSEN_FACTOR)
except Exception as e:
    print(f"Failed to process data for {YEAR2}. Exiting. Error: {e}")
    exit(1)

print("\n--- Aligning datasets ---")
try:
    ocean_coarse_y1_aligned, ocean_coarse_y2_aligned = xr.align(ocean_coarse_y1, ocean_coarse_y2, join="inner")
    print("Successfully aligned datasets using 'inner' join.")
    print(f"Shape after alignment: {ocean_coarse_y1_aligned.shape}")
    if ocean_coarse_y1_aligned.size == 0:
        print("Warning: Alignment resulted in an empty dataset. Check input data overlap and grid consistency.")
        exit(1)
except Exception as e:
    print(f"Could not automatically align datasets: {e}")
    if not (np.array_equal(ocean_coarse_y1.lat.values, ocean_coarse_y2.lat.values) and \
            np.array_equal(ocean_coarse_y1.lon.values, ocean_coarse_y2.lon.values)):
        print("CRITICAL WARNING: Coordinates of coarsened datasets do not match, and alignment failed/skipped.")
        exit(1)
    ocean_coarse_y1_aligned = ocean_coarse_y1
    ocean_coarse_y2_aligned = ocean_coarse_y2


print("\n--- Calculating elevation difference ---")
ocean_elevation_diff = ocean_coarse_y2_aligned - ocean_coarse_y1_aligned
print("Difference calculation complete.")
print(f"Shape of difference data: {ocean_elevation_diff.shape}")

# --- Exploratory Data Analysis (EDA) of the Difference Data ---
print("\n--- Analyzing the Difference Data ---")
diff_values_np = ocean_elevation_diff.values # This is the numerical data for the difference

actual_min_diff = np.nan
actual_max_diff = np.nan

if diff_values_np.size == 0:
    print("Difference data array is empty. No statistics to report.")
else:
    valid_diff_values = diff_values_np[~np.isnan(diff_values_np)]
    print(f"Total number of grid cells in difference data: {diff_values_np.size}")
    if valid_diff_values.size > 0:
        actual_min_diff = np.nanmin(diff_values_np)
        actual_max_diff = np.nanmax(diff_values_np)
    print(f"Min difference: {actual_min_diff:.4f}")
    print(f"Max difference: {actual_max_diff:.4f}")
    # Add more EDA stats here if desired (mean, median, unique counts, percentiles etc.)

# --- Define Diverging Colormap for Differences ---
print(f"\n--- Preparing Diverging Colormap ('RdYlBu_r') centered at Zero ---")
if np.isnan(actual_min_diff) or np.isnan(actual_max_diff):
    sym_max_abs = 1.0
else:
    sym_max_abs = np.nanmax([np.abs(actual_min_diff), np.abs(actual_max_diff)])
if sym_max_abs == 0: sym_max_abs = 0.1

vmin_norm_diverging = -sym_max_abs
vmax_norm_diverging = sym_max_abs
cmap_diverging = plt.get_cmap('RdYlBu_r') # Blue -> Yellow (0) -> Red
norm_diverging = colors.Normalize(vmin=vmin_norm_diverging, vmax=vmax_norm_diverging)
print(f"Colormap limits for diverging map: vmin={vmin_norm_diverging:.2f}, vmax={vmax_norm_diverging:.2f}.")

# --- Apply colormap ---
print("\n--- Applying diverging colormap to difference data ---")
# rgba_array_diff is equivalent to 'rgba_array' in the single-year script
rgba_array_diff = cmap_diverging(norm_diverging(diff_values_np))
# image_data_for_geotiff is equivalent to 'image_data' in the single-year script
image_data_for_geotiff = (rgba_array_diff * 255).astype(np.uint8)


# --- Prepare metadata and write to GeoTIFF (User's Exact Method) ---
print("\n--- Preparing GeoTIFF metadata (User's method) for difference map ---")

lon_coords = ocean_elevation_diff['lon'].values
lat_coords = ocean_elevation_diff['lat'].values

# This will be the final data array to write, possibly flipped.
# Equivalent to 'image_data_to_write' in the single-year script.
image_data_to_write_final = image_data_for_geotiff

# Critical check for valid dimensions (from user's script)
if len(lon_coords) <= 1 or len(lat_coords) <= 1 or image_data_for_geotiff.shape[1] == 0 or image_data_for_geotiff.shape[0] == 0:
    raise ValueError(f"Insufficient dimensions after processing to define transform or write data. "
                     f"Image width: {image_data_for_geotiff.shape[1]}, Image height: {image_data_for_geotiff.shape[0]}, "
                     f"Num lon: {len(lon_coords)}, Num lat: {len(lat_coords)}")

# Conditional flip based on latitude coordinate order (from user's script)
if len(lat_coords) > 1 and lat_coords[0] < lat_coords[-1]: # Check if latitudes are ascending
    print("Latitudes are ascending. Flipping image data vertically.")
    image_data_to_write_final = np.flipud(image_data_for_geotiff)
    # Note: If latitudes were ascending, lat_coords[0] is South, lat_coords[-1] is North.
    # lat_coords.min() is South, lat_coords.max() is North.
else:
    print("Latitudes are descending or single. No vertical flip of image data needed based on this condition.")
    # If latitudes are descending, lat_coords[0] is North, lat_coords[-1] is South.
    # lat_coords.min() is South, lat_coords.max() is North.

# Calculate pixel dimensions (from user's script)
calculated_pixel_width = (lon_coords[-1] - lon_coords[0]) / (len(lon_coords) - 1)
positive_pixel_height = abs((lat_coords[0] - lat_coords[-1]) / (len(lat_coords) - 1)) # Ensures positive

# Coordinates for from_origin are the outer upper-left corner (from user's script)
ul_lon_corner = lon_coords.min() - calculated_pixel_width / 2.0
ul_lat_corner = lat_coords.max() + positive_pixel_height / 2.0 # lat_coords.max() is the northernmost latitude value

# Create transform using positive pixel height as y-resolution (from user's script)
users_method_transform = from_origin(ul_lon_corner,
                                     ul_lat_corner,
                                     calculated_pixel_width,
                                     positive_pixel_height) # User's method uses POSITIVE y-resolution

out_meta = {
    "driver": "GTiff",
    "height": image_data_to_write_final.shape[0],
    "width": image_data_to_write_final.shape[1],
    "count": 4,  # RGBA
    "dtype": 'uint8',
    "crs": "EPSG:4326",
    "transform": users_method_transform # Use the user's method transform
}

print(f"Writing difference GeoTIFF to {output_filename} using user's transform logic...")
try:
    with rasterio.open(output_filename, "w", **out_meta) as dest:
        # Handle NaNs in the source data by making them transparent in the RGBA image
        # This is done when creating rgba_array_diff if NaNs propagate through norm,
        # or can be explicitly set: image_data_to_write_final[np.isnan(diff_values_np)] = [0,0,0,0]
        # However, our rgba_array_diff already handles NaNs from diff_values_np by colormapper.
        # The image_data_to_write_final is uint8, so direct NaN check isn't possible here.
        # The transparency was baked into rgba_array_diff which then became image_data_for_geotiff.

        if image_data_to_write_final.ndim == 3 and image_data_to_write_final.shape[2] == out_meta["count"]:
            for i in range(out_meta["count"]):
                dest.write(image_data_to_write_final[:, :, i], i + 1)
        else:
            raise ValueError(f"image_data_to_write_final dimensions ({image_data_to_write_final.shape}) or count ({out_meta['count']}) mismatch.")
    print(f"{output_filename} has been written successfully (using user's flip/transform logic).")
except Exception as e:
    print(f"Error writing GeoTIFF {output_filename}: {e}")

print("\n--- Script finished ---")