import xarray as xr
import rasterio
from rasterio.transform import from_origin
import matplotlib.colors as colors
import numpy as np
import cmocean
import os

# --- Configuration: SET YOUR INPUT FILE HERE ---
input_nc_filename = "gebco_2020.nc"  # Example: "gebco_2023.nc", "gebco_2022.nc", etc.

# --- Generate dynamic filenames based on the input ---
base_name_from_input = os.path.splitext(os.path.basename(input_nc_filename))[0]
intermediate_data_filepath = f"ocean_coarse_{base_name_from_input}_intermediate.nc"
output_tif_filename = f"styled_output_{base_name_from_input}.tif"
# You might also want a dynamic name for the gdal2tiles output folder later
# output_tiles_directory = f"ocean_tiles_{base_name_from_input}"

print(f"--- Processing for input: {input_nc_filename} ---")
print(f"Intermediate cache file: {intermediate_data_filepath}")
print(f"Output GeoTIFF file: {output_tif_filename}")

# --- Try to load intermediate data ---
if os.path.exists(intermediate_data_filepath):
    print(f"Loading pre-computed 'ocean_coarse' data from {intermediate_data_filepath}...")
    ocean_coarse = xr.open_dataarray(intermediate_data_filepath)
    print("'ocean_coarse' loaded successfully.")
else:
    print(f"No pre-computed data found at {intermediate_data_filepath}. Computing from scratch...")
    # --- Original Data Loading and Processing ---
    # Open the dataset with Dask chunking
    ds = xr.open_dataset(input_nc_filename, chunks={'lat': 500, 'lon': 500}) # Use parameterized input filename

    # Ensure the required variable exists
    if 'elevation' not in ds.data_vars:
        raise ValueError("Variable 'elevation' not found. Available variables: " + str(list(ds.data_vars)))
    elevation = ds['elevation']

    # Filter: keep only ocean (negative values)
    ocean_elevation = elevation.where(elevation < 0)

    # Downsample the data for easier handling
    ocean_coarse = ocean_elevation.coarsen(lat=10, lon=10, boundary="trim").mean()
    
    # Save the computed ocean_coarse for next time
    try:
        ocean_coarse.to_netcdf(intermediate_data_filepath)
        print(f"Saved 'ocean_coarse' data to {intermediate_data_filepath}")
    except Exception as e:
        print(f"Error saving 'ocean_coarse' data: {e}")

# --- Subsequent steps now use 'ocean_coarse' (either loaded or freshly computed) ---

# Define the color mapping parameters
vmin = -10919  # Minimum depth
vmax = 0       # Maximum (sea level)

# Create a normalization instance and get the desired cmocean colormap
norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
cmap = cmocean.cm.dense_r  # Reversed dense colormap for bathymetry

# Apply the colormap to the data.
rgba_array = cmap(norm(ocean_coarse.values))

# Convert the RGBA values from float (0.0-1.0) to 8-bit integers (0-255)
image_data = (rgba_array * 255).astype(np.uint8)

# Extract coordinate values from your data
lon_coords = ocean_coarse['lon'].values
lat_coords = ocean_coarse['lat'].values

# --- Data Integrity and Flipping Logic ---
if len(lon_coords) <= 1 or len(lat_coords) <= 1 or image_data.shape[1] == 0 or image_data.shape[0] == 0:
    raise ValueError(f"Insufficient dimensions after coarsening to define transform or write data. "
                     f"Image width: {image_data.shape[1]}, Image height: {image_data.shape[0]}, "
                     f"Num lon: {len(lon_coords)}, Num lat: {len(lat_coords)}")

image_data_to_write = image_data # Default to original
if len(lat_coords) > 1 and lat_coords[0] < lat_coords[-1]: # Check if latitudes are ascending (South to North)
    print("Latitude coordinates are ascending (South to North). Flipping image_data vertically for GeoTIFF writing.")
    image_data_to_write = np.flipud(image_data)
else:
    print("Latitude coordinates are descending (North to South) or single line. No vertical flip needed for image_data.")

# --- Transform Calculation ---
calculated_pixel_width = (lon_coords[-1] - lon_coords[0]) / (len(lon_coords) - 1)
positive_pixel_height = abs((lat_coords[0] - lat_coords[-1]) / (len(lat_coords) - 1))

ul_lon_corner = lon_coords.min() - calculated_pixel_width / 2
ul_lat_corner = lat_coords.max() + positive_pixel_height / 2

# IMPORTANT NOTE ON THE TRANSFORM:
# The line below uses 'positive_pixel_height'. In our previous discussions,
# we found that for standard GeoTIFFs (North-up, upper-left origin),
# a NEGATIVE y-pixel size is usually required in the transform, meaning you'd use
# '-positive_pixel_height'. You reported that removing the minus sign (using positive_pixel_height)
# eventually led to a working gdalinfo with a negative pixel size and functional gdal2tiles,
# in conjunction with the np.flipud() fix.
# Please ensure this transform results in the correct georeferencing
# (negative y-pixel size in gdalinfo and correct corner coordinates) for your specific setup.
# If gdalinfo shows a positive y-pixel size again with this, you likely need to revert to
# using -positive_pixel_height here.
my_final_transform = from_origin(ul_lon_corner,
                                 ul_lat_corner,
                                 calculated_pixel_width,
                                 positive_pixel_height) # Review this line based on your gdalinfo output

print(f"DEBUG: Calculated Affine transform to be used:\n{my_final_transform}")
print(f"DEBUG: Transform components: a={my_final_transform.a}, e={my_final_transform.e}")


# --- GeoTIFF Metadata and Writing ---
out_meta = {
    "driver": "GTiff",
    "height": image_data_to_write.shape[0],
    "width": image_data_to_write.shape[1],
    "count": 4,
    "dtype": 'uint8',
    "crs": "EPSG:4326",
    "transform": my_final_transform
}

try:
    with rasterio.open(output_tif_filename, "w", **out_meta) as dest: # Use dynamic output TIFF filename
        if image_data_to_write.ndim == 3 and image_data_to_write.shape[2] == out_meta["count"]:
            for i in range(out_meta["count"]):
                dest.write(image_data_to_write[:, :, i], i + 1)
        else:
            raise ValueError("image_data_to_write dimensions or count mismatch.")
    print(f"{output_tif_filename} has been written successfully.")
except Exception as e:
    print(f"Error writing GeoTIFF {output_tif_filename}: {e}")