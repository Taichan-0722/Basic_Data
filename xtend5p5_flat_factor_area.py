import os
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
from astropy.wcs.utils import proj_plane_pixel_scales

OBSID = os.environ.get("OBSID")
SUB = os.environ.get("SUB")
region_basename = f"xa{OBSID}xtd_{SUB}_pol_"
division_number = 20

model_image = input("\n Enter the inputting model image file name > ").strip()

if not os.path.exists(model_image):
    print(f"\n '{model_image}' does'nt exist.\n")
    sys.exit(1)

exclusion_region_file = input(" Enter the exclusion region file (or blank) > ").strip()

if exclusion_region_file and not os.path.exists(exclusion_region_file):
    print(f"\n '{exclusion_region_file}' doesn't exist.\n")
    sys.exit(1)

exclusion_regions = None

if exclusion_region_file:
    exclusion_regions = Regions.read(exclusion_region_file, format="ds9")

start = int(input(" Input the starting region number. > ").strip())
finish = int(input(" Input the finishing region number. > ").strip())

with fits.open(model_image) as hdul:
    count_data = hdul[0].data.astype(float)
    header = hdul[0].header
    wcs = WCS(header)

pixscale = proj_plane_pixel_scales(wcs) * 60
arcmin_per_pix_x = pixscale[0]
arcmin_per_pix_y = pixscale[1]
arcmin2_per_pixel = arcmin_per_pix_x * arcmin_per_pix_y

total_counts = np.nansum(count_data)

if exclusion_regions:

    exclusion_mask_total = np.zeros_like(count_data)

    for exclusion_region in exclusion_regions:
        exclusion_mask = exclusion_region.to_pixel(wcs).to_mask(mode="subpixels", subpixels=division_number)
        exclusion_mask_data = exclusion_mask.to_image(count_data.shape)
        exclusion_mask_total += exclusion_mask_data

    exclusion_mask_total = np.clip(exclusion_mask_total, 0, 1)

with open("flat_norm_factor.list", "w", newline="") as factor_list, open("region_area.list", "w") as area_list:

    for i in range(1, start):
        factor_list.write("0.0\n")
        area_list.write("0.0\n")

    for i in range(start, finish + 1):

        region_file = f"{region_basename}{i}.reg"

        if not os.path.exists(region_file):
            print(f"\n '{region_file}' does'nt exist.\n")
            sys.exit(1)

        regions = Regions.read(region_file, format="ds9")

        if len(regions) != 1:
            print(f"\n '{region_file}' doesn't contain exactly one region.\n")
            sys.exit(1)

        region = regions[0]
        mask = region.to_pixel(wcs).to_mask(mode="subpixels", subpixels=division_number)

# factor (counts in inclusion region / counts in ARF region)

        mask_data = mask.multiply(count_data, fill_value=0.0)
        region_counts = np.nansum(mask_data)
        factor = region_counts / total_counts

        factor_list.write(f"{factor}\n")

# region area (inclusion region - exclusion region)

        mask_data = mask.to_image(count_data.shape)

        if exclusion_regions:
            excluded_mask_data = mask_data * (1 - exclusion_mask_total)

        elif exclusion_regions is None:
            excluded_mask_data = mask_data

        excluded_mask_data = np.nan_to_num(excluded_mask_data, nan=0.0)
        pixels = np.nansum(excluded_mask_data)
        area_arcmin2 = pixels * arcmin2_per_pixel

        area_list.write(f"{area_arcmin2}\n")

        print(f" {region_file} has been done.")

print()