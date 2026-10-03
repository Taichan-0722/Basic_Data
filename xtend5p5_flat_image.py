import os
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions

S_const = 100.0

count_image = input("\n Enter the inputting count image file name > ").strip()

if not os.path.exists(count_image):
    print(f"\n '{count_image}' does'nt exist.\n")
    sys.exit(1)

region_file = input(" Enter the inputting region file name > ").strip()

if not os.path.exists(region_file):
    print(f"\n '{region_file}' does'nt exist.\n")
    sys.exit(1)

regions = Regions.read(region_file, format="ds9")

if len(regions) != 1:
    print(f"\n '{region_file}' doesn't contain exactly one region.\n")
    sys.exit(1)

flat_image_file = input(" Enter the outputting flat image file name > ").strip()

with fits.open(count_image) as hdul:
    count_data = hdul[0].data.astype(float)
    header = hdul[0].header
    wcs = WCS(header)

num_y, num_x = count_data.shape
flat_image = np.full((num_y, num_x), S_const)
region = regions[0]
mask = region.to_pixel(wcs).to_mask(mode="center")
mask_image = mask.to_image((num_y, num_x))
mask_image = np.nan_to_num(mask_image)
mask_data = flat_image * mask_image

header_crop = header.copy()
fits.PrimaryHDU(mask_data, header=header_crop).writeto(flat_image_file, overwrite=True)

print("\n All processes have been completed.\n")