import os
import sys
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from regions import Regions

S0 = 100.0
rc_arcmin = 2.763565   # arcmin
beta = 0.430815
ra_center  = 192.2045649   # degree
dec_center = -41.3118720   # degree

def single_beta_model(S0, rc, beta, r):
    return S0 * (1 + (r / rc) ** 2) ** (-3 * beta + 0.5)

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

model_image = input(" Enter the outputting model image file name > ").strip()

with fits.open(count_image) as hdul:
    count_data = hdul[0].data.astype(float)
    header = hdul[0].header
    wcs = WCS(header)

num_y, num_x = count_data.shape
y, x = np.indices((num_y, num_x))
ra, dec = wcs.wcs_pix2world(x, y, 0)
radec = SkyCoord(ra * u.deg, dec * u.deg, frame="fk5")
radec_center = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame="fk5")
r_arcmin = radec.separation(radec_center).arcmin
model_function = single_beta_model(S0, rc_arcmin, beta, r_arcmin)

region = regions[0]
mask = region.to_pixel(wcs).to_mask(mode="center")
mask_image = mask.to_image((num_y, num_x))
mask_image = np.nan_to_num(mask_image)
mask_data = model_function * mask_image

header_crop = header.copy()
fits.PrimaryHDU(mask_data, header=header_crop).writeto(model_image, overwrite=True)

print("\n All processes have been completed.\n")