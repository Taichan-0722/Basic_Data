import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from regions import Regions
from reproject import reproject_interp

title = "HOGE Map"
colorbar_label = "HOGE (hoge)"
colorbar_type = "inferno"
optical_image_color = "Greys_r"
count_image_color = "gnuplot2"
colorbar_min = 2.5
colorbar_max = 5.5
center_ra = 192.167   # degree
center_dec = -41.3119031   # degree
scalebar_length = 10   # arcmin
scalebar_width = 1
scalebar_start_ra = 192.54   # degree
scalebar_start_dec = -41.77   # degree
scalebar_fontsize = 20
plot_range_ra = 33   # arcmin
plot_range_dec = 30   # arcmin
figure_size = 7
optical_scale_min = 1000
count_scale_min = 1.0e-1
color_map_transparency = 1.0
count_image_transparency = 0.6
scalebar2text = 20

optical_image_file = input("\n Input the optical image file name > ").strip()

if not os.path.exists(optical_image_file):
    print(f"\n '{optical_image_file}' doesn't exist.\n")
    sys.exit(1)

with fits.open(optical_image_file) as hdul:
    optical_image_data = np.squeeze(hdul[0].data.astype(float))
    optical_wcs = WCS(hdul[0].header).celestial

count_image_file = input("\n Input the count image file name > ").strip()

if not os.path.exists(count_image_file):
    print(f"\n '{count_image_file}' doesn't exist.\n")
    sys.exit(1)

with fits.open(count_image_file) as hdul:
    count_image_data = np.squeeze(hdul[0].data.astype(float))
    count_wcs = WCS(hdul[0].header).celestial

csv_file = input(" Input the csv file name > ").strip()

if not os.path.exists(csv_file):
    print(f"\n '{csv_file}' doesn't exist.\n")
    sys.exit(1)

df = pd.read_csv(csv_file)
df = df.dropna(subset=["colorbar_value"])

optical_image_reprojected, footprint = reproject_interp((optical_image_data, optical_wcs), count_wcs, shape_out=count_image_data.shape)
optical_image_reprojected[footprint <= 0] = np.nan

plt.figure(figsize=(figure_size, figure_size))
ax = plt.subplot(projection=count_wcs)
ax.set_title(title)
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
ax.coords[1].set_ticklabel(rotation=90)
ax.invert_xaxis()

center = SkyCoord(center_ra * u.deg, center_dec * u.deg, frame="fk5")
offset_ra = plot_range_ra * u.arcmin
offset_dec = plot_range_dec * u.arcmin

bottom_left = center.spherical_offsets_by(-offset_ra, -offset_dec).to_pixel(count_wcs)
top_right = center.spherical_offsets_by(offset_ra, offset_dec).to_pixel(count_wcs)

ax.set_xlim(top_right[0], bottom_left[0])
ax.set_ylim(bottom_left[1], top_right[1])

masked_optical_image_data = np.ma.masked_less_equal(optical_image_reprojected, 0)
ax.imshow(masked_optical_image_data, origin="lower", cmap=optical_image_color, norm=Normalize(vmin=optical_scale_min, vmax=np.nanmax(optical_image_reprojected)))

masked_count_image_data = np.ma.masked_less_equal(count_image_data, 0)
ax.imshow(masked_count_image_data, origin="lower", cmap=count_image_color, norm=LogNorm(vmin=count_scale_min, vmax=np.nanmax(count_image_data)), alpha=count_image_transparency)

norm = Normalize(vmin=colorbar_min, vmax=colorbar_max)
colormap = plt.get_cmap(colorbar_type)

for index, row in df.iterrows():

    regions = Regions.read(row["region_name"], format="ds9")
    color = colormap(norm(row["colorbar_value"]))

    for reg in regions:

        pix = reg.to_pixel(count_wcs)
        ax.add_artist(pix.as_artist(facecolor=color, edgecolor="yellow", alpha=color_map_transparency, fill=False))

sm = ScalarMappable(norm=norm, cmap=colormap)
colorbar = plt.colorbar(sm, ax=ax, pad=0.02)
colorbar.set_label(colorbar_label)
colorbar.ax.tick_params(which="both", direction="out")

scalebar_start = SkyCoord(scalebar_start_ra * u.deg, scalebar_start_dec * u.deg, frame="fk5")
scalebar_end = scalebar_start.spherical_offsets_by(scalebar_length * u.arcmin, 0 * u.arcmin)

x0, y0 = scalebar_start.to_pixel(count_wcs)
x1, y1 = scalebar_end.to_pixel(count_wcs)

ax.plot([x0, x1], [y0, y1], color="white", linewidth=scalebar_width)
ax.text((x0 + x1) / 2, y0 + scalebar2text, f"{scalebar_length} arcmin", ha="center", va="bottom", fontsize=scalebar_fontsize, color="white", backgroundcolor="none")

plt.show()