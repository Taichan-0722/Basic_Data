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

title = "Temperature Map"
colorbar_label = "Temperature (keV)"
colorbar_type = "inferno"
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
image_scale_min = 1E-5
scalebar2text = 20
transparency = 0.9

image_file = input("\n Input the image file name > ").strip()

if not os.path.exists(image_file):
    print(f"\n '{image_file}' doesn't exist.\n")
    sys.exit(1)

csv_file = input(" Input the csv file name > ").strip()

if not os.path.exists(csv_file):
    print(f"\n '{csv_file}' doesn't exist.\n")
    sys.exit(1)

colorbar_scale = input(" Input the color bar scale (lin|log) > ").strip()

if colorbar_scale not in ["lin", "log"]:
    print(" Please input 'lin' or 'log'.\n")
    sys.exit(1)

image_scale = input(" Input the image scale (lin|log) > ").strip()

if image_scale not in ["lin", "log"]:
    print(" Please input 'lin' or 'log'.\n")
    sys.exit(1)

with fits.open(image_file) as hdul:
    image_data = hdul[0].data.astype(float)
    wcs = WCS(hdul[0].header)

df = pd.read_csv(csv_file)
df = df.dropna(subset=["colorbar_value"])

plt.figure(figsize=(figure_size, figure_size))
ax = plt.subplot(projection=wcs)
ax.set_title(title)
ax.set_xlabel("Right Ascension (J2000)")
ax.set_ylabel("Declination (J2000)")
ax.coords[1].set_ticklabel(rotation=90)
ax.invert_xaxis()

center = SkyCoord(center_ra * u.deg, center_dec * u.deg, frame="fk5")
offset_ra = plot_range_ra * u.arcmin
offset_dec = plot_range_dec * u.arcmin

bottom_left = center.spherical_offsets_by(-offset_ra, -offset_dec).to_pixel(wcs)
top_right = center.spherical_offsets_by(offset_ra, offset_dec).to_pixel(wcs)

ax.set_xlim(top_right[0], bottom_left[0])
ax.set_ylim(bottom_left[1], top_right[1])

masked_image_data = np.ma.masked_less_equal(image_data, 0)

if image_scale == "lin":
    ax.imshow(masked_image_data, origin="lower", cmap="gray", norm=Normalize(vmin=image_scale_min, vmax=np.nanmax(image_data)))

elif image_scale == "log":
    ax.imshow(masked_image_data, origin="lower", cmap="gray", norm=LogNorm(vmin=image_scale_min, vmax=np.nanmax(image_data)))

colormap = plt.get_cmap(colorbar_type)

if colorbar_scale == "lin":
    norm = Normalize(vmin=colorbar_min, vmax=colorbar_max)

elif colorbar_scale == "log":
    norm = LogNorm(vmin=colorbar_min, vmax=colorbar_max)

for index, row in df.iterrows():

    regions = Regions.read(row["region_name"], format="ds9")
    color = colormap(norm(row["colorbar_value"]))

    for reg in regions:

        pix = reg.to_pixel(wcs)
        ax.add_artist(pix.as_artist(facecolor=color, edgecolor="black", alpha=transparency, fill=True))

sm = ScalarMappable(norm=norm, cmap=colormap)
colorbar = plt.colorbar(sm, ax=ax, pad=0.02)
colorbar.set_label(colorbar_label)
colorbar.ax.tick_params(which="both", direction="out")

scalebar_start = SkyCoord(scalebar_start_ra * u.deg, scalebar_start_dec * u.deg, frame="fk5")
scalebar_end = scalebar_start.spherical_offsets_by(scalebar_length * u.arcmin, 0 * u.arcmin)

x0, y0 = scalebar_start.to_pixel(wcs)
x1, y1 = scalebar_end.to_pixel(wcs)

ax.plot([x0, x1], [y0, y1], color="black", linewidth=scalebar_width)
ax.text((x0 + x1) / 2, y0 + scalebar2text, f"{scalebar_length} arcmin", ha="center", va="bottom", fontsize=scalebar_fontsize)

plt.show()