import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions

image_color = "gnuplot2"
center_detx = 0.0
center_dety = 0.0
plot_range = 35
figure_size = 7
image_scale_min = 1.0e-1
image_transparency = 1.0

image_file = input("\n Input the count image file name > ").strip()

if not os.path.exists(image_file):
    print(f"\n '{image_file}' doesn't exist.\n")
    sys.exit(1)

with fits.open(image_file) as hdul:
    image_data = np.squeeze(hdul[0].data.astype(float))
    wcs = WCS(hdul[0].header, naxis=2)

plt.figure(figsize=(figure_size, figure_size))
ax = plt.subplot(projection=wcs)
ax.set_facecolor("black")
ax.set_xlabel("DETX")
ax.set_ylabel("DETY")
ax.coords[1].set_ticklabel(rotation=90)

det_bottom_left = np.array([[center_detx - plot_range, center_dety - plot_range,]])
det_top_right = np.array([[center_detx + plot_range,center_dety + plot_range,]])

pixel_bottom_left = wcs.wcs_world2pix(det_bottom_left, 0)[0]
pixel_top_right = wcs.wcs_world2pix(det_top_right, 0)[0]

ax.set_xlim(min(pixel_bottom_left[0], pixel_top_right[0]), max(pixel_bottom_left[0], pixel_top_right[0]))
ax.set_ylim(min(pixel_bottom_left[1], pixel_top_right[1]), max(pixel_bottom_left[1], pixel_top_right[1]))

norm = LogNorm(vmin=image_scale_min, vmax=np.nanmax(image_data))
masked_image_data = np.ma.masked_less_equal(image_data, 0)
ax.imshow(masked_image_data, origin="lower", cmap=image_color, norm=norm, alpha=image_transparency)

sm = ScalarMappable(norm=norm, cmap=image_color)
colorbar = plt.colorbar(sm, ax=ax, pad=0.02)
colorbar.set_label("Counts")
colorbar.ax.tick_params(which="both", direction="out")

plt.show()