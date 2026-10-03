import os
import sys
import shutil
import numpy as np
import pandas as pd
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from matplotlib.patches import Annulus
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions

colorbar_label = "Temperature (keV)"
colorbar_type = "gnuplot2"
colorbar_min = 0.0
colorbar_max = 8.0
colorbar_labelsize = 55
colorbar_fontsize = 45
contour_color = "white"
contour_linewidth = 1
contour_alpha = 1.0
velocity_step = 30
velocity_color = "grey"
velocity_width = 0.003
velocity_scale = 500
center_x = 230   # kpc
center_y = -230   # kpc
plot_range_x = 2500   # kpc
plot_range_y = 2500   # kpc
scalebar_length = 1000   # kpc
scalebar_width = 3
scalebar_start_x = -2000   # kpc
scalebar_start_y = -2500   # kpc
scalebar_fontsize = 50
scalebar2text = 20   # pixels
time_text_x = -800   # kpc
time_text_y = 2000   # kpc
time_text_fontsize = 50
figure_size = 16
dpi = 300
fps = 3
animation_name = "movie.mp4"
frame_directory_name = "frames"

def plot_contours(filename, ax, wcs, color, linewidth):

    contours = []
    coords = []
    level = None

    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:

        line = line.strip()

        if (not line or line.startswith("#") or line.startswith("global") or line == "wcs"):
            continue

        if line.startswith("level="):
            level = float(line.split("=")[1].strip())
            continue

        if line == "(":
            coords = []
            continue

        if line == ")":
            if len(coords) >= 2:
                contours.append((level, np.array(coords)))
            coords = []
            continue

        coord_values = line.split()

        if len(coord_values) >= 2:
            x = float(coord_values[0])
            y = float(coord_values[1])
            coords.append((x, y))

    for level, coords in contours:
        x = coords[:, 0]
        y = coords[:, 1]
        x_pix, y_pix = wcs.world_to_pixel_values(x, y)
        ax.plot(x_pix, y_pix, color=color, linewidth=linewidth, zorder=1, alpha=contour_alpha)


contour_or_not = input("\n Plot contours? (yes|no) > ").strip()
velocity_or_not = input(" Display the velocity field with arrows? (yes|no) > ").strip()
csv_file = input(" Input the CSV file name > ").strip()

if not os.path.exists(csv_file):
    print(f"\n '{csv_file}' doesn't exist.\n")
    sys.exit(1)

print()
df = pd.read_csv(csv_file)
fits_list = df["fits_file"]
contour_list = df["contour_file"]
time_list = df["time"]
param_index_list = df["param_index"]
Vx_index_list = df["Vx_index"]
Vy_index_list = df["Vy_index"]

os.makedirs(f"{frame_directory_name}", exist_ok=True)

for f in os.listdir(f"{frame_directory_name}/"):
    os.remove(os.path.join(f"{frame_directory_name}/", f))

with imageio.get_writer(animation_name, fps=fps) as writer:

    for fits_file, contour_file, time, param_index, Vx_index, Vy_index in zip(fits_list, contour_list, time_list, param_index_list, Vx_index_list, Vy_index_list):

        if not os.path.exists(fits_file):
            print(f"\n '{fits_file}' doesn't exist.\n")
            sys.exit(1)

        with fits.open(fits_file) as hdul:
            image_data = np.squeeze(hdul[int(param_index)].data.astype(float))
            wcs = WCS(hdul[int(param_index)].header, naxis=2)
            if velocity_or_not == "yes":
                Vx_data = np.squeeze(hdul[int(Vx_index)].data.astype(float))
                Vy_data = np.squeeze(hdul[int(Vy_index)].data.astype(float))

        plt.figure(figsize=(figure_size, figure_size))
        ax = plt.subplot(projection=wcs)
        ax.axis("off")

        bottom_left = np.array([[center_x - plot_range_x, center_y - plot_range_y]])
        top_right = np.array([[center_x + plot_range_x, center_y + plot_range_y]])

        pixel_bottom_left = wcs.wcs_world2pix(bottom_left, 0)[0]
        pixel_top_right = wcs.wcs_world2pix(top_right, 0)[0]

        x_left, y_bottom = pixel_bottom_left
        x_right, y_top = pixel_top_right

        ax.set_xlim(min(x_left, x_right), max(x_left, x_right))
        ax.set_ylim(min(y_bottom, y_top), max(y_bottom, y_top))

        norm = Normalize(vmin=colorbar_min, vmax=colorbar_max)
        # norm = Normalize(vmin=colorbar_min, vmax=np.nanmax(image_data))
        ax.imshow(image_data, origin="lower", cmap=colorbar_type, norm=norm)

        sm = ScalarMappable(norm=norm, cmap=colorbar_type)
        colorbar = plt.colorbar(sm, ax=ax, pad=0.02, shrink=0.885)
        colorbar.set_label(colorbar_label, labelpad=30, fontsize=colorbar_labelsize)
        colorbar.ax.tick_params(which="both", direction="out", pad=10, labelsize=colorbar_fontsize)

        if contour_or_not == "yes":

            if not os.path.exists(contour_file):
                print(f"\n '{contour_file}' doesn't exist.\n")
                sys.exit(1)

            plot_contours(contour_file, ax, wcs, color=contour_color, linewidth=contour_linewidth)

        if velocity_or_not == "yes":
            ny, nx = image_data.shape
            y, x = np.mgrid[0:ny:velocity_step, 0:nx:velocity_step]
            Vx = Vx_data[::velocity_step, ::velocity_step]
            Vy = Vy_data[::velocity_step, ::velocity_step]
            V = np.sqrt(Vx**2 + Vy**2)
            Vx_unit = np.divide(Vx, V, out=np.zeros_like(Vx), where=V != 0)
            Vy_unit = np.divide(Vy, V, out=np.zeros_like(Vy), where=V != 0)
            length = np.sqrt(V)
            vx_plot = Vx_unit * length
            vy_plot = Vy_unit * length
            ax.quiver(x, y, vx_plot, vy_plot, color=velocity_color, scale=velocity_scale, width=velocity_width, zorder=2)

        scalebar_start = np.array([[scalebar_start_x, scalebar_start_y]])
        scalebar_end = np.array([[scalebar_start_x + scalebar_length, scalebar_start_y]])

        scalebar_start_pixel = wcs.wcs_world2pix(scalebar_start, 0)[0]
        scalebar_end_pixel = wcs.wcs_world2pix(scalebar_end, 0)[0]

        x0, y0 = scalebar_start_pixel
        x1, y1 = scalebar_end_pixel

        ax.plot([x0, x1], [y0, y1], color="white", linewidth=scalebar_width)
        ax.text((x0 + x1) / 2, y0 + scalebar2text, f"{scalebar_length} kpc", ha="center", va="bottom", fontsize=scalebar_fontsize, backgroundcolor="none", color="white")

        time_text = np.array([[time_text_x, time_text_y]])
        time_text_pixel = wcs.wcs_world2pix(time_text, 0)[0]
        ax.text(time_text_pixel[0], time_text_pixel[1], f"{float(time):.2f} Gyr", ha="right", va="top", fontsize=time_text_fontsize, backgroundcolor="none", color="yellow")

        point = np.array([[center_x, center_y]])
        point_pixel = wcs.wcs_world2pix(point, 0)[0]
        ax.add_patch(Annulus((point_pixel[0], point_pixel[1]), r=15.0, width=3, facecolor='lime', edgecolor='black', lw=0, zorder=3))

        frame_name = f"frame_{float(time):.2f}_ctr".replace(".", "p") + ".png"
        plt.savefig(frame_name, dpi=dpi)
        plt.close()

        writer.append_data(imageio.imread(frame_name))
        shutil.move(frame_name, f"{frame_directory_name}/")
        print(f" '{fits_file}' has been done.")