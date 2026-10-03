import numpy as np
import matplotlib.pyplot as plt

rin = 0.5
rout = 1.0
label_size = 80
division_number = 100
distance_circle_to_label = -60

theta = np.linspace(0, 2*np.pi, division_number)
radius = np.array([rin, rout])
azimuth = np.linspace(0, 360, division_number)
Theta, Radius = np.meshgrid(theta, radius)
Azimuth = np.tile(azimuth, (2, 1))

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
for spine in ax.spines.values():
    spine.set_visible(False)
ax.pcolormesh(Theta, Radius, Azimuth, cmap='hsv', shading='gouraud')
ax.set_thetagrids([0, 90, 180, 270], labels=["W", "N", "E", "S"], color="black")
ax.set_yticklabels([])
ax.grid(False)
ax.tick_params(axis='x', labelsize=label_size, pad=distance_circle_to_label)

plt.show()