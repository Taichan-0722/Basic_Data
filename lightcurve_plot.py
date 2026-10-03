import os
import sys
from astropy.io import fits
import matplotlib.pyplot as plt

marker_size = 1
marker_color = "dodgerblue"

light_curve = input("\n Input the name of light curve > ").strip()

if not os.path.exists(light_curve):
    print(f"\n '{light_curve}' does'nt exist.\n")
    sys.exit(1)

print()
hdul = fits.open(light_curve)
data = hdul[1].data

time = data["TIME"]
rate = data["RATE"]
err  = data["ERROR"]

plt.errorbar(time, rate, yerr=err, fmt='o', ms=marker_size, color=marker_color)
plt.xlabel("Time (s)")
plt.ylabel("Count rate (count/s)")
plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
plt.show()