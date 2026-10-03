import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData
from scipy.optimize import differential_evolution

point_format = "o"
point_color = "c"
point_edge_color = "k"
point_edge_width = 0.5
errorbar_color = "grey"
errorbar_width = 1
errorbar_edge_length = 3

def power_model(parameters, r):
    K0, A, exp = parameters
    return K0 + A * r ** exp

csv_in = input("\n Enter the inputting csv file name. > ").strip()

if not os.path.exists(csv_in):
    print(f"\n '{csv_in}' does'nt exist.\n")
    sys.exit(1)

csv_out = input(" Enter the outputting csv file name. > ").strip()
xerr_or_not = input(" Consider the uncertainty in the x direction? (yes|no) > ").strip()
yerr_or_not = input(" Consider the uncertainty in the y direction? (yes|no) > ").strip()

df_in = pd.read_csv(csv_in)
index = df_in["Index"]
x_data = df_in["x"]
y_data = df_in["y"]
xerr_low = df_in["xerr_low"]
xerr_high = df_in["xerr_high"]
yerr_low = df_in["yerr_low"]
yerr_high = df_in["yerr_high"]

xerr_ave = (xerr_low + xerr_high) / 2
yerr_ave = (yerr_low + yerr_high) / 2

def residual(parameters):
    return np.sum((y_data - power_model(parameters, x_data)) ** 2)

bounds = [(-100, 100), (0, 100), (-10, 10)]

DE = differential_evolution(residual, bounds, maxiter=100000, popsize=50)

first_paremeters = DE.x

if xerr_or_not == "yes":

    if yerr_or_not == "yes":
        data = RealData(x_data, y_data, sx=xerr_ave, sy=yerr_ave)

    elif yerr_or_not == "no":
        data = RealData(x_data, y_data, sx=xerr_ave)

elif xerr_or_not == "no":

    if yerr_or_not == "yes":
        data = RealData(x_data, y_data, sy=yerr_ave)

    elif yerr_or_not == "no":
        data = RealData(x_data, y_data)

model = Model(power_model)
odr = ODR(data, model, beta0=first_paremeters)
output = odr.run()
second_parameters = output.beta

x_line = np.linspace(min(x_data), max(x_data), 100)
y_line = power_model(second_parameters, x_line)
y_fit = power_model(second_parameters, x_data)
ratio = y_data / y_fit

df_out = pd.DataFrame({"Index": index, "ratio_powerlaw": ratio})
df_out.to_csv(csv_out, index=False)

print("\n ----- First Fit (Differential Evolution) -----")
print(f"\n powerlaw model : K0 = {first_paremeters[0]:.5f}, A = {first_paremeters[1]:.5f}, exp = {first_paremeters[2]:.5f}")

print("\n ----- Second Fit (ODR) -----")
print(f"\n powerlaw model : K0 = {second_parameters[0]:.5f}, A = {second_parameters[1]:.5f}, exp = {second_parameters[2]:.5f}")

plt.plot(x_line, y_line, linestyle='--', color='r', label="powerlaw model")
plt.errorbar(x_data, y_data, xerr=xerr_ave, yerr=yerr_ave, fmt=point_format, color=point_color, mec=point_edge_color, mew=point_edge_width, ecolor=errorbar_color, elinewidth=errorbar_width, capsize=errorbar_edge_length)
plt.yscale('log')
plt.legend()
plt.show()