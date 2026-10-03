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

def single_beta_model(parameters, r):
    S0, rc, beta = parameters
    return S0 * (1 + (r / rc) ** 2) ** (-3 * beta + 0.5)

def double_beta_model(parameters, r):
    S1, rc1, beta1, S2, rc2, beta2 = parameters
    return S1 * (1 + (r / rc1) ** 2) ** (-3 * beta1 + 0.5) + S2 * (1 + (r / rc2) ** 2) ** (-3 * beta2 + 0.5)

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

def residual_single(parameters):
    return np.sum((y_data - single_beta_model(parameters, x_data)) ** 2)

def residual_double(parameters):
    return np.sum((y_data - double_beta_model(parameters, x_data)) ** 2)

bounds_single = [(0, 100), (0.01, 100), (0, 10)]
bounds_double = [(0, 100), (0.01, 100), (0, 10), (0, 100), (0.01, 100), (0, 10)]

de_single = differential_evolution(residual_single, bounds_single, maxiter=100000, popsize=50)
de_double = differential_evolution(residual_double, bounds_double, maxiter=100000, popsize=50)

first_p_single = de_single.x
first_p_double = de_double.x

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

model_single = Model(single_beta_model)
model_double = Model(double_beta_model)

odr_single = ODR(data, model_single, beta0=first_p_single)
odr_double = ODR(data, model_double, beta0=first_p_double)

output_single = odr_single.run()
output_double = odr_double.run()

second_p_single = output_single.beta
second_p_double = output_double.beta

x_line = np.linspace(min(x_data), max(x_data), 100)

y_line_single = single_beta_model(second_p_single, x_line)
y_line_double = double_beta_model(second_p_double, x_line)

y_fit_single = single_beta_model(second_p_single, x_data)
y_fit_double = double_beta_model(second_p_double, x_data)

ratio_single = y_data / y_fit_single
ratio_double = y_data / y_fit_double

df_out = pd.DataFrame({"Index": index, "ratio_single_model": ratio_single, "ratio_double_model": ratio_double})
df_out.to_csv(csv_out, index=False)

print("\n ----- First Fit (Differential Evolution) -----")
print(f"\n single beta model : S0 = {first_p_single[0]:.5f}, rc = {first_p_single[1]:.5f}, beta = {first_p_single[2]:.5f}")
print(f"\n double beta model : S1 = {first_p_double[0]:.5f}, rc1 = {first_p_double[1]:.5f}, beta1 = {first_p_double[2]:.5f}")
print(f"                     S2 = {first_p_double[3]:.5f}, rc2 = {first_p_double[4]:.5f}, beta2 = {first_p_double[5]:.5f}\n")

print("\n ----- Second Fit (ODR) -----")
print(f"\n single beta model : S0 = {second_p_single[0]:.5f}, rc = {second_p_single[1]:.5f}, beta = {second_p_single[2]:.5f}")
print(f"\n double beta model : S1 = {second_p_double[0]:.5f}, rc1 = {second_p_double[1]:.5f}, beta1 = {second_p_double[2]:.5f}")
print(f"                     S2 = {second_p_double[3]:.5f}, rc2 = {second_p_double[4]:.5f}, beta2 = {second_p_double[5]:.5f}\n")

plt.plot(x_line, y_line_single, linestyle='--', color='b', label="single beta model")
plt.plot(x_line, y_line_double, linestyle='--', color='r', label="double beta model")
plt.errorbar(x_data, y_data, xerr=xerr_ave, yerr=yerr_ave, fmt=point_format, color=point_color, mec=point_edge_color, mew=point_edge_width, ecolor=errorbar_color, elinewidth=errorbar_width, capsize=errorbar_edge_length)
plt.yscale('log')
plt.legend()
plt.show()