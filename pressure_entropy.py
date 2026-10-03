import os
import sys
import numpy as np
import pandas as pd

def pressure(kT, norm):
    return kT * norm ** (1/2)

def entropy(kT, norm):
    return kT / norm ** (1/3)

csv_in = input("\n Enter the inputting csv file name. > ").strip()

if not os.path.exists(csv_in):
    print(f"\n '{csv_in}' doesn't exist.\n")
    sys.exit(1)

csv_out = input(" Enter the outputting csv file name. > ").strip()

df_in = pd.read_csv(csv_in)
index_list = df_in["Index"]
kT_list = df_in["kT"]
kT_error_low_list = df_in["kT_error_low"]
kT_error_high_list = df_in["kT_error_high"]
norm_list = df_in["norm"]
norm_error_low_list = df_in["norm_error_low"]
norm_error_high_list = df_in["norm_error_high"]

pressure_list = pressure(kT_list, norm_list)
pressure_error_low_list = np.sqrt((pressure(kT_list - kT_error_low_list, norm_list) - pressure_list) ** 2 + (pressure(kT_list, norm_list - norm_error_low_list) - pressure_list ) ** 2)
pressure_error_high_list = np.sqrt((pressure(kT_list + kT_error_high_list, norm_list) - pressure_list) ** 2 + (pressure(kT_list, norm_list + norm_error_high_list) - pressure_list ) ** 2)

entropy_list = entropy(kT_list, norm_list)
entropy_error_low_list = np.sqrt((entropy(kT_list - kT_error_low_list, norm_list) - entropy_list) ** 2 + (entropy(kT_list, norm_list + norm_error_high_list) - entropy_list ) ** 2)
entropy_error_high_list = np.sqrt((entropy(kT_list + kT_error_high_list, norm_list) - entropy_list) ** 2 + (entropy(kT_list, norm_list - norm_error_low_list) - entropy_list ) ** 2)

df_out = pd.DataFrame({"Index": index_list, "pressure": pressure_list, "pressure_error_low": pressure_error_low_list, "pressure_error_high": pressure_error_high_list, "entropy": entropy_list, "entropy_error_low": entropy_error_low_list, "entropy_error_high": entropy_error_high_list})
df_out.to_csv(csv_out, index=False)