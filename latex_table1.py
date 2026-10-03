import os
import sys
import pandas as pd

csv_file = input("\n Enter the csv file name. > ").strip()

if not os.path.exists(csv_file):
    print(f"\n '{csv_file}' does'nt exist.\n")
    sys.exit(1)

out_file = input(" Enter the outputting file name. > ").strip()

decimal_places = int(input(" How many decimal places do you want to display? > ").strip())

multiplier = float(input(" How many times do you want to multiply the values? > ").strip())

df = pd.read_csv(csv_file)
value = df["value"] * multiplier
error_low = df["error_low"] * multiplier
error_high = df["error_high"] * multiplier
print()


with open(out_file, "w") as f:

    for value_i, error_low_i, error_high_i in zip(value, error_low, error_high):

        if pd.isna(error_low_i) or pd.isna(error_high_i):
            f.write(f"${value_i:.{decimal_places}f}$\n")

        else:
            f.write(f"${value_i:.{decimal_places}f}_{{-{error_low_i:.{decimal_places}f}}}^{{+{error_high_i:.{decimal_places}f}}}$\n")


if not multiplier == 1.0:
    print(f" !! Warning : values have been multiplied by '{multiplier}'. Be careful. !!\n")