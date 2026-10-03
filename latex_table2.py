import os
import sys
import pandas as pd

csv_file = input("\n Enter the csv file name. > ").strip()

if not os.path.exists(csv_file):
    print(f"\n '{csv_file}' does'nt exist.\n")
    sys.exit(1)

out_file = input(" Enter the outputting file name. > ").strip()

df = pd.read_csv(csv_file)
table = df.to_latex(index=False, float_format="%.1f")
print()

with open(out_file, "w") as f:
    f.write(table)