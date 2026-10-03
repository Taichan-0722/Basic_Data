import os
import sys

infile = input("\n Enter the panda region name. > ").strip()

if not os.path.exists(infile):
    print(f" '{infile}' does'nt exist.\n")
    sys.exit(1)

coordinate = input(" Enter the coordinate system. (physical|fk5) > ").strip()

if coordinate not in ["physical", "fk5"]:
    print(" Please input 'physical' or 'fk5'.\n")
    sys.exit(1)

start_number = int(input(" Enter the starting number. > ").strip())
sort_number = int(input(" Enter the sort number (If this is fine as it is, please enter '1'). > ").strip())
print("")

OBSID = os.environ.get("OBSID")
SUB = os.environ.get("SUB")

with open(infile, "r") as file:
    lines = file.readlines()

reg_lines = [line.strip() for line in lines if line.strip().startswith("panda")]
total = len(reg_lines)
number_rows = total // sort_number
sort_lines = []

for column in range(sort_number):
    for row in range(number_rows):
        id = row * sort_number + column
        sort_lines.append(reg_lines[id])

for i, region in enumerate(sort_lines, start_number):

    if "#" in region:
        region = region.split("#")[0].strip()

    region_with_number = f"{region} #text={{{i}}}"

    outfile = f"xa{OBSID}xtd_{SUB}_src_{i}.reg"

    with open(outfile, "w") as file:
        file.write(f"{coordinate}")
        file.write("\n")
        file.write(region_with_number)

    print(f' "{outfile}" has been created.')

print("")