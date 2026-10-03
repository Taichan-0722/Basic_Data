import os
import re
import sys

def extract_parameters(line):

    panda_parameters = re.findall(r"[-+]?\d*\.\d+|\d+", line)

    if len(panda_parameters) != 8:
        raise ValueError(" the number of parameters is not correct (!= 8).\n")

    cx, cy, a1, a2, na, r1, r2, nr = map(float, panda_parameters)

    return cx, cy, a1, a2, int(na), r1, r2, int(nr)

def generate_panda_by_azimuth(a1, a2, na, r1, r2, nr):

    da = (a2 - a1) / na
    dr = (r2 - r1) / nr

    for i in range(nr):
        rr1 = r1 + dr * i
        rr2 = r1 + dr * (i + 1)

        for j in range(na):
            aa1 = a1 + da * j
            aa2 = a1 + da * (j + 1)

            yield aa1, aa2, rr1, rr2

def generate_panda_by_radius(a1, a2, na, r1, r2, nr):

    da = (a2 - a1) / na
    dr = (r2 - r1) / nr

    for i in range(na):
        aa1 = a1 + da * i
        aa2 = a1 + da * (i + 1)

        for j in range(nr):
            rr1 = r1 + dr * j
            rr2 = r1 + dr * (j + 1)

            yield aa1, aa2, rr1, rr2

if __name__ == "__main__":

    infile = input("\n Enter the inputting panda region file. > ").strip()

    if not os.path.exists(infile):
        print(f" '{infile}' does'nt exist.\n")
        sys.exit(1)

    outfile = input(" Enter the outputting panda region file. > ").strip()

    coordinate = input(" Enter the coordinate system. (physical|fk5) > ").strip()

    if coordinate not in ["physical", "fk5"]:
        print(" Please input 'physical' or 'fk5'.\n")
        sys.exit(1)

    sort_way = input(" Sort by azimuth or radius? (azimuth|radius) > ").strip()

    if sort_way not in ["azimuth", "radius"]:
        print(" Please input 'azimuth' or 'radius'.\n")
        sys.exit(1)

    with open(infile, "r") as f:
        lines = f.readlines()

    panda_line = None

    for line in lines:
        if line.strip().startswith("panda"):
            panda_line = line.split("#")[0].strip()
            break

    if panda_line is None:
        print(f" Cannot find panda row in '{infile}'.\n")
        sys.exit(1)

    cx, cy, a1, a2, na, r1, r2, nr = extract_parameters(panda_line)

    if sort_way == "azimuth":
        regions = list(generate_panda_by_azimuth(a1, a2, na, r1, r2, nr))

    elif sort_way == "radius":
        regions = list(generate_panda_by_radius(a1, a2, na, r1, r2, nr))

    with open(outfile, "w") as f:
        f.write(f"{coordinate}\n")

        if coordinate == "physical":
            for aa1, aa2, rr1, rr2 in regions:
                f.write(f'panda({cx},{cy},{aa1},{aa2},1,{rr1},{rr2},1)\n')

        elif coordinate == "fk5":
            for aa1, aa2, rr1, rr2 in regions:
                f.write(f'panda({cx},{cy},{aa1},{aa2},1,{rr1}",{rr2}",1)\n')

    print(f"\n '{len(regions)}' regions have written in '{outfile}' sorted by '{sort_way}' in '{coordinate}' coordinate.\n")