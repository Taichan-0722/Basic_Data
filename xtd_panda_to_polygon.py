import os
import re
import sys
import math
import subprocess

def panda_to_polygon(line, division_number):

    panda_parameters = re.findall(r"[-+]?\d*\.\d+|\d+", line)

    if len(panda_parameters) != 8:
        raise ValueError(" the number of parameters is not correct (!= 8).\n")

    cx, cy, a1_deg, a2_deg, na, r1_arcsec, r2_arcsec, nr = map(float, panda_parameters)

    r1_deg = r1_arcsec * 1.0 / 3600.0
    r2_deg = r2_arcsec * 1.0 / 3600.0
    a1_rad = math.radians(a1_deg)
    a2_rad = math.radians(a2_deg)

    angles = [a1_rad + (a2_rad - a1_rad) * i / division_number for i in range(division_number + 1)]

    vertices = []

    for angle in angles:
        y = cy + r1_deg * math.sin(angle)
        dec_rad = math.radians(y)
        x = cx - r1_deg * math.cos(angle) / math.cos(dec_rad)
        vertices.append((x, y))

    for angle in reversed(angles):
        y = cy + r2_deg * math.sin(angle)
        dec_rad = math.radians(y)
        x = cx - r2_deg * math.cos(angle) / math.cos(dec_rad)
        vertices.append((x, y))

    return "polygon(" + ",".join(f"{x:.8f},{y:.8f}" for x, y in vertices) + ")"

if __name__ == "__main__":

    OBSID = os.environ.get("OBSID")
    SUB = os.environ.get("SUB")

    input_file_name = f"xa{OBSID}xtd_{SUB}_src_"
    output_file_name = f"xa{OBSID}xtd_{SUB}_pol_"
    start_number = 1
    end_number = 131
    division_number = 100

    print()

    for i in range(start_number, end_number+1):

        if not os.path.exists(f"{input_file_name}{i}.reg"):
            print(f" '{input_file_name}{i}.reg' does'nt exist.\n")
            sys.exit(1)

        with open(f"{input_file_name}{i}.reg", "r") as f:
            lines = f.readlines()

        panda_line = None

        for line in lines:
            if line.strip().startswith("panda"):
                panda_line = line.split("#")[0].strip()
                break

        if panda_line is None:
            subprocess.run(["cp", f"{input_file_name}{i}.reg", f"{output_file_name}{i}.reg"])
            print(f" Cannot find panda row in '{input_file_name}{i}.reg'. Copied the original file to '{output_file_name}{i}.reg'.")
            continue

        with open(f"{output_file_name}{i}.reg", "w") as f:
            f.write("fk5\n")
            f.write(panda_to_polygon(panda_line, division_number))

        print(f" '{input_file_name}{i}.reg' has been converted to '{output_file_name}{i}.reg'")