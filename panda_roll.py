import os
import re
import sys

def extract_parameters(line):

    panda_parameters = re.findall(r"[-+]?\d*\.\d+|\d+", line)

    if len(panda_parameters) != 8:
        raise ValueError("\n the number of parameters is not correct (!= 8).\n")

    cx, cy, a1, a2, na, r1, r2, nr = map(float, panda_parameters)

    return cx, cy, a1, a2, int(na), r1, r2, int(nr)

if __name__ == "__main__":

    OBSID = os.environ.get("OBSID")
    SUB = os.environ.get("SUB")
    PA_NOM = float(os.environ.get("PA_NOM"))

    input_file_name = f"xa{OBSID}xtd_{SUB}_src_"
    output_file_name = f"xa{OBSID}xtd_{SUB}_roll_"
    start_number = 1
    end_number = 80

    print()

    for i in range(start_number, end_number+1):

        if not os.path.exists(f"{input_file_name}{i}.reg"):
            print(f"\n '{input_file_name}{i}.reg' does'nt exist.\n")
            sys.exit(1)

        with open(f"{input_file_name}{i}.reg", "r") as f:
            lines = f.readlines()

        panda_line = None

        for line in lines:
            if line.strip().startswith("panda"):
                panda_line = line.split("#")[0].strip()
                break

        if panda_line is None:
            print(f" Cannot find panda row in '{input_file_name}{i}.reg'.")
            continue

        cx, cy, a1, a2, na, r1, r2, nr = extract_parameters(panda_line)

        if na != 1 or nr != 1:
            print(f"\n The number of azimuthal or radial divisions in '{input_file_name}{i}.reg' is not 1.\n")
            sys.exit(1)

        a1 = (a1 - PA_NOM) % 360.0
        a2 = (a2 - PA_NOM) % 360.0

        with open(f"{output_file_name}{i}.reg", "w") as f:
            f.write("fk5\n")
            f.write(f'panda({cx},{cy},{a1},{a2},1,{r1}",{r2}",1)')

        print(f" '{input_file_name}{i}.reg' has been rotated by {PA_NOM} degrees to produce '{output_file_name}{i}.reg'")