import os
import sys
import re
import io
import csv
import glob
import subprocess
import shutil
import math
import numpy as np
import pandas as pd
from xspec import *

nH_value = 0.0778
Redshift_value = 0.0114
Redshift_low = 0.008
Redshift_high = 0.014
CXB_norm_400pi = 1.0E-3
minSig_value = 5
maxBins_value = 10
y_range_value = 1.0E-5
raytrace_source_area = 44 * 44   # arcmin^2
title = "None"
model = "1vapec"
version = "v5"

OBSID = os.environ.get("OBSID")
SUB = os.environ.get("SUB")

start = int(input("\n Input the starting region number. > ").strip())
finish = int(input(" Input the finishing region number. > ").strip())
error_do_or_not = input(" Execute the error command? (yes|no) > ").strip()

if error_do_or_not not in ["yes", "no"]:
    print(" Please input 'yes' or 'no'.\n")
    sys.exit(1)

if error_do_or_not == "yes":
    csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
        "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
        "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
        "src_Abund_He",
        "src_Abund_C",
        "src_Abund_N",
        "src_Abund_O",
        "src_Abund_Ne",
        "src_Abund_Mg",
        "src_Abund_Al",
        "src_Abund_Si",
        "src_Abund_S",
        "src_Abund_Ar",
        "src_Abund_Ca",
        "src_Abund_Fe",
        "src_Abund_Ni",
        "NXB_gain_slope",
        "C-stat",
        "d.o.f",
        "C-stat/d.o.f"
    ]

elif error_do_or_not == "no":
    csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "src_kT(keV)",
        "src_norm(/400pi_arcmin2)",
        "src_Redshift",
        "src_Abund_He",
        "src_Abund_C",
        "src_Abund_N",
        "src_Abund_O",
        "src_Abund_Ne",
        "src_Abund_Mg",
        "src_Abund_Al",
        "src_Abund_Si",
        "src_Abund_S",
        "src_Abund_Ar",
        "src_Abund_Ca",
        "src_Abund_Fe",
        "src_Abund_Ni",
        "NXB_gain_slope",
        "C-stat",
        "d.o.f",
        "C-stat/d.o.f"
    ]

with open(csv_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=parameter_names)
    writer.writeheader()

Fit.statMethod = "cstat"
Fit.query = "yes"
Xset.abund = "lpgs"
Xset.xsect = "vern"

with open("flat_norm_factor.list", "r") as f:
    lines = f.readlines()
factor_lines = [float(line.strip()) for line in lines]

with open("region_area.list", "r") as f:
    lines = f.readlines()
region_area_lines = [float(line.strip()) for line in lines]

for i in range(start, finish + 1):

    AllData.clear()
    AllModels.clear()

    AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
    AllData.ignore("bad")
    s1 = AllData(1)
    s1.ignore("**-2.0 15.0-**")

    s1.multiresponse[89] = "xtend_nxb.rmf"
    s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

    s1.multiresponse[89].gain.slope = 1.0
    s1.multiresponse[89].gain.offset.frozen = True

    m1 = Model("TBabs*vapec", "src", 1)
    m1.TBabs.nH = nH_value
    m1.vapec.kT = 4.0
    m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
    m1.vapec.norm = 3.0

    m1.TBabs.nH.frozen = True
    m1.vapec.kT.frozen = False
    m1.vapec.Redshift.frozen = False
    m1.vapec.norm.frozen = False

    m1.vapec.He = 1.0
    m1.vapec.C = 1.0
    m1.vapec.N = 1.0
    m1.vapec.O = 1.0
    m1.vapec.Ne = 1.0
    m1.vapec.Mg = 1.0
    m1.vapec.Al = 1.0
    m1.vapec.Si = 1.0
    m1.vapec.S = 1.0
    m1.vapec.Ar = 1.0
    m1.vapec.Ca = 1.0
    m1.vapec.Fe = 1.0
    m1.vapec.Ni = 1.0

    m1.vapec.He.frozen = True
    m1.vapec.C.frozen = True
    m1.vapec.N.frozen = True
    m1.vapec.O.frozen = True
    m1.vapec.Ne.frozen = True
    m1.vapec.Mg.frozen = True
    m1.vapec.Al.frozen = True
    m1.vapec.Si.frozen = False
    m1.vapec.S.frozen = False
    m1.vapec.Ar.frozen = False
    m1.vapec.Ca.frozen = False
    m1.vapec.Fe.frozen = False
    m1.vapec.Ni.frozen = False

    Xset.restore("xtd_nxb.xcm")
    Xset.restore("xtd_cxb.xcm")

    m90 = AllModels(1, "det")
    m92 = AllModels(1, "cxb")

    def reverse_Normalize(norm_400pi):
        return norm_400pi * raytrace_source_area / (400*math.pi)

    m92.TBabs.nH = nH_value
    m92.powerlaw.norm = reverse_Normalize(CXB_norm_400pi)

    Fit.renorm('prefit')
    Fit.perform()

    if error_do_or_not == "yes":
        Fit.error("1.0 src:2 src:16 src:17")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.xcm"):
        os.remove(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.xcm")

    Xset.save(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.xcm", info='m')

    Xset.openLog(f"fit_flat_{i}_{model}_{version}.log")
    AllData.show()
    AllModels.show()
    Fit.show()
    Xset.closeLog()

    Plot.commands = ()
    Plot.device = "/xw"
    Plot.xAxis = "keV"
    Plot.xLog = True
    Plot.yLog = True
    Plot.setRebin(minSig=minSig_value, maxBins=maxBins_value)
    Plot.addCommand(f"r y {y_range_value} auto")
    Plot.addCommand(f'label top "Region {i} ({title})"')
    Plot("ldata delc")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.qdp"):
        os.remove(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.qdp")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.pco"):
        os.remove(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.pco")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.ps"):
        os.remove(f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.ps")

    Plot.device = f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.ps/cps"
    Plot("ldata delc")
    Plot.device = "/xw"
    Plot.addCommand(f"we xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.qdp")
    Plot("ldata delc")

    subprocess.run(["ps2pdf", f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.ps", f"xa{OBSID}xtd_{SUB}_flat_{i}_{model}_{version}.pdf"])

    def Normalize(norm):
        return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

    if error_do_or_not == "yes":
        csv_line = {
            "Region_number": i,
            "src_kT(keV)": m1.vapec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
            "src_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
            "src_Redshift": m1.vapec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.vapec.Redshift.error[0] == 0.0 else m1.vapec.Redshift.values[0] - m1.vapec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.vapec.Redshift.error[1] == 0.0 else m1.vapec.Redshift.error[1] - m1.vapec.Redshift.values[0],
            "src_Abund_He": m1.vapec.He.values[0],
            "src_Abund_C": m1.vapec.C.values[0],
            "src_Abund_N": m1.vapec.N.values[0],
            "src_Abund_O": m1.vapec.O.values[0],
            "src_Abund_Ne": m1.vapec.Ne.values[0],
            "src_Abund_Mg": m1.vapec.Mg.values[0],
            "src_Abund_Al": m1.vapec.Al.values[0],
            "src_Abund_Si": m1.vapec.Si.values[0],
            "src_Abund_S": m1.vapec.S.values[0],
            "src_Abund_Ar": m1.vapec.Ar.values[0],
            "src_Abund_Ca": m1.vapec.Ca.values[0],
            "src_Abund_Fe": m1.vapec.Fe.values[0],
            "src_Abund_Ni": m1.vapec.Ni.values[0],
            "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
            "C-stat": Fit.statistic,
            "d.o.f": Fit.dof,
            "C-stat/d.o.f": Fit.statistic/Fit.dof}

    elif error_do_or_not == "no":
        csv_line = {
            "Region_number": i,
            "src_kT(keV)": m1.vapec.kT.values[0],
            "src_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
            "src_Redshift": m1.vapec.Redshift.values[0],
            "src_Abund_He": m1.vapec.He.values[0],
            "src_Abund_C": m1.vapec.C.values[0],
            "src_Abund_N": m1.vapec.N.values[0],
            "src_Abund_O": m1.vapec.O.values[0],
            "src_Abund_Ne": m1.vapec.Ne.values[0],
            "src_Abund_Mg": m1.vapec.Mg.values[0],
            "src_Abund_Al": m1.vapec.Al.values[0],
            "src_Abund_Si": m1.vapec.Si.values[0],
            "src_Abund_S": m1.vapec.S.values[0],
            "src_Abund_Ar": m1.vapec.Ar.values[0],
            "src_Abund_Ca": m1.vapec.Ca.values[0],
            "src_Abund_Fe": m1.vapec.Fe.values[0],
            "src_Abund_Ni": m1.vapec.Ni.values[0],
            "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
            "C-stat": Fit.statistic,
            "d.o.f": Fit.dof,
            "C-stat/d.o.f": Fit.statistic/Fit.dof}

    with open(csv_file, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writerow(csv_line)

for log_file in glob.glob("*.log"):
    shutil.copy(log_file, f"log/{log_file}")
    os.remove(log_file)

print("\n All processes have been completed.\n")