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
model = "2apec"
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
        "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
        "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
        "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
        "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
        "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
        "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
        "NXB_gain_slope",
        "C-stat",
        "d.o.f",
        "C-stat/d.o.f"
    ]

elif error_do_or_not == "no":
    csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "src1_kT(keV)",
        "src1_norm(/400pi_arcmin2)",
        "src2_kT(keV)",
        "src2_norm(/400pi_arcmin2)",
        "src_Abundance",
        "src_Redshift",
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
Plot.add = True

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

    m1 = Model("TBabs*(apec+apec)", "src", 1)
    m1.TBabs.nH = nH_value
    m1.apec.kT = 4.0
    m1.apec.Abundanc = 1.0
    m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
    m1.apec.norm = 3.0
    m1.apec_3.kT = 2.0
    m1.apec_3.Abundanc.link = m1.apec.Abundanc
    m1.apec_3.Redshift.link = m1.apec.Redshift
    m1.apec_3.norm = 3.0

    m1.TBabs.nH.frozen = True
    m1.apec.kT.frozen = False
    m1.apec.Abundanc.frozen = False
    m1.apec.Redshift.frozen = False
    m1.apec.norm.frozen = False
    m1.apec_3.kT.frozen = False
    m1.apec_3.Abundanc.frozen = False
    m1.apec_3.Redshift.frozen = False
    m1.apec_3.norm.frozen = False

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
        Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9")

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
            "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
            "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
            "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
            "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
            "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
            "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
            "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
            "C-stat": Fit.statistic,
            "d.o.f": Fit.dof,
            "C-stat/d.o.f": Fit.statistic/Fit.dof}

    elif error_do_or_not == "no":
        csv_line = {
            "Region_number": i,
            "src1_kT(keV)": m1.apec.kT.values[0],
            "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
            "src2_kT(keV)": m1.apec_3.kT.values[0],
            "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]),
            "src_Abundance": m1.apec.Abundanc.values[0],
            "src_Redshift": m1.apec.Redshift.values[0],
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