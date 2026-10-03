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
GH_norm_400pi = 1.0
GH_norm_low = 1.0
GH_norm_high = 1.0
LHB_norm_400pi = 1.0
LHB_norm_low = 1.0
LHB_norm_high = 1.0
keV_norm_400pi = 1.0
keV_norm_low = 1.0
keV_norm_high = 1.0
OI_norm_400pi = 1.0
OI_norm_low = 1.0
OI_norm_high = 1.0
minSig_value = 5
maxBins_value = 10
y_range_value = 1.0E-5
title = "None"
model = "2apec"
version = "v4Fe55"

OBSID = os.environ.get("OBSID")
SUB = os.environ.get("SUB")

start = int(input("\n Input the starting region number. > ").strip())
finish = int(input(" Input the finishing region number. > ").strip())
error_do_or_not = input(" Execute the error command? (yes|no) > ").strip()

if error_do_or_not not in ["yes", "no"]:
    print(" Please input 'yes' or 'no'.\n")
    sys.exit(1)

if error_do_or_not == "yes":
    csv_file = f"{start}_to_{finish}_with_error_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "ARF_radius(arcmin)",
        "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
        "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
        "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
        "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
        "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
        "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
        "NXB_gain_slope",
        "Fe55_5.89keV_norm(/400pi_arcmin2)",
        "Fe55_6.49keV_norm(/400pi_arcmin2)",
        "C-stat",
        "d.o.f",
        "C-stat/d.o.f"
    ]

elif error_do_or_not == "no":
    csv_file = f"{start}_to_{finish}_without_error_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "ARF_radius(arcmin)",
        "src1_kT(keV)",
        "src1_norm(/400pi_arcmin2)",
        "src2_kT(keV)",
        "src2_norm(/400pi_arcmin2)",
        "src_Abundance",
        "src_Redshift",
        "NXB_gain_slope",
        "Fe55_5.89keV_norm(/400pi_arcmin2)",
        "Fe55_6.49keV_norm(/400pi_arcmin2)",
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

with open("radius.list", "r") as f:
    lines = f.readlines()
arf_radius_lines = [float(line.strip()) for line in lines]

for i in range(start, finish + 1):

    AllData.clear()
    AllModels.clear()

    AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
    AllData.ignore("bad")
    s1 = AllData(1)
    s1.ignore("**-0.5 15.0-**")

    s1.multiresponse[89] = "xtend_nxb.rmf"
    s1.multiresponse[90] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[90].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"
    s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"
    s1.multiresponse[92] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[92].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"
    s1.multiresponse[93] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[93].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"
    s1.multiresponse[94] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[94].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"
    s1.multiresponse[95] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[95].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

    s1.multiresponse[89].gain.slope = 1.0
    s1.multiresponse[89].gain.offset.frozen = True
    s1.multiresponse[90].gain.slope = "1.0,,0.98,0.98,1.02,1.02"
    s1.multiresponse[90].gain.offset.frozen = True

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
    Xset.restore("xtd_Fe55.xcm")
    Xset.restore("xtd_cxb.xcm")
    Xset.restore("xtd_gh.xcm")
    Xset.restore("xtd_lhb.xcm")
    Xset.restore("xtd_0.8keV.xcm")
    Xset.restore("xtd_OI.xcm")


    m90 = AllModels(1, "det")
    m91 = AllModels(1, "Fe55")
    m92 = AllModels(1, "cxb")
    m93 = AllModels(1, "gh")
    m94 = AllModels(1, "lhb")
    m95 = AllModels(1, "0.8keV")
    m96 = AllModels(1, "OI")

    def reverse_Normalize(norm_400pi):
        return norm_400pi * (arf_radius_lines[i-1]*arf_radius_lines[i-1]*math.pi) / (400*math.pi)

    m92.TBabs.nH = nH_value
    m92.powerlaw.norm = reverse_Normalize(CXB_norm_400pi)
    m93.TBabs.nH = nH_value
    m93.apec.norm = [reverse_Normalize(GH_norm_400pi), 0.01, reverse_Normalize(GH_norm_low), reverse_Normalize(GH_norm_low), reverse_Normalize(GH_norm_high), reverse_Normalize(GH_norm_high)]
    m94.apec.norm = [reverse_Normalize(LHB_norm_400pi), 0.01, reverse_Normalize(LHB_norm_low), reverse_Normalize(LHB_norm_low), reverse_Normalize(LHB_norm_high), reverse_Normalize(LHB_norm_high)]
    m95.TBabs.nH = nH_value
    m95.apec.norm = [reverse_Normalize(keV_norm_400pi), 0.01, reverse_Normalize(keV_norm_low), reverse_Normalize(keV_norm_low), reverse_Normalize(keV_norm_high), reverse_Normalize(keV_norm_high)]
    m96.gaussian.norm = [reverse_Normalize(OI_norm_400pi), 0.01, reverse_Normalize(OI_norm_low), reverse_Normalize(OI_norm_low), reverse_Normalize(OI_norm_high), reverse_Normalize(OI_norm_high)]

    Fit.renorm('prefit')
    Fit.perform()

    if error_do_or_not == "yes":
        Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.xcm"):
        os.remove(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.xcm")

    Xset.save(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.xcm", info='m')

    Xset.openLog(f"fit_{i}_{model}_{version}.log")
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

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.qdp"):
        os.remove(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.qdp")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.pco"):
        os.remove(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.pco")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.ps"):
        os.remove(f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.ps")

    Plot.device = f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.ps/cps"
    Plot("ldata delc")
    Plot.device = "/xw"
    Plot.addCommand(f"we xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.qdp")
    Plot("ldata delc")

    subprocess.run(["ps2pdf", f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.ps", f"xa{OBSID}xtd_{SUB}_src_{i}_{model}_{version}.pdf"])

    def Normalize(norm):
        return norm * 400*math.pi / (arf_radius_lines[i-1]*arf_radius_lines[i-1]*math.pi)

    if error_do_or_not == "yes":
        csv_line = {
            "Region_number": i,
            "ARF_radius(arcmin)": arf_radius_lines[i-1],
            "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
            "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
            "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
            "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
            "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
            "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
            "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
            "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
            "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
            "C-stat": Fit.statistic,
            "d.o.f": Fit.dof,
            "C-stat/d.o.f": Fit.statistic/Fit.dof}

    elif error_do_or_not == "no":
        csv_line = {
            "Region_number": i,
            "ARF_radius(arcmin)": arf_radius_lines[i-1],
            "src1_kT(keV)": m1.apec.kT.values[0],
            "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
            "src2_kT(keV)": m1.apec_3.kT.values[0],
            "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]),
            "src_Abundance": m1.apec.Abundanc.values[0],
            "src_Redshift": m1.apec.Redshift.values[0],
            "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
            "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
            "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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