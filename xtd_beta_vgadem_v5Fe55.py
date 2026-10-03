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
model = "vgadem"
version = "v5Fe55"

OBSID = os.environ.get("OBSID")
SUB = os.environ.get("SUB")

start = int(input("\n Input the starting region number. > ").strip())
finish = int(input(" Input the finishing region number. > ").strip())
error_do_or_not = input(" Execute the error command? (yes|no) > ").strip()

if error_do_or_not not in ["yes", "no"]:
    print(" Please input 'yes' or 'no'.\n")
    sys.exit(1)

if error_do_or_not == "yes":
    csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
        "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
        "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
        "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
        "src_Abund_He",
        "src_Abund_C",
        "src_Abund_N",
        "src_Abund_O",
        "src_Abund_Ne",
        "src_Abund_Na",
        "src_Abund_Mg",
        "src_Abund_Al",
        "src_Abund_Si",
        "src_Abund_S",
        "src_Abund_Ar",
        "src_Abund_Ca",
        "src_Abund_Fe",
        "src_Abund_Ni",
        "NXB_gain_slope",
        "Fe55_5.89keV_norm(/400pi_arcmin2)",
        "Fe55_6.49keV_norm(/400pi_arcmin2)",
        "C-stat",
        "d.o.f",
        "C-stat/d.o.f"
    ]

elif error_do_or_not == "no":
    csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
    parameter_names = [
        "Region_number",
        "src_Tmean(keV)",
        "src_Tsigma(keV)",
        "src_Redshift",
        "src_norm(/400pi_arcmin2)",
        "src_Abund_He",
        "src_Abund_C",
        "src_Abund_N",
        "src_Abund_O",
        "src_Abund_Ne",
        "src_Abund_Na",
        "src_Abund_Mg",
        "src_Abund_Al",
        "src_Abund_Si",
        "src_Abund_S",
        "src_Abund_Ar",
        "src_Abund_Ca",
        "src_Abund_Fe",
        "src_Abund_Ni",
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

with open("beta_norm_factor.list", "r") as f:
    lines = f.readlines()
factor_lines = [float(line.strip()) for line in lines]

with open("region_area.list", "r") as f:
    lines = f.readlines()
region_area_lines = [float(line.strip()) for line in lines]

for i in range(start, finish + 1):

    AllData.clear()
    AllModels.clear()

    AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
    AllData.ignore("bad")
    s1 = AllData(1)
    s1.ignore("**-2.0 15.0-**")

    s1.multiresponse[89] = "xtend_nxb.rmf"
    s1.multiresponse[90] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[90].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"
    s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
    s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

    s1.multiresponse[89].gain.slope = 1.0
    s1.multiresponse[89].gain.offset.frozen = True
    s1.multiresponse[90].gain.slope = "1.0,,0.98,0.98,1.02,1.02"
    s1.multiresponse[90].gain.offset.frozen = True

    m1 = Model("TBabs*vgadem", "src", 1)
    m1.TBabs.nH = nH_value
    m1.vgadem.Tmean = 4.0
    m1.vgadem.Tsigma = 1.0
    m1.vgadem.nH = 1.0
    m1.vgadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
    m1.vgadem.switch = 2
    m1.vgadem.norm = 1.0

    m1.TBabs.nH.frozen = True
    m1.vgadem.Tmean.frozen = False
    m1.vgadem.Tsigma.frozen = False
    m1.vgadem.nH.frozen = True
    m1.vgadem.Redshift.frozen = False
    m1.vgadem.switch.frozen = True
    m1.vgadem.norm.frozen = False

    m1.vgadem.He = 1.0
    m1.vgadem.C = 1.0
    m1.vgadem.N = 1.0
    m1.vgadem.O = 1.0
    m1.vgadem.Ne = 1.0
    m1.vgadem.Na = 1.0
    m1.vgadem.Mg = 1.0
    m1.vgadem.Al = 1.0
    m1.vgadem.Si = 1.0
    m1.vgadem.S = 1.0
    m1.vgadem.Ar = 1.0
    m1.vgadem.Ca = 1.0
    m1.vgadem.Fe = 1.0
    m1.vgadem.Ni = 1.0

    m1.vgadem.He.frozen = True
    m1.vgadem.C.frozen = True
    m1.vgadem.N.frozen = True
    m1.vgadem.O.frozen = True
    m1.vgadem.Ne.frozen = True
    m1.vgadem.Na.frozen = True
    m1.vgadem.Mg.frozen = True
    m1.vgadem.Al.frozen = True
    m1.vgadem.Si.frozen = False
    m1.vgadem.S.frozen = False
    m1.vgadem.Ar.frozen = False
    m1.vgadem.Ca.frozen = False
    m1.vgadem.Fe.frozen = False
    m1.vgadem.Ni.frozen = False

    Xset.restore("xtd_nxb.xcm")
    Xset.restore("xtd_Fe55.xcm")
    Xset.restore("xtd_cxb.xcm")

    m90 = AllModels(1, "det")
    m91 = AllModels(1, "Fe55")
    m92 = AllModels(1, "cxb")

    def reverse_Normalize(norm_400pi):
        return norm_400pi * raytrace_source_area / (400*math.pi)

    m92.TBabs.nH = nH_value
    m92.powerlaw.norm = reverse_Normalize(CXB_norm_400pi)

    Fit.renorm('prefit')
    Fit.perform()

    if error_do_or_not == "yes":
        Fit.error("1.0 src:2 src:3 src:19 src:21")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.xcm"):
        os.remove(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.xcm")

    Xset.save(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.xcm", info='m')

    Xset.openLog(f"fit_beta_{i}_{model}_{version}.log")
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

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.qdp"):
        os.remove(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.qdp")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.pco"):
        os.remove(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.pco")

    if os.path.exists(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.ps"):
        os.remove(f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.ps")

    Plot.device = f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.ps/cps"
    Plot("ldata delc")
    Plot.device = "/xw"
    Plot.addCommand(f"we xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.qdp")
    Plot("ldata delc")

    subprocess.run(["ps2pdf", f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.ps", f"xa{OBSID}xtd_{SUB}_beta_{i}_{model}_{version}.pdf"])

    def Normalize(norm):
        return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

    if error_do_or_not == "yes":
        csv_line = {
            "Region_number": i,
            "src_Tmean(keV)": m1.vgadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.vgadem.Tmean.error[0] == 0.0 else m1.vgadem.Tmean.values[0] - m1.vgadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.vgadem.Tmean.error[1] == 0.0 else m1.vgadem.Tmean.error[1] - m1.vgadem.Tmean.values[0],
            "src_Tsigma(keV)": m1.vgadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.vgadem.Tsigma.error[0] == 0.0 else m1.vgadem.Tsigma.values[0] - m1.vgadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.vgadem.Tsigma.error[1] == 0.0 else m1.vgadem.Tsigma.error[1] - m1.vgadem.Tsigma.values[0],
            "src_Redshift": m1.vgadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.vgadem.Redshift.error[0] == 0.0 else m1.vgadem.Redshift.values[0] - m1.vgadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.vgadem.Redshift.error[1] == 0.0 else m1.vgadem.Redshift.error[1] - m1.vgadem.Redshift.values[0],
            "src_norm(/400pi_arcmin2)": Normalize(m1.vgadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.vgadem.norm.error[0] == 0.0 else Normalize(m1.vgadem.norm.values[0] - m1.vgadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.vgadem.norm.error[1] == 0.0 else Normalize(m1.vgadem.norm.error[1] - m1.vgadem.norm.values[0]),
            "src_Abund_He": m1.vgadem.He.values[0],
            "src_Abund_C": m1.vgadem.C.values[0],
            "src_Abund_N": m1.vgadem.N.values[0],
            "src_Abund_O": m1.vgadem.O.values[0],
            "src_Abund_Ne": m1.vgadem.Ne.values[0],
            "src_Abund_Na": m1.vgadem.Na.values[0],
            "src_Abund_Mg": m1.vgadem.Mg.values[0],
            "src_Abund_Al": m1.vgadem.Al.values[0],
            "src_Abund_Si": m1.vgadem.Si.values[0],
            "src_Abund_S": m1.vgadem.S.values[0],
            "src_Abund_Ar": m1.vgadem.Ar.values[0],
            "src_Abund_Ca": m1.vgadem.Ca.values[0],
            "src_Abund_Fe": m1.vgadem.Fe.values[0],
            "src_Abund_Ni": m1.vgadem.Ni.values[0],
            "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
            "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
            "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
            "C-stat": Fit.statistic,
            "d.o.f": Fit.dof,
            "C-stat/d.o.f": Fit.statistic/Fit.dof}

    elif error_do_or_not == "no":
        csv_line = {
            "Region_number": i,
            "src_Tmean(keV)": m1.vgadem.Tmean.values[0],
            "src_Tsigma(keV)": m1.vgadem.Tsigma.values[0],
            "src_Redshift": m1.vgadem.Redshift.values[0],
            "src_norm(/400pi_arcmin2)": Normalize(m1.vgadem.norm.values[0]),
            "src_Abund_He": m1.vgadem.He.values[0],
            "src_Abund_C": m1.vgadem.C.values[0],
            "src_Abund_N": m1.vgadem.N.values[0],
            "src_Abund_O": m1.vgadem.O.values[0],
            "src_Abund_Ne": m1.vgadem.Ne.values[0],
            "src_Abund_Na": m1.vgadem.Na.values[0],
            "src_Abund_Mg": m1.vgadem.Mg.values[0],
            "src_Abund_Al": m1.vgadem.Al.values[0],
            "src_Abund_Si": m1.vgadem.Si.values[0],
            "src_Abund_S": m1.vgadem.S.values[0],
            "src_Abund_Ar": m1.vgadem.Ar.values[0],
            "src_Abund_Ca": m1.vgadem.Ca.values[0],
            "src_Abund_Fe": m1.vgadem.Fe.values[0],
            "src_Abund_Ni": m1.vgadem.Ni.values[0],
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