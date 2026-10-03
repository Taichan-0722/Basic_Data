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
GH_norm_400pi = 1.5E-3
GH_norm_low = 0
GH_norm_high = 1.0E-2
LHB_norm_400pi = 1.6E-2
LHB_norm_low = 0
LHB_norm_high = 1.9E-2
keV_norm_400pi = 1.0E-4
keV_norm_low = 0
keV_norm_high = 1.0E-3
OI_norm_400pi = 5E-4
OI_norm_low = 0
OI_norm_high = 1.0E-3

minSig_value = 5
maxBins_value = 10
y_range_value = 1.0E-5
raytrace_source_area = 44 * 44   # arcmin^2

title = "TITLE"
arf_type = "flat"
model = "1apec"
version = "v4"

Fit.statMethod = "cstat"
Fit.query = "yes"
Xset.abund = "lpgs"
Xset.xsect = "vern"
Plot.add = False

OBSID = os.environ.get("OBSID")
SUB = os.environ.get("SUB")

start = int(input("\n Input the starting region number. > ").strip())
finish = int(input(" Input the finishing region number. > ").strip())
error_do_or_not = input(" Execute the error command? (yes|no) > ").strip()

if error_do_or_not not in ["yes", "no"]:
    print(" Please input 'yes' or 'no'.\n")
    sys.exit(1)

with open("norm_factor.list", "r") as f:
    lines = f.readlines()
factor_lines = [float(line.strip()) for line in lines]

with open("region_area.list", "r") as f:
    lines = f.readlines()
region_area_lines = [float(line.strip()) for line in lines]


if arf_type == "flat" and model == "1apec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
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
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "1apec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "1apec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
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
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "1apec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "1vapec" and version == "v4":

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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17")

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


if arf_type == "flat" and model == "1vapec" and version == "v4Fe55":

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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17")

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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "1vapec" and version == "v5":

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
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


if arf_type == "flat" and model == "1vapec" and version == "v5Fe55":

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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

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
            Fit.error("1.0 src:2 src:16 src:17")

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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "2apec" and version == "v4":

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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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


if arf_type == "flat" and model == "2apec" and version == "v4Fe55":

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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "2apec" and version == "v5":

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


if arf_type == "flat" and model == "2apec" and version == "v5Fe55":

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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "2vapec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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


if arf_type == "flat" and model == "2vapec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "2vapec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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


if arf_type == "flat" and model == "2vapec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "3apec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "3apec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
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


if arf_type == "flat" and model == "3apec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "3apec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
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


if arf_type == "flat" and model == "3vapec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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


if arf_type == "flat" and model == "3vapec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "3vapec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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


if arf_type == "flat" and model == "3vapec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
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
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "gadem" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)",
            "src_Tsigma(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "gadem" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)",
            "src_Tsigma(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "gadem" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)",
            "src_Tsigma(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "gadem" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)",
            "src_Tsigma(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "vgadem" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
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
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
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
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
        m1.vgadem.Al.frozen = True
        m1.vgadem.Si.frozen = False
        m1.vgadem.S.frozen = False
        m1.vgadem.Ar.frozen = False
        m1.vgadem.Ca.frozen = False
        m1.vgadem.Fe.frozen = False
        m1.vgadem.Ni.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:19 src:21")

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
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "vgadem" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
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
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:19 src:21")

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


if arf_type == "flat" and model == "vgadem" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
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
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
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
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
        m1.vgadem.Al.frozen = True
        m1.vgadem.Si.frozen = False
        m1.vgadem.S.frozen = False
        m1.vgadem.Ar.frozen = False
        m1.vgadem.Ca.frozen = False
        m1.vgadem.Fe.frozen = False
        m1.vgadem.Ni.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:19 src:21")

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
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "flat" and model == "vgadem" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_flat_{model}_{version}.csv"
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
        csv_file = f"{start}_to_{finish}_without_error_flat_{model}_{version}.csv"
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_flat_{i}.grp.gz")
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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
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


if arf_type == "beta" and model == "1apec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "1apec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
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
            "src_kT(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "1apec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-2.0 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
        s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
        s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

        s1.multiresponse[89].gain.slope = 1.0
        s1.multiresponse[89].gain.offset.frozen = True

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "1apec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_kT(keV)", "src_kT_error_bottom", "src_kT_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
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
            "src_kT(keV)",
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*apec", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0], "src_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.apec.Abundanc.error[0] == 0.0 else m1.apec.Abundanc.values[0] - m1.apec.Abundanc.error[0], "src_Abund_error_top": 0.0 if m1.apec.Abundanc.error[1] == 0.0 else m1.apec.Abundanc.error[1] - m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.apec.Redshift.error[0] == 0.0 else m1.apec.Redshift.values[0] - m1.apec.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.apec.Redshift.error[1] == 0.0 else m1.apec.Redshift.error[1] - m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_kT(keV)": m1.apec.kT.values[0],
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "1vapec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17")

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


if arf_type == "beta" and model == "1vapec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17")

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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "1vapec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
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


if arf_type == "beta" and model == "1vapec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

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
            Fit.error("1.0 src:2 src:16 src:17")

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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "2apec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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


if arf_type == "beta" and model == "2apec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "2apec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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


if arf_type == "beta" and model == "2apec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "2vapec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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


if arf_type == "beta" and model == "2vapec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "2vapec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-2.0 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
        s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
        s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

        s1.multiresponse[89].gain.slope = 1.0
        s1.multiresponse[89].gain.offset.frozen = True

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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


if arf_type == "beta" and model == "2vapec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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

        m1 = Model("TBabs*(vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 3.0
        m1.vapec_3.kT = 2.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "3apec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "3apec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
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


if arf_type == "beta" and model == "3apec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    elif error_do_or_not == "no":
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-2.0 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
        s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
        s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

        s1.multiresponse[89].gain.slope = 1.0
        s1.multiresponse[89].gain.offset.frozen = True

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
                "src_Abundance": m1.apec.Abundanc.values[0],
                "src_Redshift": m1.apec.Redshift.values[0],
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "3apec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*(apec+apec+apec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.apec.kT = 4.0
        m1.apec.Abundanc = 1.0
        m1.apec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.apec.norm = 1.0
        m1.apec_3.kT = 3.0
        m1.apec_3.Abundanc.link = m1.apec.Abundanc
        m1.apec_3.Redshift.link = m1.apec.Redshift
        m1.apec_3.norm = 1.0
        m1.apec_4.kT = 2.0
        m1.apec_4.Abundanc.link = m1.apec.Abundanc
        m1.apec_4.Redshift.link = m1.apec.Redshift
        m1.apec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.apec.kT.frozen = False
        m1.apec.Abundanc.frozen = False
        m1.apec.Redshift.frozen = False
        m1.apec.norm.frozen = False
        m1.apec_3.kT.frozen = False
        m1.apec_3.Abundanc.frozen = False
        m1.apec_3.Redshift.frozen = False
        m1.apec_3.norm.frozen = False
        m1.apec_4.kT.frozen = False
        m1.apec_4.Abundanc.frozen = False
        m1.apec_4.Redshift.frozen = False
        m1.apec_4.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:4 src:5 src:6 src:9 src:10 src:13")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.apec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.apec.kT.error[0] == 0.0 else m1.apec.kT.values[0] - m1.apec.kT.error[0], "src1_kT_error_top": 0.0 if m1.apec.kT.error[1] == 0.0 else m1.apec.kT.error[1] - m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.apec.norm.error[0] == 0.0 else Normalize(m1.apec.norm.values[0] - m1.apec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.apec.norm.error[1] == 0.0 else Normalize(m1.apec.norm.error[1] - m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.apec_3.kT.error[0] == 0.0 else m1.apec_3.kT.values[0] - m1.apec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.apec_3.kT.error[1] == 0.0 else m1.apec_3.kT.error[1] - m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.apec_3.norm.error[0] == 0.0 else Normalize(m1.apec_3.norm.values[0] - m1.apec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.apec_3.norm.error[1] == 0.0 else Normalize(m1.apec_3.norm.error[1] - m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.apec_4.kT.error[0] == 0.0 else m1.apec_4.kT.values[0] - m1.apec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.apec_4.kT.error[1] == 0.0 else m1.apec_4.kT.error[1] - m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.apec_4.norm.error[0] == 0.0 else Normalize(m1.apec_4.norm.values[0] - m1.apec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.apec_4.norm.error[1] == 0.0 else Normalize(m1.apec_4.norm.error[1] - m1.apec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.apec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.apec.norm.values[0]),
                "src2_kT(keV)": m1.apec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.apec_3.norm.values[0]),
                "src3_kT(keV)": m1.apec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.apec_4.norm.values[0]),
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


if arf_type == "beta" and model == "3vapec" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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


if arf_type == "beta" and model == "3vapec" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "3vapec" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
        csv_file = f"{start}_to_{finish}_without_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-2.0 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
        s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
        s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

        s1.multiresponse[89].gain.slope = 1.0
        s1.multiresponse[89].gain.offset.frozen = True

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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


if arf_type == "beta" and model == "3vapec" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src1_kT(keV)", "src1_kT_error_bottom", "src1_kT_error_top",
            "src1_norm(/400pi_arcmin2)", "src1_norm_error_bottom", "src1_norm_error_top",
            "src2_kT(keV)", "src2_kT_error_bottom", "src2_kT_error_top",
            "src2_norm(/400pi_arcmin2)", "src2_norm_error_bottom", "src2_norm_error_top",
            "src3_kT(keV)", "src3_kT_error_bottom", "src3_kT_error_top",
            "src3_norm(/400pi_arcmin2)", "src3_norm_error_bottom", "src3_norm_error_top",
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
            "src1_kT(keV)",
            "src1_norm(/400pi_arcmin2)",
            "src2_kT(keV)",
            "src2_norm(/400pi_arcmin2)",
            "src3_kT(keV)",
            "src3_norm(/400pi_arcmin2)",
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
            "Fe55_5.89keV_norm(/400pi_arcmin2)",
            "Fe55_6.49keV_norm(/400pi_arcmin2)",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

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

        m1 = Model("TBabs*(vapec+vapec+vapec)", "src", 1)
        m1.TBabs.nH = nH_value
        m1.vapec.kT = 4.0
        m1.vapec.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.vapec.norm = 1.0
        m1.vapec_3.kT = 3.0
        m1.vapec_3.Redshift.link = m1.vapec.Redshift
        m1.vapec_3.norm = 1.0
        m1.vapec_4.kT = 2.0
        m1.vapec_4.Redshift.link = m1.vapec.Redshift
        m1.vapec_4.norm = 1.0

        m1.TBabs.nH.frozen = True
        m1.vapec.kT.frozen = False
        m1.vapec.Redshift.frozen = False
        m1.vapec.norm.frozen = False
        m1.vapec_3.kT.frozen = False
        m1.vapec_3.Redshift.frozen = False
        m1.vapec_3.norm.frozen = False
        m1.vapec_4.kT.frozen = False
        m1.vapec_4.Redshift.frozen = False
        m1.vapec_4.norm.frozen = False

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
        m1.vapec.O.frozen = False
        m1.vapec.Ne.frozen = False
        m1.vapec.Mg.frozen = False
        m1.vapec.Al.frozen = True
        m1.vapec.Si.frozen = False
        m1.vapec.S.frozen = False
        m1.vapec.Ar.frozen = False
        m1.vapec.Ca.frozen = False
        m1.vapec.Fe.frozen = False
        m1.vapec.Ni.frozen = False

        m1.vapec_3.He.link = m1.vapec_4.He.link = m1.vapec.He
        m1.vapec_3.C.link = m1.vapec_4.C.link = m1.vapec.C
        m1.vapec_3.N.link = m1.vapec_4.N.link = m1.vapec.N
        m1.vapec_3.O.link = m1.vapec_4.O.link = m1.vapec.O
        m1.vapec_3.Ne.link = m1.vapec_4.Ne.link = m1.vapec.Ne
        m1.vapec_3.Mg.link = m1.vapec_4.Mg.link = m1.vapec.Mg
        m1.vapec_3.Al.link = m1.vapec_4.Al.link = m1.vapec.Al
        m1.vapec_3.Si.link = m1.vapec_4.Si.link = m1.vapec.Si
        m1.vapec_3.S.link = m1.vapec_4.S.link = m1.vapec.S
        m1.vapec_3.Ar.link = m1.vapec_4.Ar.link = m1.vapec.Ar
        m1.vapec_3.Ca.link = m1.vapec_4.Ca.link = m1.vapec.Ca
        m1.vapec_3.Fe.link = m1.vapec_4.Fe.link = m1.vapec.Fe
        m1.vapec_3.Ni.link = m1.vapec_4.Ni.link = m1.vapec.Ni

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
            Fit.error("1.0 src:2 src:16 src:17 src:18 src:33 src:34 src:49")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0], "src1_kT_error_bottom": 0.0 if m1.vapec.kT.error[0] == 0.0 else m1.vapec.kT.values[0] - m1.vapec.kT.error[0], "src1_kT_error_top": 0.0 if m1.vapec.kT.error[1] == 0.0 else m1.vapec.kT.error[1] - m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]), "src1_norm_error_bottom": 0.0 if m1.vapec.norm.error[0] == 0.0 else Normalize(m1.vapec.norm.values[0] - m1.vapec.norm.error[0]), "src1_norm_error_top": 0.0 if m1.vapec.norm.error[1] == 0.0 else Normalize(m1.vapec.norm.error[1] - m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0], "src2_kT_error_bottom": 0.0 if m1.vapec_3.kT.error[0] == 0.0 else m1.vapec_3.kT.values[0] - m1.vapec_3.kT.error[0], "src2_kT_error_top": 0.0 if m1.vapec_3.kT.error[1] == 0.0 else m1.vapec_3.kT.error[1] - m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]), "src2_norm_error_bottom": 0.0 if m1.vapec_3.norm.error[0] == 0.0 else Normalize(m1.vapec_3.norm.values[0] - m1.vapec_3.norm.error[0]), "src2_norm_error_top": 0.0 if m1.vapec_3.norm.error[1] == 0.0 else Normalize(m1.vapec_3.norm.error[1] - m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0], "src3_kT_error_bottom": 0.0 if m1.vapec_4.kT.error[0] == 0.0 else m1.vapec_4.kT.values[0] - m1.vapec_4.kT.error[0], "src3_kT_error_top": 0.0 if m1.vapec_4.kT.error[1] == 0.0 else m1.vapec_4.kT.error[1] - m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]), "src3_norm_error_bottom": 0.0 if m1.vapec_4.norm.error[0] == 0.0 else Normalize(m1.vapec_4.norm.values[0] - m1.vapec_4.norm.error[0]), "src3_norm_error_top": 0.0 if m1.vapec_4.norm.error[1] == 0.0 else Normalize(m1.vapec_4.norm.error[1] - m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src1_kT(keV)": m1.vapec.kT.values[0],
                "src1_norm(/400pi_arcmin2)": Normalize(m1.vapec.norm.values[0]),
                "src2_kT(keV)": m1.vapec_3.kT.values[0],
                "src2_norm(/400pi_arcmin2)": Normalize(m1.vapec_3.norm.values[0]),
                "src3_kT(keV)": m1.vapec_4.kT.values[0],
                "src3_norm(/400pi_arcmin2)": Normalize(m1.vapec_4.norm.values[0]),
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
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "gadem" and version == "v4":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
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
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "gadem" and version == "v4Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
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
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "gadem" and version == "v5":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
            "NXB_gain_slope",
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
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
            "NXB_gain_slope",
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-2.0 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
        s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
        s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

        s1.multiresponse[89].gain.slope = 1.0
        s1.multiresponse[89].gain.offset.frozen = True

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "gadem" and version == "v5Fe55":

    if error_do_or_not == "yes":
        csv_file = f"{start}_to_{finish}_with_error_beta_{model}_{version}.csv"
        parameter_names = [
            "Region_number",
            "src_Tmean(keV)", "src_Tmean_error_bottom", "src_Tmean_error_top",
            "src_Tsigma(keV)", "src_Tsigma_error_bottom", "src_Tsigma_error_top",
            "src_Abundance", "src_Abund_error_bottom", "src_Abund_error_top",
            "src_Redshift", "src_Redshift_error_bottom", "src_Redshift_error_top",
            "src_norm(/400pi_arcmin2)", "src_norm_error_bottom", "src_norm_error_top",
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
            "src_Abundance",
            "src_Redshift",
            "src_norm(/400pi_arcmin2)",
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

        m1 = Model("TBabs*gadem", "src", 1)
        m1.TBabs.nH = nH_value
        m1.gadem.Tmean = 4.0
        m1.gadem.Tsigma = 1.0
        m1.gadem.nH = 1.0
        m1.gadem.abundanc = 1.0
        m1.gadem.Redshift = [Redshift_value, 0.01, Redshift_low, Redshift_low, Redshift_high, Redshift_high]
        m1.gadem.switch = 2
        m1.gadem.norm = 3.0

        m1.TBabs.nH.frozen = True
        m1.gadem.Tmean.frozen = False
        m1.gadem.Tsigma.frozen = False
        m1.gadem.nH.frozen = True
        m1.gadem.abundanc.frozen = False
        m1.gadem.Redshift.frozen = False
        m1.gadem.switch.frozen = True
        m1.gadem.norm.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:5 src:6 src:8")

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
            return norm * factor_lines[i-1] * 400*math.pi / (region_area_lines[i-1])

        if error_do_or_not == "yes":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0], "src_Tmean_error_bottom": 0.0 if m1.gadem.Tmean.error[0] == 0.0 else m1.gadem.Tmean.values[0] - m1.gadem.Tmean.error[0], "src_Tmean_error_top": 0.0 if m1.gadem.Tmean.error[1] == 0.0 else m1.gadem.Tmean.error[1] - m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0], "src_Tsigma_error_bottom": 0.0 if m1.gadem.Tsigma.error[0] == 0.0 else m1.gadem.Tsigma.values[0] - m1.gadem.Tsigma.error[0], "src_Tsigma_error_top": 0.0 if m1.gadem.Tsigma.error[1] == 0.0 else m1.gadem.Tsigma.error[1] - m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0], "src_Abund_error_bottom": 0.0 if m1.gadem.abundanc.error[0] == 0.0 else m1.gadem.abundanc.values[0] - m1.gadem.abundanc.error[0], "src_Abund_error_top": 0.0 if m1.gadem.abundanc.error[1] == 0.0 else m1.gadem.abundanc.error[1] - m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0], "src_Redshift_error_bottom": 0.0 if m1.gadem.Redshift.error[0] == 0.0 else m1.gadem.Redshift.values[0] - m1.gadem.Redshift.error[0], "src_Redshift_error_top": 0.0 if m1.gadem.Redshift.error[1] == 0.0 else m1.gadem.Redshift.error[1] - m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]), "src_norm_error_bottom": 0.0 if m1.gadem.norm.error[0] == 0.0 else Normalize(m1.gadem.norm.values[0] - m1.gadem.norm.error[0]), "src_norm_error_top": 0.0 if m1.gadem.norm.error[1] == 0.0 else Normalize(m1.gadem.norm.error[1] - m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        elif error_do_or_not == "no":
            csv_line = {
                "Region_number": i,
                "src_Tmean(keV)": m1.gadem.Tmean.values[0],
                "src_Tsigma(keV)": m1.gadem.Tsigma.values[0],
                "src_Abundance": m1.gadem.abundanc.values[0],
                "src_Redshift": m1.gadem.Redshift.values[0],
                "src_norm(/400pi_arcmin2)": Normalize(m1.gadem.norm.values[0]),
                "NXB_gain_slope": s1.multiresponse[89].gain.slope.values[0],
                "Fe55_5.89keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian.norm.values[0]),
                "Fe55_6.49keV_norm(/400pi_arcmin2)": Normalize(m91.gaussian_2.norm.values[0]),
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "vgadem" and version == "v4":

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
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-0.5 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
        m1.vgadem.Al.frozen = True
        m1.vgadem.Si.frozen = False
        m1.vgadem.S.frozen = False
        m1.vgadem.Ar.frozen = False
        m1.vgadem.Ca.frozen = False
        m1.vgadem.Fe.frozen = False
        m1.vgadem.Ni.frozen = False

        Xset.restore("xtd_nxb.xcm")
        Xset.restore("xtd_cxb.xcm")
        Xset.restore("xtd_gh.xcm")
        Xset.restore("xtd_lhb.xcm")
        Xset.restore("xtd_0.8keV.xcm")
        Xset.restore("xtd_OI.xcm")

        m90 = AllModels(1, "det")
        m92 = AllModels(1, "cxb")
        m93 = AllModels(1, "gh")
        m94 = AllModels(1, "lhb")
        m95 = AllModels(1, "0.8keV")
        m96 = AllModels(1, "OI")

        def reverse_Normalize(norm_400pi):
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:19 src:21")

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
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "vgadem" and version == "v4Fe55":

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

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
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
            return norm_400pi * raytrace_source_area / (400*math.pi)

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
            Fit.error("1.0 src:2 src:3 src:19 src:21")

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


if arf_type == "beta" and model == "vgadem" and version == "v5":

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
            "C-stat",
            "d.o.f",
            "C-stat/d.o.f"
        ]

    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=parameter_names)
        writer.writeheader()

    for i in range(start, finish + 1):

        AllData.clear()
        AllModels.clear()

        AllData(f"1:1 xa{OBSID}xtd_{SUB}_beta_{i}.grp.gz")
        AllData.ignore("bad")
        s1 = AllData(1)
        s1.ignore("**-2.0 15.0-**")

        s1.multiresponse[89] = "xtend_nxb.rmf"
        s1.multiresponse[91] = f"xa{OBSID}xtd_{SUB}_src_{i}.rmf.gz"
        s1.multiresponse[91].arf = f"xa{OBSID}xtd_{SUB}_flat_{i}.arf.gz"

        s1.multiresponse[89].gain.slope = 1.0
        s1.multiresponse[89].gain.offset.frozen = True

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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
        m1.vgadem.Al.frozen = True
        m1.vgadem.Si.frozen = False
        m1.vgadem.S.frozen = False
        m1.vgadem.Ar.frozen = False
        m1.vgadem.Ca.frozen = False
        m1.vgadem.Fe.frozen = False
        m1.vgadem.Ni.frozen = False

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
            Fit.error("1.0 src:2 src:3 src:19 src:21")

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
                "C-stat": Fit.statistic,
                "d.o.f": Fit.dof,
                "C-stat/d.o.f": Fit.statistic/Fit.dof}

        with open(csv_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=parameter_names)
            writer.writerow(csv_line)


if arf_type == "beta" and model == "vgadem" and version == "v5Fe55":

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
        m1.vgadem.O.frozen = False
        m1.vgadem.Ne.frozen = False
        m1.vgadem.Na.frozen = True
        m1.vgadem.Mg.frozen = False
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

print(f"\n All processes have been completed. (ARF type = {arf_type}, model = {model}, version = {version})\n")