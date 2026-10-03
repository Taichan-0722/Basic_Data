import os
import sys
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import collections
import lmfit as lf
import lmfit.models as lfm
from scipy.stats import poisson

def poisson_distribution(x, mu):
    return poisson.pmf(x, mu)

def FileExistence(dir_path, keyword):
    files = os.listdir(dir_path)
    check = False

    for file in files:
        if keyword in file:
            check = True
    
    return check

# get parameters
print("Input working directory ->", end=" ")
DIR = input()
if os.path.isdir(DIR) == False:
    print("Error: Such directry does not exist.")
    sys.exit(1)

print("Input image file ->", end=" ")
infile = input()
temp = os.path.join(DIR, infile)
if os.path.isfile(temp) == False:
    print("Error:", infile, "does not exist in", DIR)
    sys.exit(1)

print("Input exposure map image file ->", end=" ")
expmapfile = input()
temp = os.path.join(DIR, expmapfile)
if os.path.isfile(temp) == False:
    print("Error:", expmapfile, "does not exist in", DIR)
    sys.exit(1)

print("Output image file name ->", end=" ")
outimgfile = input()
temp = os.path.join(DIR, outimgfile)
if os.path.exists(temp) == True:
    print("Error:", outimgfile, "already exists in", DIR)
    sys.exit(1)

print("Output region file name ->", end=" ")
outregfile = input()
temp = os.path.join(DIR, outregfile)
if os.path.exists(temp) == True:
    print("Error:", outregfile, "already exists in", DIR)
    sys.exit(1)

print("Threshold for simple counts cut ( > 0 ) ->", end=" ")
maxcount = int(input())
if maxcount <= 0:
    print("Error: this value must be > 0.", maxcount, "<= 0")
    sys.exit(1)

print("Poisson clean cell size (odd integer > 1) ->", end=" ")
N = int(input())
if N % 2 == 0:
    print("Error: cell size must be odd integer.", N, "is even.")
    sys.exit(1)
elif N <= 1:
    print("Error: cell size must be > 1.", N, "<= 1")
    sys.exit(1)

print("Poisson probability threshold (multiple of sigma > 0) ->", end=" ")
sigma = float(input())
if sigma <= 0:
    print("Error: this value must be > 0.", sigma, "<= 0")
    sys.exit(1)

print("")

bigval_x = np.array([])
bigval_y = np.array([])
under = int((N - 1) / 2.0)
upper = int((N + 1) / 2.0)

os.chdir(DIR)

# get event image
hdu = fits.open(infile)
data = hdu[0].data
hdu.close()

# get exposuremap image 
hdu = fits.open(expmapfile)
expdata = hdu[0].data
hdu.close()

# Apply a simple counts cut threshold
print("Searching pixels with counts over", maxcount, "...")
for i in range(0, 1810):
    for j in range(0, 1810):
        if data[j][i] > maxcount:
            bigval_x = np.append(bigval_x, i)
            bigval_y = np.append(bigval_y, j)
            data[j][i] = 0
print("Finished.\n")
simple_x = bigval_x
simple_y = bigval_y

# reject with poisson fitting
print("Searching flickering pixels with counts under", maxcount, "...")
for i in range(0, 1810):
    if i % 50 == 0:
        print("Now", 1+i*1810, "/ 3276100")

    for j in range(0, 1810):
        selected_data = []

        if data[j][i] > 1:
            for x in range(i-under,i+upper):
                for y in range(j-under,j+upper):
                    if expdata[y][x] != 0:
                        if x in simple_x and y in simple_y:
                            continue
                        else:
                            selected_data.append(data[y][x])

            if len(selected_data) != 0:
                average = sum(selected_data) / len(selected_data)
                temp = collections.Counter(selected_data)
                hist = np.asarray(sorted(list(temp.items())))

                data_x = np.array([])
                data_y = np.array([])
                for k in range(len(hist)):
                    data_x = np.append(data_x, hist[k][0])
                    data_y = np.append(data_y, hist[k][1]/len(selected_data))

                model = lf.Model(poisson_distribution)
                params = model.make_params()
                params.add("mu", value=average, min=0.01, max=average*2.0)

                result = model.fit(x=data_x, data=data_y, params=params, method='leastsq')
                mu = result.params['mu'].value

                threshold = mu + sigma * np.sqrt(mu)

                if data[j][i] > threshold:
                    bigval_x = np.append(bigval_x, i)
                    bigval_y = np.append(bigval_y, j)
print("Finished.\n")

# output file
print("Extracting region and image file ...")
with open(outregfile, "w") as f:
    for i in range(len(bigval_x)):
        x = int(bigval_x[i])
        y = int(bigval_y[i])
        data[y][x] = 0
        text = "-box(" + str(x+1) + "," + str(y+1) + ",1,1,0)\n"
        f.write(text)

hdu = fits.PrimaryHDU(data)
img = fits.HDUList([hdu])
img.writeto(outimgfile)
img.close()
print("Completed.")
print(len(bigval_x), "pixels are removed as flickering pixels.")
