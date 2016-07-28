# Relation between salt concetration and freshwater inlet #
# -*- coding: utf-8 -*-
from scipy.interpolate import interp1d
import sys, string
from matplotlib import rc
import numpy as np
import pylab as pl
import netCDF4
import time as t
import datetime
from dateutil.parser import parse
from pylab import load, meshgrid, title, arange, show
from netcdftime import utime
import scipy.io
import matplotlib as mpl
import argparse
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import datetime as dt
from scipy import stats

salt_rate_list = []
frischwasser_list = []
for i in range(1850,2009):
    year = str(i)
    year_before = str(i+1)
    year_3_before = str(i)
    nc_ice = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ice_day.nc')
    nc_data_year = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung' + year + '.nc')
    nc_data_year_before = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung' + year_before + '.nc')
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')

    # rates
    # shape [t][y][x] in [kg m^-2 s^-1]
    snow = nc_ice.variables['SNOWFL'][:]
    rain = nc_ice.variables['RAIN'][:]
    runoff = nc_ice.variables['RUNOFF'][:]
    evap = nc_ice.variables['EVAP'][:]


    area_t = nc_data.variables['area_t'][:]

    dayspermonth = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


    snow_int = snow * area_t
    rain_int = rain * area_t
    runoff_int = runoff * area_t
    evap_int = evap * area_t

    # preprocessed salt
    salt_year_before = nc_data_year_before.variables['salt_field'][:]
    salt_year = nc_data_year.variables['salt_field'][:]

    #print 'salt', salt_year

    # horizontal sum freshwater
    snow_x = np.sum(snow, axis = 1)
    rain_x = np.sum(rain, axis = 1)
    runoff_x = np.sum(runoff, axis = 1)
    evap_x = np.sum(evap, axis = 1)

    snow_yx = np.sum(snow_x, axis = 1)
    rain_yx = np.sum(rain_x, axis = 1)
    runoff_yx = np.sum(runoff_x, axis = 1)
    evap_yx = np.sum(evap_x, axis = 1)

    # multiply with day/month * 24h/day * 60 min/h * 60 s/min to obtain kg/month
    for i in range(len(snow_yx)):
        snow_yxt = snow_yx * dayspermonth[i] * 24 * 60 * 60
        rain_yxt = rain_yx * dayspermonth[i] * 24 * 60 * 60
        runoff_yxt = runoff_yx * dayspermonth[i] * 24 * 60 * 60
        evap_yxt = evap_yx * dayspermonth[i] * 24 * 60 * 60

    # average over year
    sum_frischwasser = snow_yxt + rain_yxt + runoff_yxt - evap_yxt
    #print 'sum_frischwasser', sum_frischwasser.shape
    sum_frischwasser_mean = np.mean(sum_frischwasser)

    # salt changed ... mean to avoid error in broadcast for 11 months
    salt_year_before_mean = np.mean(salt_year_before)
    salt_year_mean = np.mean(salt_year)

    salt_rate_mean = salt_year_mean - salt_year_before_mean

    # append to list in order to evaluate later
    salt_rate_list.append(salt_year)
    frischwasser_list.append(sum_frischwasser)

    #print 'Done with year ' + year

# normalize Dataset
#salt_normed = salt_rate_list / np.sum(salt_rate_list)
#freshwater_normed = frischwasser_list / np.sum(frischwasser_list)
freshwater_full = []
salt_full = []
for i in range(len(frischwasser_list)):
    freshwater_full = np.append(freshwater_full, frischwasser_list[i][:])
    salt_full = np.append(salt_full,salt_rate_list[i][:])

import scipy.signal as signal


# First, design the Buterworth filter
N  = 2    # Filter order
Wn = 0.01 # Cutoff frequency
B, A = signal.butter(N, Wn, output='ba')

# Second, apply the filter
saltf = signal.filtfilt(B,A, salt_full)
freshwaterf = signal.filtfilt(B,A, freshwater_full)

check = 0
for i in range(len(freshwaterf)):
    corrcoff = np.corrcoef(freshwaterf, np.roll(saltf,i-1))
    print corrcoff[0,1]
    if abs(corrcoff[0,1]) > abs(check):
        print abs(corrcoff[0,1]), "is bigger then " , check
        check = corrcoff[0,1]
        label = i

print check, i

fig = pl.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot(saltf, '--r', label = 'Salt concentration')
ax2 = ax.twinx()
lns2 = ax2.plot(freshwaterf, '-b', label = 'Freshwater flow')

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.grid()
ax.set_xlabel("Time (months)")
ax.set_ylabel(r"Salt concentration")
ax2.set_ylabel(r"Freshwater flow")

fig = pl.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot(salt_full, '--r', label = 'Salt concentration')
ax2 = ax.twinx()
lns2 = ax2.plot(freshwater_full, '-b', label = 'Freshwater flow')

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.grid()
ax.set_xlabel("Time (months)")
ax.set_ylabel(r"Salt concentration")
ax2.set_ylabel(r"Freshwater flow")

pl.show()



