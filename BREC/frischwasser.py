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
    temp_year = nc_data_year.varialbes['temp_field'][:]

    print 'salt', salt_year

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
    print 'sum_frischwasser', sum_frischwasser.shape
    sum_frischwasser_mean = np.mean(sum_frischwasser)

    # salt changed ... mean to avoid error in broadcast for 11 months
    salt_year_before_mean = np.mean(salt_year_before)
    salt_year_mean = np.mean(salt_year)

    salt_rate_mean = salt_year_mean - salt_year_before_mean

    # append to list in order to evaluate later
    salt_rate_list.append(salt_year_mean)
    frischwasser_list.append(sum_frischwasser_mean)

    print 'Done with year ' + year

# print results
slope, intercept, r_value, p_value, std_err = stats.linregress(salt_rate_list, frischwasser_list)
pl.scatter(salt_rate_list, frischwasser_list)
#print slope, intercept, r_value, p_value, std_err
#print 'R^2', r_value**2

# running mean
salt_run_mean = np.ones((len(salt_rate_list)))
freshwater_run_mean = np.ones((len(frischwasser_list)))
for j in range(len(salt_rate_list)):
    salt_run_mean[j] = np.mean(salt_rate_list[j:j+4])
    freshwater_run_mean[j] = np.mean(frischwasser_list[j:j+4])

# normalize Dataset
salt_normed = salt_rate_list / np.sum(salt_rate_list)
freshwater_normed = frischwasser_list / np.sum(frischwasser_list)

xcross = np.correlate(freshwater_normed,salt_normed, "full")

salt_normed_cross = np.ones((len(salt_normed)-1))
freshwater_normed_cross = np.ones((len(salt_normed)-1))

shift = len(salt_normed)-xcross.argmax()

try:
    for i in range(len(salt_normed_cross)):
        print i
        salt_normed_cross[i] = salt_normed[i+shift]
        freshwater_normed_cross[i] = freshwater_normed[i]
except:
    print 'run ipython and check shift size'

crosscoeff = np.corrcoef(salt_normed_cross, freshwater_normed_cross)

print xcross.argmax()

# Fit Functions
x = np.arange((len(frischwasser_list)))
x = x+1850
#temp_mean = np.mean(temp_field)
#temp_std = np.std(temp_field)
fit_salt = np.polyfit(x,salt_rate_list,6)
fit_fn_salt = np.poly1d(fit_salt)
fit_frisch = np.polyfit(x,frischwasser_list,6)
fit_fn_frisch = np.poly1d(fit_frisch)

#pl.plot(x,fit_fn_salt(x))
#pl.plot(x,fit_fn_frisch(x))


a = pl.figure()
ax1 = pl.subplot(211)
pl.plot(x,salt_rate_list)
pl.plot(x,salt_run_mean, '--r')
pl.setp(ax1.get_xticklabels(), fontsize=6)
pl.setp(ax1.get_yticklabels(), fontsize=6)

ax2 = pl.subplot(212, sharex=ax1)
pl.plot(x,frischwasser_list)
pl.plot(x, freshwater_run_mean, '--r')
pl.setp(ax2.get_xticklabels(), fontsize=6)
pl.setp(ax2.get_yticklabels(), fontsize=6)
a.subplots_adjust(hspace=0.2)

pl.show()
