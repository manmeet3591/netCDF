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


nc_data = netCDF4.Dataset('/silos/boergel/BREC/files/mean/flux3d_means.nc')
nc_flux_surf = netCDF4.Dataset('/silos/boergel/BREC/files/mean/flux_surf_means.nc')
variable_list = []
variable_list_surf = []
variable_list_run = []
# Load flux3d
for i in nc_data.variables:
	variable_list.append(nc_data.variables[str(i)])
# load flux_surf
for i in nc_flux_surf.variables:
	variable_list_surf.append(nc_flux_surf.variables[str(i)])

mean = np.ones((len(variable_list[0][:])))
ref_year_surf = np.ones((len(variable_list_surf[:])))
diff_anomalie = np.ones((len(variable_list_surf),len(variable_list_surf[0][:])))


x = np.arange((len(variable_list_surf[0][:])))
x = x+1850
# Anomalie
for i in range(len(variable_list_surf)):
	ref_year_surf[i] = np.mean(variable_list_surf[i][126:155])
	for j in range(len(variable_list_surf[0][:])):
		diff_anomalie[i][j] = (variable_list_surf[i][j] - ref_year_surf[i])/ref_year_surf[i]	

# Running mean
for i in range(len(variable_list)):
	for j in range(len(variable_list[i][:])):
		mean[j] = np.mean(variable_list[i][j:j+10])
	variable_list_run.append(mean)
# plot
for i in range(len(variable_list)):
        fig = pl.figure(i)
	pl.cla()
	pl.clf()
        pl.plot(x,variable_list[i][:])
	pl.title(str(variable_list[i].name))
	pl.xlabel('Years')
	pl.savefig('../figures/mean/process/' + str(variable_list[i].name) + '.png')


for i in range(len(variable_list)):
        fig = pl.figure(i)
        pl.cla()
        pl.clf()
        pl.plot(x,variable_list_run[i][:])
        pl.title(str(variable_list[i].name))
        pl.xlabel('Years')
        pl.savefig('../figures/mean/process/' + str(variable_list[i].name) + '_run.png')

for i in range(len(variable_list_surf)):
	figs = pl.figure()
	pl.plot(x,diff_anomalie[i][:], '--g')
	pl.axhline(y=0, xmin=0, xmax=1, hold=None, color='red')
	pl.savefig('../figures/mean/process/' + str(variable_list[i].name) + '_anomalie.png')

pos_values = diff_anomalie[5].copy()
neg_values = diff_anomalie[5].copy()

pos_values[pos_values < 0.0000000] = np.nan
neg_values[neg_values > 0.00000000] = np.nan

figs = pl.figure()
pl.plot(diff_anomalie[5][:], '--g')
pl.plot(pos_values, color='r')
pl.plot(neg_values, color='b')
pl.axhline(y=0, xmin=0, xmax=1, hold=None, color='red')
figs.show()
