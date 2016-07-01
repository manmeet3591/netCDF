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
nc_flux_sed = netCDF4.Dataset('/silos/boergel/BREC/files/mean/flux_sed_means.nc')

variable_list = []
variable_list_surf = []
variable_list_sed = []
variable_list_run = []
variable_list_surf_run = []
variable_list_sed_run = []
# Load flux3d
for i in nc_data.variables:
	variable_list.append(nc_data.variables[str(i)])
# load flux_surf
for i in nc_flux_surf.variables:
	variable_list_surf.append(nc_flux_surf.variables[str(i)])
# load flux_sed
for i in nc_flux_sed.variables:
	variable_list_sed.append(nc_flux_sed.variables[str(i)])

mean = np.ones((len(variable_list[0][:])))
mean_surf = np.ones((len(variable_list_surf[0][:])))
mean_sed = np.ones((len(variable_list_sed[0][:])))
ref_year_sed = np.ones((len(variable_list_sed[:])))
ref_year_surf = np.ones((len(variable_list_surf[:])))
ref_year_3d = np.ones((len(variable_list[:])))
diff_anomalie_sed = np.ones((len(variable_list_sed),len(variable_list_sed[0][:])))
diff_anomalie = np.ones((len(variable_list_surf),len(variable_list_surf[0][:])))
diff_anomalie_3d = np.ones((len(variable_list), len(variable_list[0][:])))

x3d = np.arange((len(variable_list[0][:])))
x3d = x3d+1850
xsurf = np.arange((len(variable_list_surf[0][:])))
xsurf = xsurf + 1850
xsed = np.arange((len(variable_list_sed[0][:])))
xsed = xsed + 1850

# Anomalie
#	SURF
for i in range(len(variable_list_surf)):
	ref_year_surf[i] = np.mean(variable_list_surf[i][126:155])
	for j in range(len(variable_list_surf[0][:])):
		diff_anomalie[i][j] = (variable_list_surf[i][j] - ref_year_surf[i])/ref_year_surf[i]	
# 	3D
for i in range(len(variable_list)):
	ref_year_3d[i] = np.mean(variable_list[i][126:155])
	for j in range(len(variable_list[0][:])):
		diff_anomalie_3d[i][j] = (variable_list[i][j] - ref_year_3d[i])/ref_year_3d[i]
#	SED
for i in range(len(variable_list_sed)):
        ref_year_sed[i] = np.mean(variable_list_sed[i][126:155])
        for j in range(len(variable_list_sed[0][:])):
                diff_anomalie_sed[i][j] = (variable_list_sed[i][j] - ref_year_sed[i])/ref_year_sed[i] 

# Running mean
#	3D
for i in range(len(variable_list)):
	for j in range(len(variable_list[i][:])):
		mean[j] = np.mean(variable_list[i][j:j+10])
	variable_list_run.append(mean)
#	SURF
for i in range(len(variable_list_surf)):
	for j in range(len(variable_list_surf[i][:])):
		mean_surf[j] = np.mean(variable_list_surf[i][j:j+10])
	variable_list_surf_run.append(mean_surf)
#	SED
for i in range(len(variable_list_sed)):
        for j in range(len(variable_list_sed[i][:])):
                mean_sed[j] = np.mean(variable_list_sed[i][j:j+10])
        variable_list_sed_run.append(mean_sed)



# plot
for i in range(len(variable_list)):
        fig = pl.figure(i)
	pl.cla()
	pl.clf()
        pl.plot(x3d,variable_list[i][:])
	pl.title(str(variable_list[i].name))
	pl.xlabel('Years')
	pl.savefig('../figures/mean/process/' + str(variable_list[i].name) + '_3d.png')

for i in range(len(variable_list)):
        fig = pl.figure(i)
        pl.cla()
        pl.clf()
        pl.plot(x3d,variable_list_run[i][:])
        pl.title(str(variable_list[i].name))
        pl.xlabel('Years')
        pl.savefig('../figures/mean/process/' + str(variable_list[i].name) + '_run.png')

for i in range(len(variable_list)):
        fig = pl.figure(i)
        pl.cla()
        pl.clf()
        pl.plot(x3d,diff_anomalie_3d[i][:], '--g')
        pl.axhline(y=0, xmin=0, xmax=1, hold=None, color='red')
        pl.title(str(variable_list[i].name))
        pl.xlabel('Years')
        pl.savefig('../figures/mean/process/' + str(variable_list[i].name) + '_anomalie.png')

for i in range(len(variable_list_surf)):
        fig = pl.figure(i)
        pl.cla()
        pl.clf()
        pl.plot(xsurf,variable_list_surf[i][:])
        pl.title(str(variable_list_surf[i].name))
        pl.xlabel('Years')
        pl.savefig('../figures/mean/process/' + str(variable_list_surf[i].name) + '_surf.png')

for i in range(len(variable_list_surf)):
	figs = pl.figure()
	pl.plot(xsurf,diff_anomalie[i][:], '--g')
	pl.title(str(variable_list_surf[i].name))
	pl.axhline(y=0, xmin=0, xmax=1, hold=None, color='red')
	pl.savefig('../figures/mean/process/' + str(variable_list_surf[i].name) + '_surf_anomalie.png')

for i in range(len(variable_list_surf)):
        figs = pl.figure()
        pl.plot(xsurf,variable_list_surf_run[i][:])
        pl.title(str(variable_list_surf[i].name))
	pl.xlabel('Years')
	pl.savefig('../figures/mean/process/' + str(variable_list_surf[i].name) + '_surf_run.png')


for i in range(len(variable_list_sed)):
        fig = pl.figure(i)
        pl.cla()
        pl.clf()
        pl.plot(xsed,variable_list_sed[i][:])
        pl.title(str(variable_list_sed[i].name))
        pl.xlabel('Years')
        pl.savefig('../figures/mean/process/' + str(variable_list_sed[i].name) + '_sed.png')

for i in range(len(variable_list_sed)):
        figs = pl.figure()
        pl.plot(xsed,diff_anomalie_sed[i][:], '--g')
        pl.axhline(y=0, xmin=0, xmax=1, hold=None, color='red')
        pl.savefig('../figures/mean/process/' + str(variable_list_sed[i].name) + '_sed_anomalie.png')

for i in range(len(variable_list_sed)):
        figs = pl.figure()
        pl.plot(xsed,variable_list_sed_run[i][:])
	pl.title(str(variable_list_sed[i].name))
        pl.savefig('../figures/mean/process/' + str(variable_list_sed[i].name) + '_sed_run.png')





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
