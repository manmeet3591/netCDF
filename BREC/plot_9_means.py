import sys, string
from matplotlib import rc
import numpy
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


def do_plot(n, X, title):
    #plt.clf()
    pl.subplot(3, 3, n)
    pl.pcolor(X, vmin=-min_value, vmax=max_value)
    pl.title(title, fontsize=8)

nc_list = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung1850.nc')
variable_list = []
for i in nc_list.variables:
        variable_list.append(nc_list.variables[str(i)])
placeholder = variable_list[4:16]

for j in placeholder:
    count = j.name
    list_variables = []
    max_value = 0.0
    min_value = 0.0
    values = 0
    print 'Plotting', str(count)
    for i in range (1851, 2010):
        year = str(i)
        variables_list = []
        nc_data = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung'+year+'.nc')
        variable = nc_data.variables[str(count)]
        mean = numpy.mean(variable[:], axis = 0)
        list_variables.append(mean)
        if list_variables[values].max() > max_value:
            max_value = list_variables[values].max()
        if list_variables[values].min() < min_value:
            min_value = list_variables[values].min()
        values = values + 1


    for i in range(len(list_variables)/9):
        year_begin = 1850+i*9
        year_end = 1850+i*9+9
        fig1 = pl.figure(1)
        fig1.suptitle(str(variable.name) +str(year_begin)+'-'+str(year_end), fontsize=16)
        fig1.figsize=(10,10)
        do_plot(1, list_variables[i*9+0][:,:],str(year_begin+0))
        pl.subplots_adjust(hspace = .5)
        pl.subplot(3,3,2)
        pl.pcolor(list_variables[i*9+1][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+1), fontsize=8)
        pl.subplots_adjust(hspace = .5)
        pl.subplot(3,3,3)
        pl.pcolor(list_variables[i*9+2][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+2), fontsize=8)
        pl.colorbar()
        pl.subplots_adjust(hspace = .5)
        pl.subplot(3,3,4)
        pl.pcolor(list_variables[i*9+3][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+3), fontsize=8)
        pl.subplots_adjust(hspace = .5)
        pl.subplot(3,3,5)
        pl.pcolor(list_variables[i*9+4][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+4), fontsize=8)
        pl.subplots_adjust(hspace = .5)
        pl.subplot(3,3,6)
        pl.pcolor(list_variables[i*9+5][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+5), fontsize=8)
        pl.colorbar()
        pl.subplots_adjust(hspace = .5)
        pl.subplot(3,3,7)
        pl.pcolor(list_variables[i*9+6][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+6), fontsize=8)
        pl.subplot(3,3,8)
        pl.pcolor(list_variables[i*9+7][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+7), fontsize=8)
        pl.subplot(3,3,9)
        pl.pcolor(list_variables[i*9+8][:,:], vmin=min_value, vmax=max_value)
        pl.title(str(year_begin+8), fontsize=8)
        pl.colorbar()
        print 'saving figure to figures/plot', str(year_begin), '.png ...'
        pl.savefig('/silos/boergel/BREC/figures/plot_average'+ str(variable.name) + str(year_begin) +'-'+str(year_end)+ '.png')
        pl.cla()
        pl.clf()

