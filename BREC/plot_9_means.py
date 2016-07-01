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

list_cya = []
for i in range (1851, 2010):
	year = str(i)
	nc_data = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung'+year+'.nc')
	cya_vert = (nc_data.variables['t_cya_vert'][:])
	cya_mean = numpy.mean(cya_vert, axis = 0)
	list_cya.append(cya_mean)

for i in range(len(list_cya)/9):
    year_begin = 1850+i*9
    year_end = 1850+i*9+9
    fig1 = pl.figure(1)
    fig1.suptitle('Cya-vert'+str(year_begin)+'-'+str(year_end), fontsize=16)
    fig1.figsize=(10,10)
    pl.subplot(3,3,1)
    pl.pcolor(list_cya[i*9+0][:,:])
    pl.title(str(year_begin+0),fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,2)
    pl.pcolor(list_cya[i*9+1][:,:])
    pl.title(str(year_begin+1), fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,3)
    pl.pcolor(list_cya[i*9+2][:,:])
    pl.title(str(year_begin+2), fontsize=8)
    pl.colorbar()
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,4)
    pl.pcolor(list_cya[i*9+3][:,:],)
    pl.title(str(year_begin+3), fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,5)
    pl.pcolor(list_cya[i*9+4][:,:])
    pl.title(str(year_begin+4), fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,6)
    pl.pcolor(list_cya[i*9+5][:,:])
    pl.title(str(year_begin+5), fontsize=8)
    pl.colorbar()
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,7)
    pl.pcolor(list_cya[i*9+6][:,:])
    pl.title(str(year_begin+6), fontsize=8)
    pl.subplot(3,3,8)
    pl.pcolor(list_cya[i*9+7][:,:])
    pl.title(str(year_begin+7), fontsize=8)
    pl.subplot(3,3,9)
    pl.pcolor(list_cya[i*9+8][:,:])
    pl.title(str(year_begin+8), fontsize=8)
    pl.colorbar()
    print 'saving figure to figures/plot', str(year_begin), '.png ...'
    pl.savefig('/silos/boergel/BREC/figures/plot_average' + str(year_begin) +'-'+str(year_end)+ '.png')
    pl.cla()
    pl.clf()

