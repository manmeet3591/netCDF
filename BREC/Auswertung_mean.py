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

def writeFile():
    location = ('files/mean/mean_ocean.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset



variable_means = []
for i in range (1850,2010):
    year = str(i)
    nc_data = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung'+year+'.nc')
    variable_list = []
    for i in nc_data.variables:
	variable_list.append(nc_data.variables[str(i)])
    variable_mean = numpy.ones((14))
    variable_mean = numpy.ma.array(variable_mean)
    j = 0
    for i in range(16,30):
	variable_mean[j] = numpy.mean(variable_list[i])
	j = j+1
	print i,j
    variable_means.append(variable_mean)
    print 'done with year' + year
   
a = numpy.ones((len(variable_means[:])))
dataset_mean = writeFile()
for j in range (len(variable_means[0][:])):
    for i in range (len(variable_means[:])):
        a[i] = variable_means[i][j]
    means = dataset_mean.createVariable(str(variable_list[16+j].name)+'_mea',numpy.float32, 'time')
    means[:] = a
