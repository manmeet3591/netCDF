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



dataset = netCDF4.Dataset('calculation_test.nc','w',format='NETCDF4')
time = dataset.createDimension('time', None)
st_ocean = dataset.createDimension('st_ocean', 134)
yt_ocean = dataset.createDimension('yt_ocean', 242)
xt_ocean = dataset.createDimension('xt_ocean', 224)
# DATA
nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_1850/ocean_day3d.nc')
time = nc_data.variables['time']
st_ocean = nc_data.variables['st_ocean']
yt_ocean = nc_data.variables['yt_ocean']
xt_ocean = nc_data.variables['xt_ocean']
# other ...
dzt = nc_data.variables['dzt']
area_t = nc_data.variables['area_t']
temp = nc_data.variables['temp']
# Tracer variables
t_cya = nc_data.variables['t_cya']
t_lpp = nc_data.variables['t_lpp']
t_spp = nc_data.variables['t_spp']
t_n2 = nc_data.variables['t_n2']
t_o2 = nc_data.variables['t_o2']
t_n2 = nc_data.variables['t_n2']
t_nh4 = nc_data.variables['t_nh4']
t_no3 = nc_data.variables['t_no3']
t_po4 = nc_data.variables['t_po4']
t_zoo = nc_data.variables['t_zoo']

# CALC
print numpy.shape(dzt)
print numpy.shape(area_t)
print len(dzt[0,0,0,:])
for j in range(0, len(dzt[:,0,0,0])):
    print j
for i in range(0,len(dzt[0,0,0,:])):
    print i
#
dataset.close()
