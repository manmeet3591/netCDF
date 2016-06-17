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
import calendar

def writeFile(year):
    location = ('files/flux_sed/auswertung_flux_' + year + '.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset

def writeFilemean():
    location = ('files/mean/flux_sed_means.nc')
    dataset_mean = netCDF4.Dataset(location, 'w', format='NETCDF4')
    time = dataset_mean.createDimension('time', None)
    st_ocean = dataset_mean.createDimension('st_ocean', 134)
    yt_ocean = dataset_mean.createDimension('yt_ocean', 242)
    xt_ocean = dataset_mean.createDimension('xt_ocean', 224)
    return dataset_mean


variable_means = []
for i in range (1850,2010):
    year = str(i)
    dataset = writeFile(year)
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ergom_flux_sed.nc')
    nc_math = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    # variables
    #variable_names = {}
    #variable_names[str(i)+"s"] = nc_data.variables[str(i)]
    #for i in variable_names:
    times = dataset.createVariable('time','f8',('time'))
    xt_oceans = dataset.createVariable('xt_ocean', 'f4',('xt_ocean'))
    yt_oceans = dataset.createVariable('yt_ocean', 'f4',('yt_ocean'))
    st_oceans = dataset.createVariable('st_ocean', 'i4', ('st_ocean'))
    
    area_t_data = nc_math.variables['area_t'][:]
    j = 0
    variable_list = []

    for i in nc_data.variables:
        variable_list.append(nc_data.variables[str(i)])

    shaping = variable_list[4].shape
    shaping = numpy.insert(shaping,0,13)
    variable_int = numpy.ones((shaping))
    variable_int = numpy.ma.array(variable_int)

    daysinyear = 365
    if calendar.isleap(int(year)) == True:
	daysinyear = daysinyear + 1
    for i in range(4,17):
        variable_int[j] = (variable_list[i][:]*area_t_data)
	print variable_list[i].name
        j = j + 1
 
    # fieldsum
    variable_field = numpy.sum(variable_int, axis = 2)
    variable_field = numpy.sum(variable_field, axis = 2)

    variable_mean = numpy.mean(variable_field, axis = 1)
    print variable_mean.shape
    variable_means.append(variable_mean)
    #write Data to file
    
    attrs = [variable_list[2],variable_list[1], variable_list[0]]
    attrs_n = [times,yt_oceans, xt_oceans]
    i = 0
    for var in attrs:
        for j in var.ncattrs():
            attrs_n[i].setncattr(j,var.getncattr(j))
        i = i+1

    times[:] = variable_list[2]
    xt_oceans[:] = variable_list[0]
    yt_oceans[:] = variable_list[1]

    count = 0
    for i in variable_list[4:17]:
	k = dataset.createVariable(str(i.name)+'_field', numpy.float32,('time'))
	k.setncattr('units','mol * d-1')
        k.setncattr('long_name',i.long_name + '+ horizontal')
        k.setncattr('_FillValue',i._FillValue)
        k.setncattr('missing_value', i.missing_value)
        k.setncattr('cell_methods', i.cell_methods)
        k.setncattr('time_avg_info', i.time_avg_info)
        k.setncattr('coordinates', i.coordinates)
	k[:] = variable_field[count]
	count = count + 1

    print 'Done with year' + year
    dataset.close()

a = numpy.ones((len(variable_means[:])))
dataset_mean = writeFilemean()
for j in range (len(variable_means[0][:])):
    for i in range (len(variable_means[:])):
	a[i] = variable_means[i][j]
    means = dataset_mean.createVariable(str(variable_list[4+j].name)+'_mea',numpy.float32, 'time')
    means[:] = a       
