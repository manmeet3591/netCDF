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

def writeFile(year):
    location = ('/silos/boergel/BREC/files/flux3d/auswertung_flux_' + year + '.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset

def writeFilemean():
    location = ('/silos/boergel/BREC/files/mean/flux3d_means.nc')
    dataset_mean = netCDF4.Dataset(location, 'w', format='NETCDF4')
    time = dataset_mean.createDimension('time', None)
    st_ocean = dataset_mean.createDimension('st_ocean', 134)
    yt_ocean = dataset_mean.createDimension('yt_ocean', 242)
    xt_ocean = dataset_mean.createDimension('xt_ocean', 224)
    return dataset_mean


variable_means = []
for i in range (1850,2009):
    # initialize
    year = str(i)
    print year
    dataset = writeFile(year)
    times = dataset.createVariable('time','f8',('time'))
    xt_oceans = dataset.createVariable('xt_ocean', 'f4',('xt_ocean'))
    yt_oceans = dataset.createVariable('yt_ocean', 'f4',('yt_ocean'))
    st_oceans = dataset.createVariable('st_ocean', 'i4', ('st_ocean'))
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ergom_flux3d.nc')
    nc_math = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    dzt_data = nc_math.variables['dzt'][:]
    area_t_data = nc_math.variables['area_t'][:]

    # Load all variables from nc_data
    j = 0
    variable_list = []

    for i in nc_data.variables:
        variable_list.append(nc_data.variables[str(i)])
   
    # create array for all variables from nc_data
    shaping = variable_list[6].shape
    shaping = numpy.insert(shaping,0,30)
    variable_int = numpy.ones((shaping))
    variable_int = numpy.ma.array(variable_int)

    for i in range(6,36):
        variable_int[j] = (variable_list[i][:]*dzt_data*area_t_data)
        j = j + 1

    # vertsum
    variable_vert = numpy.sum(variable_int, axis = 2)
    # fieldsum
    variable_field = numpy.sum(variable_vert, axis = 2)
    variable_field = numpy.sum(variable_field, axis = 2)
    # mean
    variable_mean = numpy.mean(variable_field, axis = 1)
    variable_means.append(variable_mean)
    print variable_field.shape 
    
    #write data to file
    attrs = [variable_list[4],variable_list[2],variable_list[1], variable_list[0]]
    attrs_n = [times,st_oceans,yt_oceans, xt_oceans]
    i = 0
    for var in attrs:
        for j in var.ncattrs():
            attrs_n[i].setncattr(j,var.getncattr(j))
        i = i+1

    times[:] = variable_list[4]
    xt_oceans[:] = variable_list[0]
    yt_oceans[:] = variable_list[1]
    st_oceans[:] = variable_list[2] 

    # set attributes for .nc file
    count = 0
    for i in variable_list[6:36]:
        j = dataset.createVariable(str(i.name)+'_vert',  numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
        j.setncattr('units','mol * m^3/kg')
        j.setncattr('long_name',i.long_name + '+ vertical sum')
        j.setncattr('_FillValue',i._FillValue)
        j.setncattr('missing_value', i.missing_value)
        j.setncattr('cell_methods', i.cell_methods)
        j.setncattr('time_avg_info', i.time_avg_info)
        j.setncattr('coordinates', i.coordinates)
        j[:] = variable_vert[count]
	k = dataset.createVariable(str(i.name)+'_field', numpy.float32,('time'))
	k.setncattr('units','mol * m^3/kg')
        k.setncattr('long_name',i.long_name + '+ horizontal + vertical sum')
        k.setncattr('_FillValue',i._FillValue)
        k.setncattr('missing_value', i.missing_value)
        k.setncattr('cell_methods', i.cell_methods)
        k.setncattr('time_avg_info', i.time_avg_info)
        k.setncattr('coordinates', i.coordinates)
	k[:] = variable_field[count]
	count = count + 1
         
   
    print 'Done with year' + year
    dataset.close()

# Write mean
a = numpy.ones((len(variable_means[:])))
dataset_mean = writeFilemean()
for j in range (len(variable_means[0][:])):
    for i in range (len(variable_means[:])):
        a[i] = variable_means[i][j]
    means = dataset_mean.createVariable(str(variable_list[6+j].name)+'_mea',numpy.float32, 'time')
    means[:] = a

