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
    location = ('../files/ocean_trps/auswertung_ocean_trps_' + year + '.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset

def writeFilemean():
    location = ('../files/mean/ocean_trps_means.nc')
    dataset_mean = netCDF4.Dataset(location, 'w', format='NETCDF4')
    time = dataset_mean.createDimension('time', None)
    st_ocean = dataset_mean.createDimension('st_ocean', 134)
    yt_ocean = dataset_mean.createDimension('yt_ocean', 242)
    xt_ocean = dataset_mean.createDimension('xt_ocean', 224)
    return dataset_mean


variable_means = []
for i in range (1850,2009):
    # Intialize
    year = str(i)
    print year
    dataset = writeFile(year)
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' +year + '/ocean_trps.nc')
    nc_math = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')

    time = nc_data.variables['time']
    st_ocean = nc_data.variables['st_ocean']
    yt_ocean = nc_data.variables['yt_ocean']
    xt_ocean = nc_data.variables['xt_ocean']

    times = dataset.createVariable('time','f8',('time'))
    xt_oceans = dataset.createVariable('xt_ocean', 'f4',('xt_ocean'))
    yt_oceans = dataset.createVariable('yt_ocean', 'f4',('yt_ocean'))
    st_oceans = dataset.createVariable('st_ocean', 'i4', ('st_ocean'))
    attrs = [time,st_ocean,yt_ocean, xt_ocean]
    attrs_n = [times,st_oceans,yt_oceans, xt_oceans]
    i = 0
    for var in attrs:
        for j in var.ncattrs():
            attrs_n[i].setncattr(j,var.getncattr(j))
        i = i+1


    # Load variables
    t_no3_xflux_adv = nc_data.variables['t_no3_xflux_adv']
    t_no3_xflux_dif = nc_data.variables['t_no3_xflux_dif']
    t_nh4_xflux_adv = nc_data.variables['t_nh4_xflux_adv']
    t_nh4_xflux_dif = nc_data.variables['t_nh4_xflux_dif']
    t_det_xflux_adv = nc_data.variables['t_det_xflux_adv']
    t_det_xflux_dif = nc_data.variables['t_det_xflux_dif']

    area_t_data = nc_math.variables['area_t'][:]
    j = 0

    t_no3_xflux_adv_int = t_no3_xflux_adv[:]
    t_no3_xflux_dif_int = t_no3_xflux_dif[:]
    t_nh4_xflux_adv_int = t_nh4_xflux_adv[:]
    t_nh4_xflux_dif_int = t_nh4_xflux_dif[:]
    t_det_xflux_adv_int = t_det_xflux_adv[:]
    t_det_xflux_dif_int = t_det_xflux_dif[:]

    # Calculation
    # vertsum
    t_no3_xflux_adv_vert = numpy.sum(t_no3_xflux_adv_int, axis = 1)
    t_no3_xflux_dif_vert = numpy.sum(t_no3_xflux_dif_int, axis = 1)
    t_nh4_xflux_adv_vert = numpy.sum(t_nh4_xflux_adv_int, axis = 1)
    t_nh4_xflux_dif_vert = numpy.sum(t_nh4_xflux_dif_int, axis = 1)
    t_det_xflux_adv_vert = numpy.sum(t_det_xflux_adv_int, axis = 1)
    t_det_xflux_dif_vert = numpy.sum(t_det_xflux_dif_int, axis = 1)
    print 'Shape of t_no3_xflux_adv_vert: ', t_no3_xflux_adv_vert.shape
    # y-sum
    t_no3_xflux_adv_vy = numpy.sum(t_no3_xflux_adv_vert, axis = 1)
    t_no3_xflux_dif_vy = numpy.sum(t_no3_xflux_dif_vert, axis = 1)
    t_nh4_xflux_adv_vy = numpy.sum(t_nh4_xflux_adv_vert, axis = 1)
    t_nh4_xflux_dif_vy = numpy.sum(t_nh4_xflux_dif_vert, axis = 1)
    t_det_xflux_adv_vy = numpy.sum(t_det_xflux_adv_vert, axis = 1)
    t_det_xflux_dif_vy = numpy.sum(t_det_xflux_dif_vert, axis = 1)
    print 'Shape of t_no3_xflux_adv_vy: ', t_no3_xflux_adv_vy.shape
    # Boundary at x = 1
    t_no3_xflux_adv_vyx = t_no3_xflux_adv_vy[:,1]
    t_no3_xflux_dif_vyx = t_no3_xflux_dif_vy[:,1]
    t_nh4_xflux_adv_vyx = t_nh4_xflux_adv_vy[:,1]
    t_nh4_xflux_dif_vyx = t_nh4_xflux_dif_vy[:,1]
    t_det_xflux_adv_vyx = t_det_xflux_adv_vy[:,1]
    t_det_xflux_dif_vyx = t_det_xflux_dif_vy[:,1]
    print 'Shape of t_no3_xflux_adv_vyx: ', t_no3_xflux_adv_vyx.shape

    # write data to file
    t_no3_xflux_adv_vyxs = dataset.createVariable('t_no3_xflux_adv_vyx', numpy.float32,('time'))
    t_no3_xflux_dif_vyxs = dataset.createVariable('t_no3_xflux_dif_vyx', numpy.float32,('time'))
    t_nh4_xflux_adv_vyxs = dataset.createVariable('t_nh4_xflux_adv_vyx', numpy.float32,('time'))
    t_nh4_xflux_dif_vyxs = dataset.createVariable('t_nh4_xflux_dif_vyx', numpy.float32,('time'))
    t_det_xflux_adv_vyxs = dataset.createVariable('t_det_xflux_adv_vyx', numpy.float32,('time'))
    t_det_xflux_dif_vyxs = dataset.createVariable('t_det_xflux_dif_vyx', numpy.float32,('time'))
    # set units
    attrs_vars = [t_no3_xflux_adv, t_no3_xflux_dif, t_nh4_xflux_adv, t_nh4_xflux_dif, t_det_xflux_adv, t_det_xflux_dif]
    attrs_calc = [t_no3_xflux_adv_vyxs, t_no3_xflux_dif_vyxs, t_nh4_xflux_adv_vyxs, t_nh4_xflux_dif_vyxs, t_det_xflux_adv_vyxs, t_det_xflux_dif_vyxs]
    i = 0
    for var in attrs_calc:
       var.setncattr('units','mol *s^-1')
       var.setncattr('long_name',attrs_vars[i].long_name + 'vertical sum, sum over y and x = 2')
       var.setncattr('_FillValue', attrs_vars[i]._FillValue)
       var.setncattr('missing_value', attrs_vars[i].missing_value)
       var.setncattr('cell_methods', attrs_vars[i].cell_methods)
       var.setncattr('time_avg_info', attrs_vars[i].time_avg_info)
       var.setncattr('coordinates', attrs_vars[i].coordinates)
       i = i + 1

    # write data
    times[:] = time
    xt_oceans[:] = xt_ocean
    yt_oceans[:] = yt_ocean
    st_oceans[:] = st_ocean

    t_no3_xflux_adv_vyxs[:] = t_no3_xflux_adv_vyx
    t_no3_xflux_dif_vyxs[:] = t_no3_xflux_dif_vyx
    t_nh4_xflux_adv_vyxs[:] = t_nh4_xflux_adv_vyx
    t_nh4_xflux_dif_vyxs[:] = t_nh4_xflux_dif_vyx
    t_det_xflux_adv_vyxs[:] = t_det_xflux_adv_vyx
    t_det_xflux_dif_vyxs[:] = t_det_xflux_dif_vyx
    print 'Done with year ' + year
    dataset.close()
