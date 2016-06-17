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
    location = ('/silos/boergel/BREC/files/ocean_trps/auswertung_add_ocean_trps' + year + '.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset



rho = 1000
variable_means = []
for i in range (1850,2010):
    year = str(i)
    dataset = writeFile(year)
    nc_ocean = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    nc_ocean_trps = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_trps.nc')

    time = nc_ocean.variables['time']
    st_ocean = nc_ocean.variables['st_ocean']
    yt_ocean = nc_ocean.variables['yt_ocean']
    xt_ocean = nc_ocean.variables['xt_ocean']

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

    tx_trans  = nc_ocean_trps.variables['tx_trans'][:]
    t_cya = nc_ocean.variables['t_cya']
    t_lpp = nc_ocean.variables['t_lpp']
    t_spp = nc_ocean.variables['t_spp']
    t_zoo = nc_ocean.variables['t_zoo']
    t_det = nc_ocean.variables['t_det']

    flux_cya = tx_trans*t_cya[:]
    flux_lpp = tx_trans*t_lpp[:]
    flux_spp = tx_trans*t_spp[:]
    flux_zoo = tx_trans*t_zoo[:]
    flux_det = tx_trans*t_det[:]


    # vertsum
    flux_cya_vert = numpy.sum(flux_cya, axis = 1)
    flux_lpp_vert = numpy.sum(flux_lpp, axis = 1)
    flux_spp_vert = numpy.sum(flux_spp, axis = 1)
    flux_zoo_vert = numpy.sum(flux_zoo, axis = 1)
    flux_det_vert = numpy.sum(flux_det, axis = 1)

    # y-sum
    flux_cya_vy = numpy.sum(flux_cya_vert, axis = 1)
    flux_lpp_vy = numpy.sum(flux_lpp_vert, axis = 1)
    flux_spp_vy = numpy.sum(flux_spp_vert, axis = 1)
    flux_zoo_vy = numpy.sum(flux_zoo_vert, axis = 1)
    flux_det_vy = numpy.sum(flux_det_vert, axis = 1)

    # Boundary at x = 1
    flux_cya_vyx = flux_cya_vy[:,1]
    flux_lpp_vyx = flux_lpp_vy[:,1]
    flux_spp_vyx = flux_spp_vy[:,1]
    flux_zoo_vyx = flux_zoo_vy[:,1]
    flux_det_vyx = flux_det_vy[:,1]


    flux_cya_vyxs = dataset.createVariable('flux_cya_vyx', numpy.float32,('time'))
    flux_lpp_vyxs = dataset.createVariable('flux_lpp_vyx', numpy.float32,('time'))
    flux_spp_vyxs = dataset.createVariable('flux_spp_vyx', numpy.float32,('time'))
    flux_zoo_vyxs = dataset.createVariable('flux_zoo_vyx', numpy.float32,('time'))            
    flux_det_vyxs = dataset.createVariable('flux_det_vyx', numpy.float32,('time'))


    attrs_vars = [flux_cya_vyxs,flux_lpp_vyxs,flux_spp_vyxs,flux_zoo_vyxs,flux_det_vyxs]
    i = 0
    for var in attrs_vars:
        var.setncattr('units','mol * s^-1')
        var.setncattr('long_name', 'flux at ...')
        i = i +1
    #write Data to file

    #write Data
    times[:] = time
    xt_oceans[:] = xt_ocean
    yt_oceans[:] = yt_ocean
    st_oceans[:] = st_ocean

    flux_cya_vyxs[:] = flux_cya_vyx
    flux_lpp_vyxs[:] = flux_lpp_vyx
    flux_spp_vyxs[:] = flux_spp_vyx
    flux_zoo_vyxs[:] = flux_zoo_vyx
    flux_det_vyxs[:] = flux_det_vyx


    print 'Done with year' + year
    dataset.close()
