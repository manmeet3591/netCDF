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
    location = ('files/auswertung' + year + '.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4_CLASSIC')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset

def vertsum(data):
    data = numpy.sum(data, axis=1)
    return data

def yearsum(data):
    data = numpy.sum(data, axis=0)

for i in range (1850, 2009):
    # Create File and intialize dimensions
    year = str(i)
    dataset = writeFile(year)

    # Load Data
    print 'Loading data ...'
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    # Load variables into memory
    print 'Loading variables into memory'
    t_cya = nc_data.variables['t_cya'][:]
    t_lpp = nc_data.variables['t_lpp'][:]
    t_spp = nc_data.variables['t_spp'][:]
    dzt = nc_data.variables['dzt'][:]
    area_t = nc_data.variables['area_t'][:]
    # Calculation
    print 'Starting calculation'
    t_cya_int = t_cya*dzt*area_t
    t_cya_vert = vertsum(t_cya_int)
    t_lpp_int = t_lpp*dzt*area_t
    t_lpp_vert = vertsum(t_lpp_int)
    t_spp_int = t_spp*dzt*area_t
    t_spp_vert = vertsum(t_spp_int)

    # write data to file
    print 'writing data to file ...'
    #   create new variables
    print 'creating new variables ...'
    t_cya_ints = dataset.createVariable('t_cya_int', numpy.float64,('time', 'st_ocean', 'yt_ocean', 'xt_ocean'))
    t_cya_verts = dataset.createVariable('t_cya_vert', numpy.float64,('time', 'yt_ocean', 'xt_ocean'))
    t_lpp_ints = dataset.createVariable('t_lpp_int', numpy.float64,('time', 'st_ocean', 'yt_ocean', 'xt_ocean'))
    t_lpp_verts = dataset.createVariable('t_lpp_vert', numpy.float64,('time', 'yt_ocean', 'xt_ocean'))
    t_spp_ints = dataset.createVariable('t_spp_int', numpy.float64,('time', 'st_ocean', 'yt_ocean', 'xt_ocean'))
    t_spp_verts = dataset.createVariable('t_spp_vert', numpy.float64,('time', 'yt_ocean', 'xt_ocean'))


    #   write data
    t_cya_ints[:] = t_cya_int
    t_cya_verts[:] = t_cya_vert
    t_lpp_ints[:] = t_lpp_int
    t_lpp_verts[:] = t_lpp_vert
    t_spp_ints[:] = t_spp_int
    t_spp_verts[:] = t_spp_vert

    # plotting figures
    print 'Plotting figure ...'
    fig1 = pl.figure(1)
    fig1.suptitle(year, fontsize=16)
    pl.subplot(3,3,1)
    pl.pcolor(t_cya_vert[0,:,:])
    pl.title('Cya January')
    pl.subplot(3,3,2)
    pl.pcolor(t_cya_vert[4,:,:])
    pl.title('Cya May')
    pl.subplot(3,3,3)
    pl.pcolor(t_cya_vert[7,:,:])
    pl.title('Cya August')
    pl.colorbar()
    pl.subplot(3,3,4)
    pl.pcolor(t_lpp_vert[0,:,:])
    pl.title('Lpp January')
    pl.subplot(3,3,5)
    pl.pcolor(t_lpp_vert[4,:,:])
    pl.title('Lpp May')
    pl.subplot(3,3,6)
    pl.pcolor(t_lpp_vert[7,:,:])
    pl.title('Lpp August')
    pl.colorbar()
    pl.subplot(3,3,7)
    pl.pcolor(t_spp_vert[0,:,:])
    pl.title('Spp January')
    pl.subplot(3,3,8)
    pl.pcolor(t_spp_vert[4,:,:])
    pl.title('Spp May')
    pl.subplot(3,3,9)
    pl.pcolor(t_spp_vert[7,:,:])
    pl.title('Spp August')
    pl.colorbar()
    print('Saving figure ...')
    pl.savefig('figures/plot' + year + '.png')


    # close file
    dataset.close()
    print 'Done with year:', year
