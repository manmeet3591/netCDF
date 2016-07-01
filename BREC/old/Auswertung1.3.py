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
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset

def vertsum(data):
    data = numpy.sum(data, axis=1)
    return data

for i in range (1850, 1851):
    # Create File and intialize dimensions    
    year = str(i)
    dataset = writeFile(year)

    # Load Data
    print 'loading data ...'
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    # Load variables into memory
    print 'loading variables into memory ...'
    # Coordinates
    time = nc_data.variables['time']
    st_ocean = nc_data.variables['st_ocean']
    yt_ocean = nc_data.variables['yt_ocean']
    xt_ocean = nc_data.variables['xt_ocean']
    # other ...
    dzt = nc_data.variables['dzt'][:]
    area_t = nc_data.variables['area_t'][:]
    temp = nc_data.variables['temp'][:]
    # Tracer variables
    t_cya = nc_data.variables['t_cya'][:]
    t_lpp = nc_data.variables['t_lpp'][:]
    t_spp = nc_data.variables['t_spp'][:]
    t_n2 = nc_data.variables['t_n2'][:]
    t_o2 = nc_data.variables['t_o2'][:]
    t_n2 = nc_data.variables['t_n2'][:]
    t_nh4 = nc_data.variables['t_nh4'][:]
    t_no3 = nc_data.variables['t_no3'][:]
    t_po4 = nc_data.variables['t_po4'][:]
    t_zoo = nc_data.variables['t_zoo'][:]
    

    # Calculation
    print 'starting calculation ...'
    t_cya_int = t_cya*dzt*area_t
    t_cya_vert = vertsum(t_cya_int)
    t_lpp_int = t_lpp*dzt*area_t
    t_lpp_vert = vertsum(t_lpp_int)
    t_spp_int = t_spp*dzt*area_t
    t_spp_vert = vertsum(t_spp_int)
    t_cya_field = numpy.sum(t_cya_vert, axis=1)
    t_cya_field = numpy.sum(t_cya_field, axis=1)
    t_spp_field = numpy.sum(t_spp_vert, axis=1)
    t_spp_field = numpy.sum(t_spp_field, axis=1)
    t_lpp_field = numpy.sum(t_lpp_vert, axis=1)
    t_lpp_field = numpy.sum(t_lpp_field, axis=1)
    print numpy.shape(t_cya_field)

    

    # write data to file
    print 'writing data to file ...'
    #   create new variables
    print 'creating new variables ...'
    
    times = dataset.createVariable('time','f8',('time'))
    xt_oceans = dataset.createVariable('xt_ocean', 'f4',('xt_ocean'))
    yt_oceans = dataset.createVariable('yt_ocean', 'f4',('yt_ocean'))
    st_oceans = dataset.createVariable('st_ocean', 'i4', ('st_ocean'))
    # get attributes
    attrs = [time,st_ocean,yt_ocean, xt_ocean]
    attrs_n = [times,st_oceans,yt_oceans, xt_oceans]
    i = 0
    for var in attrs:
        for j in var.ncattrs():
            attrs_n[i].setncattr(j,var.getncattr(j))
        i = i+1

    t_cya_verts = dataset.createVariable('t_cya_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_lpp_verts = dataset.createVariable('t_lpp_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_spp_verts = dataset.createVariable('t_spp_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_cya_fields = dataset.createVariable('t_cya_field', numpy.float32,('time'))
    t_lpp_fields = dataset.createVariable('t_lpp_field', numpy.float32,('time'))
    t_spp_fields = dataset.createVariable('t_spp_field', numpy.float32,('time'))
    # set units

    #   write data
    t_cya_verts[:] = t_cya_vert
    t_lpp_verts[:] = t_lpp_vert
    t_spp_verts[:] = t_spp_vert
    t_cya_fields[:] = t_cya_field
    t_lpp_fields[:] = t_lpp_field
    t_spp_fields[:] = t_spp_field
   
    times[:] = time  
    xt_oceans[:] = xt_ocean
    yt_oceans[:] = yt_ocean
    st_oceans[:] = st_ocean
    
    # plotting figures
    print 'plotting figure ...'
    fig1 = pl.figure(1)
    fig1.suptitle(year, fontsize=16)
    fig1.figsize=(10,10)
    pl.subplot(3,3,1)
    pl.pcolor(t_cya_vert[0,:,:])
    pl.title('Cya January',fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,2)
    pl.pcolor(t_cya_vert[4,:,:])
    pl.title('Cya May', fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,3)
    pl.pcolor(t_cya_vert[7,:,:])
    pl.title('Cya August', fontsize=8)
    pl.colorbar()
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,4)
    pl.pcolor(t_lpp_vert[0,:,:],)
    pl.title('Lpp January', fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,5)
    pl.pcolor(t_lpp_vert[4,:,:])
    pl.title('Lpp May', fontsize=8)
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,6)
    pl.pcolor(t_lpp_vert[7,:,:])
    pl.title('Lpp August', fontsize=8)
    pl.colorbar()
    pl.subplots_adjust(hspace = .5)
    pl.subplot(3,3,7)
    pl.pcolor(t_spp_vert[0,:,:])
    pl.title('Spp January', fontsize=8)
    pl.subplot(3,3,8)
    pl.pcolor(t_spp_vert[4,:,:])
    pl.title('Spp May', fontsize=8)
    pl.subplot(3,3,9)
    pl.pcolor(t_spp_vert[7,:,:])
    pl.title('Spp August', fontsize=8)
    pl.colorbar()
    print 'saving figure to figures/plot', year, '.png ...'
    pl.savefig('figures/plot' + year + '.png')
    pl.cla()
    pl.clf()
 
   # close file
    dataset.close()
    print 'Done with year:', year


