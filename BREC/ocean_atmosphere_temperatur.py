# -*- coding: utf-8 -*-
from scipy.interpolate import interp1d
import sys, string
from matplotlib import rc
import numpy as np
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
    location = ('/silos/boergel/BREC/files/atmosphere/' + year + '.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 1)
    yt_ocean = dataset.createDimension('yt_ocean', 49)
    xt_ocean = dataset.createDimension('xt_ocean', 89)
    return dataset

temp_atmosphere_list =[]
temp_ocean_list = []
temp_ocean_list_gotland = []
temp_atmos_list_gotland = []
for i in range (1851, 2010):
    year = str(i)
    variables_list = []

    dataset = writeFile(year)
    times = dataset.createVariable('time','f8',('time'))
    xt_oceans = dataset.createVariable('xt_ocean', 'f4',('xt_ocean'))
    yt_oceans = dataset.createVariable('yt_ocean', 'f4',('yt_ocean'))
    st_oceans = dataset.createVariable('st_ocean', 'i4', ('st_ocean'))
    temp_atmos = dataset.createVariable('temp', np.float32,('time','st_ocean', 'yt_ocean','xt_ocean'))

    nc_data = netCDF4.Dataset('/silod1/meteo/BalticSea/reconstructions/hiresaff_v02/MOM4_FORCING/'+year+'/tairK.mom.dta.nc')
    temp = nc_data.variables['tairK']
    lon = nc_data.variables['lon']
    lat = nc_data.variables['lat']
    temp_gotland = temp[:,:,37,100]
    temp_gotland = np.mean(temp_gotland, axis = 1)
    temp_gotland = np.mean(temp_gotland, axis = 0)
    temp_gotland = temp_gotland - 273.15

    temp = temp[:,:,23:72]
    temp = temp[:,:,:,53:142]

    print temp.shape
    temp_atmos[:] = temp
    temp_field = np.mean(temp, axis = 3)
    temp_field = np.mean(temp_field, axis = 2)
    temp_field = np.mean(temp_field, axis = 1)
    temp_field = np.mean(temp_field, axis = 0)
    temp_field = temp_field - 273.15
    temp_atmosphere_list.append(temp_field)
    temp_atmos_list_gotland.append(temp_gotland)

    nc_temp = netCDF4.Dataset('/silos/boergel/BREC/files/auswertung' + year + '.nc')
    temp_ocean = nc_temp.variables['temp_field'][:]
    temp_ocean_vert = nc_temp.variables['temp_vert'][:]
    temp_ocean_vert = temp_ocean_vert[:,59,118]
    temp_ocean_vert = np.mean(temp_ocean_vert, axis = 0)
    temp_ocean = np.mean(temp_ocean)
    temp_ocean_list.append(temp_ocean)
    temp_ocean_list_gotland.append(temp_ocean_vert)

print temp_atmosphere_list

pl.figure()
pl.plot(temp_atmosphere_list, label='atmosphere')
pl.plot(temp_ocean_list,label='ocean')
pl.legend()
pl.xlabel('Years since 1850')
pl.ylabel('Temperature in Degree')
pl.savefig('/silos/boergel/BREC/figures/ocean_atmo.png')

pl.figure()
pl.plot(temp_atmos_list_gotland,label='atmosphere')
pl.plot(temp_ocean_list_gotland, label='ocean')
pl.legend()
pl.xlabel('Years since 1850')
pl.ylabel('Temperature in Degree')
pl.title('At Gotland deep')
pl.savefig('/silos/boergel/BREC/figures/ocean_atmo_gotland_deep.png')
pl.show()
