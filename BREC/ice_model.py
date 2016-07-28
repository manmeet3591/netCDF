# ICE MODEL - SUM OVER ICE CONCENTRATION CN

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
import calendar

def writeFile():
    location = ('/silos/boergel/BREC/files/ice_model/ice_model_'+ year +'.nc')
    dataset = netCDF4.Dataset(location,'w',format='NETCDF4')
    time = dataset.createDimension('time', None)
    st_ocean = dataset.createDimension('st_ocean', 134)
    yt_ocean = dataset.createDimension('yt_ocean', 242)
    xt_ocean = dataset.createDimension('xt_ocean', 224)
    return dataset

def writeFilemean():
    location = ('/silos/boergel/BREC/files/mean/ice_model.nc')
    dataset_mean = netCDF4.Dataset(location, 'w', format='NETCDF4')
    time = dataset_mean.createDimension('time', None)
    st_ocean = dataset_mean.createDimension('st_ocean', 134)
    yt_ocean = dataset_mean.createDimension('yt_ocean', 242)
    xt_ocean = dataset_mean.createDimension('xt_ocean', 224)
    return dataset_mean

def xcrosscoeff(a,b):
    # calculates biggest crosscoefficient and shift
    check = 0
    for i in range(len(a)):
        corrcoff = np.corrcoef(a, np.roll(b,i-1))
        print corrcoff[0,1]
        if abs(corrcoff[0,1]) > abs(check):
            print abs(corrcoff[0,1]), "is bigger then " , check
            check = corrcoff[0,1]
            label = i
        return check, label


variable_means = []
salt_list_2 = []
salt_list_200 = []
field_ice = []
salt_full_200_list = []
ice_full_list = []

for i in range (1850,2009):
    year = str(i)
    nc_data = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_'+ year +'/ice_day.nc')
    dataset = writeFile()
    # Load variables into memory
    print 'loading variables into memory ...'
    # Coordinates
    time = nc_data.variables['time']
    st_ocean = nc_data.variables['ct']
    yt_ocean = nc_data.variables['yt']
    xt_ocean = nc_data.variables['xt']


    ice_concentration = nc_data.variables['CN'][:]
    ice_mass = nc_data.variables['MI'][:]
    ice_concentration_vert = np.sum(ice_concentration, axis = 1)
    ice_concentration_field = np.mean(ice_concentration_vert, axis = 1)
    ice_concentration_field = np.mean(ice_concentration_field, axis = 1)
    print ice_concentration_field.shape
    ice_full_list.append(ice_concentration_field)
    ice_concentration_field = np.mean(ice_concentration_field, axis = 0)
    ice_concentration_mean = np.mean(ice_concentration_vert, axis = 0)
    variable_means.append(ice_concentration_mean)
    field_ice.append(ice_concentration_field)

    # additional
    nc_raw = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    salt_2 = nc_raw.variables['salt'][:,1,23,53]
    salt = np.mean(salt_2, axis = 0)
    salt_list_2.append(salt)
    salt_200 = nc_raw.variables['salt'][:,22,23,53]
    #salt_200 = salt_200[1],salt_200[2],salt_200[3],salt_200[4],salt_200[8],salt_200[9],salt_200[10]
    salt_full_200_list.append(salt_200)
    salt_200 = np.mean(salt_200, axis = 0)
    salt_list_200.append(salt_200)

    # write data to file
    print 'writing data to file ...'
    #   create new variables
    print 'creating new variables ...'
    times = dataset.createVariable('time','f8',('time'))
    xt_oceans = dataset.createVariable('xt_ocean', 'f4',('xt_ocean'))
    yt_oceans = dataset.createVariable('yt_ocean', 'f4',('yt_ocean'))
    st_oceans = dataset.createVariable('st_ocean', 'i4', ('st_ocean'))
    times[:] = time
    xt_oceans[:] = xt_ocean
    yt_oceans[:] = yt_ocean

    #
    salt_2s = dataset.createVariable('salt_2', np.float32,('time'))
    salt_200s = dataset.createVariable('salt_200', np.float32,('time'))
    ice_concentration_verts = dataset.createVariable('CN_vert', np.float32,('time', 'yt_ocean', 'xt_ocean'))
    ice_concentration_fields = dataset.createVariable('CN_field', np.float32,('time'))
    ice_masss = dataset.createVariable('ice mass', np.float32,('time', 'yt_ocean', 'xt_ocean'))
    # write Data
    ice_concentration_verts[:] = ice_concentration_vert
    ice_concentration_fields[:] = ice_concentration_field
    salt_200s[:] = salt_200
    salt_2s[:] = salt_2
    ice_masss[:] = ice_mass
    print 'Done with year' + year
    dataset.close()

# Write mean
shaping = variable_means[1].shape
shaping = np.insert(shaping,0,len(variable_means[:]))
a = np.ones((shaping))
a = np.ma.array(a)
b = np.ones((len(salt_list_200[:])))
c = np.ones((len(variable_means[:])))
dataset_mean = writeFilemean()
for i in range (len(variable_means[:])):
    a[i,:,:] = variable_means[i][:,:] # ice concetration vert
    b[i] = salt_list_200[i] # salt at z = 22
    c[i] = field_ice[i] #ice concetration field
means = dataset_mean.createVariable('ice_concentration_mean',np.float32,('time', 'yt_ocean', 'xt_ocean'))
means_ice_fields = dataset_mean.createVariable('ice_concentration_field_mean',np.float32,('time'))
means_salt = dataset_mean.createVariable('salt_concentration_mean',np.float32,('time'))
means[:] = a
means_salt[:] = b
means_ice_fields[:] = c

## Evaluation

x = np.arange((len(field_ice)))
x = x+1850

fig = pl.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot(x,field_ice, '--r', label = 'Ice concentration')
ax2 = ax.twinx()
lns2 = ax2.plot(x,salt_list_200, '-b', label = 'Salt concentration')

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.grid()
ax.set_xlabel("Time (h)")
ax.set_ylabel(r"Ice Concentration")
ax2.set_ylabel(r"Salt Concentration Gotland")
pl.show()


dat =  np.column_stack((field_ice, salt_list_200, salt_list_2))
np.savetxt('test.out', dat, delimiter=" ")


# calculation
data = np.loadtxt('test.out')
field_ice = data[:,0]
salt_list_200 = data[:,1]

check = 0
for i in range(len(field_ice)):
    corrcoff = np.corrcoef(field_ice, np.roll(salt_list_2,i-1))
    print corrcoff[0,1]
    if abs(corrcoff[0,1]) > abs(check):
        print abs(corrcoff[0,1]), "is bigger then " , check
        check = corrcoff[0,1]
        label = i

ice_full = []
salt_full = []
for i in range(len(ice_full_list)):
    ice_full = np.append(ice_full, ice_full_list[i][:])
    salt_full = np.append(salt_full,salt_full_200_list[i][:])

check_2 = 0
for i in range(len(salt_full)):
    corrcoff_2 = np.corrcoef(ice_full, np.roll(ice_full,i-1))
    print corrcoff_2[0,1]
    if abs(corrcoff_2[0,1]) > abs(check_2):
        print abs(corrcoff_2[0,1]), "is bigger than " , check_2
        check_2 = corrcoff_2[0,1]
        label_2 = i


fig = pl.figure()
ax = fig.add_subplot(111)

lns1 = ax.plot(ice_full[0:100], '--r', label = 'Ice concentration')
ax2 = ax.twinx()
lns2 = ax2.plot(salt_full[0:100], '-b', label = 'Salt concentration')

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.grid()
ax.set_xlabel("Time (h)")
ax.set_ylabel(r"Ice Concentration")
ax2.set_ylabel(r"Salt Concentration Gotland")
pl.show()

