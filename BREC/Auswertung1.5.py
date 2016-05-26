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

def writeFilemean():
    location = ('files/mean/ocean_day3d_means.nc')
    dataset_mean = netCDF4.Dataset(location, 'w', format='NETCDF4')
    time = dataset_mean.createDimension('time', None)
    st_ocean = dataset_mean.createDimension('st_ocean', 134)
    yt_ocean = dataset_mean.createDimension('yt_ocean', 242)
    xt_ocean = dataset_mean.createDimension('xt_ocean', 224)
    return dataset_mean


variable_mean = []
for i in range (1850, 2010):

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
    dzt = nc_data.variables['dzt']
    area_t = nc_data.variables['area_t']
    temp = nc_data.variables['temp']
    salt = nc_data.variables['salt']
    # Tracer variables
    t_cya = nc_data.variables['t_cya']
    t_lpp = nc_data.variables['t_lpp']
    t_spp = nc_data.variables['t_spp']
    t_nh4 = nc_data.variables['t_nh4']
    t_no3 = nc_data.variables['t_no3']
    t_po4 = nc_data.variables['t_po4']
    t_zoo = nc_data.variables['t_zoo']
    t_ipw = nc_data.variables['t_ipw']
    t_det = nc_data.variables['t_det']
    t_don = nc_data.variables['t_don']
    t_sed_1 = nc_data.variables['t_sed_1']
    t_ips_1 = nc_data.variables['t_ips_1']
    dzt_data = dzt[:]
    area_t_data = area_t[:]

    # Calculation
    print 'starting calculation ...'
    t_cya_int = t_cya[:]*dzt_data*area_t_data
    t_lpp_int = t_lpp[:]*dzt_data*area_t_data
    t_spp_int = t_spp[:]*dzt_data*area_t_data
    t_nh4_int = t_nh4[:]*dzt_data*area_t_data
    t_no3_int = t_no3[:]*dzt_data*area_t_data
    t_po4_int = t_po4[:]*dzt_data*area_t_data
    t_zoo_int = t_zoo[:]*dzt_data*area_t_data
    t_ipw_int = t_ipw[:]*dzt_data*area_t_data
    t_det_int = t_det[:]*dzt_data*area_t_data
    t_don_int = t_don[:]*dzt_data*area_t_data
    temp_int =  temp[:]*dzt_data*area_t_data
    salt_int =  salt[:]*dzt_data*area_t_data
    t_sed_1_int = t_sed_1[:]*area_t_data
    t_ips_1_int = t_ips_1[:]*area_t_data


    t_cya_vert = numpy.sum(t_cya_int, axis=1)
    t_lpp_vert = numpy.sum(t_lpp_int, axis=1)
    t_spp_vert = numpy.sum(t_spp_int, axis=1)
    t_nh4_vert = numpy.sum(t_nh4_int, axis=1)
    t_no3_vert = numpy.sum(t_no3_int, axis=1)
    t_po4_vert = numpy.sum(t_po4_int, axis=1)
    t_zoo_vert = numpy.sum(t_zoo_int, axis=1)
    t_ipw_vert = numpy.sum(t_ipw_int, axis=1)
    t_det_vert = numpy.sum(t_det_int, axis=1)
    t_don_vert = numpy.sum(t_don_int, axis=1)
    temp_vert  = numpy.sum(temp_int, axis=1)
    salt_vert  = numpy.sum(salt_int, axis=1)

    
    t_nh4_field = numpy.sum(t_nh4_vert, axis=1)
    t_no3_field = numpy.sum(t_no3_vert, axis=1)
    t_po4_field = numpy.sum(t_po4_vert, axis=1)
    t_zoo_field = numpy.sum(t_zoo_vert, axis=1)
    t_ipw_field = numpy.sum(t_ipw_vert, axis=1)
    t_det_field = numpy.sum(t_det_vert, axis=1)
    t_don_field = numpy.sum(t_don_vert, axis=1)
    t_nh4_field = numpy.sum(t_nh4_field, axis=1)
    t_no3_field = numpy.sum(t_no3_field, axis=1)
    t_po4_field = numpy.sum(t_po4_field, axis=1)
    t_zoo_field = numpy.sum(t_zoo_field, axis=1)
    t_ipw_field = numpy.sum(t_ipw_field, axis=1)
    t_det_field = numpy.sum(t_det_field, axis=1)
    t_don_field = numpy.sum(t_don_field, axis=1)
    t_cya_field = numpy.sum(t_cya_vert, axis=1)
    t_cya_field = numpy.sum(t_cya_field, axis=1)
    t_spp_field = numpy.sum(t_spp_vert, axis=1)
    t_spp_field = numpy.sum(t_spp_field, axis=1)
    t_lpp_field = numpy.sum(t_lpp_vert, axis=1)
    t_lpp_field = numpy.sum(t_lpp_field, axis=1)
    t_sed_1_field = numpy.sum(t_sed_1_int, axis=1)
    t_sed_1_field = numpy.sum(t_sed_1_field, axis=1)
    t_ips_1_field = numpy.sum(t_ips_1_int, axis=1)
    t_ips_1_field = numpy.sum(t_ips_1_field, axis=1)
    temp_field = numpy.sum(temp_vert, axis=1)
    temp_field = numpy.sum(temp_field, axis=1)
    salt_field = numpy.sum(salt_vert, axis=1)
    salt_field = numpy.sum(salt_field, axis=1)


    #	mean of field
    variable_means = numpy.ones((14))
    variable_means[0] = numpy.mean(t_no3_field)
    variable_means[1] = numpy.mean(t_po4_field)
    variable_means[2] = numpy.mean(t_zoo_field)
    variable_means[3] = numpy.mean(t_ipw_field)
    variable_means[4] = numpy.mean(t_det_field)
    variable_means[5] = numpy.mean(t_don_field)
    variable_means[6] = numpy.mean(t_sed_1_field)
    variable_means[7] = numpy.mean(t_ips_1_field)
    variable_means[8] = numpy.mean(t_cya_field)
    variable_means[9] = numpy.mean(t_lpp_field)
    variable_means[10] = numpy.mean(t_spp_field)
    variable_means[11] = numpy.mean(t_nh4_field)
    variable_means[12] = numpy.mean(temp_field)
    variable_means[13] = numpy.mean(salt_field)
    variable_mean.append(variable_means)
    print len(variable_mean)

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
    t_nh4_verts = dataset.createVariable('t_nh4_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_no3_verts = dataset.createVariable('t_no3_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_po4_verts = dataset.createVariable('t_po4_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_zoo_verts = dataset.createVariable('t_zoo_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_ipw_verts = dataset.createVariable('t_ipw_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_det_verts = dataset.createVariable('t_det_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    t_don_verts = dataset.createVariable('t_don_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    temp_verts = dataset.createVariable('temp_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
    salt_verts = dataset.createVariable('salt_vert', numpy.float32,('time', 'yt_ocean', 'xt_ocean'))
   

    t_cya_fields = dataset.createVariable('t_cya_field', numpy.float32,('time'))
    t_lpp_fields = dataset.createVariable('t_lpp_field', numpy.float32,('time'))
    t_spp_fields = dataset.createVariable('t_spp_field', numpy.float32,('time'))
    t_nh4_fields = dataset.createVariable('t_nh4_field', numpy.float32,('time'))
    t_no3_fields = dataset.createVariable('t_no3_field', numpy.float32,('time'))
    t_po4_fields = dataset.createVariable('t_po4_field', numpy.float32,('time'))
    t_zoo_fields = dataset.createVariable('t_zoo_field', numpy.float32,('time'))
    t_ipw_fields = dataset.createVariable('t_ipw_field', numpy.float32,('time'))
    t_det_fields = dataset.createVariable('t_det_field', numpy.float32,('time'))
    t_don_fields = dataset.createVariable('t_don_field', numpy.float32,('time'))
    t_sed_1_fields = dataset.createVariable('t_sed_1_field', numpy.float32,('time'))
    t_ips_1_fields = dataset.createVariable('t_ips_1_field', numpy.float32,('time'))
    temp_fields = dataset.createVariable('temp_field', numpy.float32,('time'))
    salt_fields = dataset.createVariable('salt_field', numpy.float32,('time'))


    # set units
    attrs_vars = [t_cya,t_lpp,t_spp,t_nh4,t_no3,t_po4,t_zoo,t_ipw,t_det,t_don, t_sed_1, t_ips_1,temp,salt]
    attrs_vert = [t_cya_verts,t_lpp_verts,t_spp_verts,t_nh4_verts,t_no3_verts,t_po4_verts,t_zoo_verts,t_ipw_verts,t_det_verts,t_don_verts,temp_verts,salt_verts]
    attrs_field = [t_cya_fields,t_lpp_fields,t_spp_fields,t_nh4_fields,t_no3_fields,t_po4_fields,t_zoo_fields,t_ipw_fields,t_det_fields,t_don_fields, t_sed_1_fields, t_ips_1_fields,temp_fields,salt_fields]
    i = 0   
    for var in attrs_vert:
        var.setncattr('units','mol * m^3/kg')
        var.setncattr('long_name',attrs_vars[i].long_name + ' vertical sum')
        var.setncattr('_FillValue', attrs_vars[i]._FillValue)
        var.setncattr('missing_value', attrs_vars[i].missing_value)
        var.setncattr('cell_methods', attrs_vars[i].cell_methods)
        var.setncattr('time_avg_info', attrs_vars[i].time_avg_info)
        var.setncattr('coordinates', attrs_vars[i].coordinates)
        i = i + 1

    i = 0
    for var in attrs_vert:
        var.setncattr('units','mol * m^3/kg')
        var.setncattr('long_name',attrs_vars[i].long_name + ' horizontal + vertical sum')
        var.setncattr('_FillValue', attrs_vars[i]._FillValue)
        var.setncattr('missing_value', attrs_vars[i].missing_value)
        var.setncattr('cell_methods', attrs_vars[i].cell_methods)
        var.setncattr('time_avg_info', attrs_vars[i].time_avg_info)
        var.setncattr('coordinates', attrs_vars[i].coordinates)
        i = i + 1

    #   write data
    t_cya_verts[:] = t_cya_vert
    t_lpp_verts[:] = t_lpp_vert
    t_spp_verts[:] = t_spp_vert
    t_nh4_verts[:] = t_nh4_vert
    t_no3_verts[:] = t_no3_vert
    t_po4_verts[:] = t_po4_vert
    t_zoo_verts[:] = t_zoo_vert
    t_ipw_verts[:] = t_ipw_vert
    t_det_verts[:] = t_det_vert
    t_don_verts[:] = t_don_vert
    temp_verts[:] = temp_vert
    salt_verts[:] = salt_vert

    t_cya_fields[:] = t_cya_field
    t_lpp_fields[:] = t_lpp_field
    t_spp_fields[:] = t_spp_field
    t_nh4_fields[:] = t_nh4_field
    t_no3_fields[:] = t_no3_field
    t_po4_fields[:] = t_po4_field
    t_zoo_fields[:] = t_zoo_field
    t_ipw_fields[:] = t_ipw_field
    t_det_fields[:] = t_det_field
    t_don_fields[:] = t_don_field
    t_sed_1_fields[:] = t_sed_1_field
    t_ips_1_fields[:] = t_ips_1_field
    temp_fields[:] = temp_field
    salt_fields[:] = salt_field

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
print len(variable_mean[:]),len(variable_mean[0][:])
a = numpy.ones((len(variable_mean[:])))
b = numpy.ones((len(variable_mean[:])))
c = numpy.ones((len(variable_mean[:])))
d = numpy.ones((len(variable_mean[:])))
e = numpy.ones((len(variable_mean[:])))
f = numpy.ones((len(variable_mean[:])))
g = numpy.ones((len(variable_mean[:])))
h = numpy.ones((len(variable_mean[:])))
k = numpy.ones((len(variable_mean[:])))
l = numpy.ones((len(variable_mean[:])))
m = numpy.ones((len(variable_mean[:])))
n = numpy.ones((len(variable_mean[:])))
o = numpy.ones((len(variable_mean[:])))
p = numpy.ones((len(variable_mean[:])))


dataset_mean = writeFilemean()
t_cya_means = dataset_mean.createVariable('t_cya_mean', numpy.float32,('time'))
t_lpp_means = dataset_mean.createVariable('t_lpp_mean', numpy.float32,('time'))
t_spp_means = dataset_mean.createVariable('t_spp_mean', numpy.float32,('time'))
t_nh4_means = dataset_mean.createVariable('t_nh4_mean', numpy.float32,('time'))
t_no3_means = dataset_mean.createVariable('t_no3_mean', numpy.float32,('time'))
t_po4_means = dataset_mean.createVariable('t_po4_mean', numpy.float32,('time'))
t_zoo_means = dataset_mean.createVariable('t_zoo_mean', numpy.float32,('time'))
t_ipw_means = dataset_mean.createVariable('t_ipw_mean', numpy.float32,('time'))
t_det_means = dataset_mean.createVariable('t_det_mean', numpy.float32,('time'))
t_don_means = dataset_mean.createVariable('t_don_mean', numpy.float32,('time'))
t_sed_1_means = dataset_mean.createVariable('t_sed_1_mean', numpy.float32,('time'))
t_ips_1_means = dataset_mean.createVariable('t_ips_1_mean', numpy.float32,('time'))
temp_means = dataset_mean.createVariable('temp_mean', numpy.float32,('time'))
salt_means = dataset_mean.createVariable('salt_mean', numpy.float32,('time'))

for i in range(len(variable_mean[:])):
    a[i] = variable_mean[i][8]
    b[i] = variable_mean[i][9]
    c[i] = variable_mean[i][10]
    d[i] = variable_mean[i][11]
    e[i] = variable_mean[i][0]
    f[i] = variable_mean[i][1]
    g[i] = variable_mean[i][2]
    h[i] = variable_mean[i][3]
    k[i] = variable_mean[i][4]
    l[i] = variable_mean[i][5]
    m[i] = variable_mean[i][6]
    n[i] = variable_mean[i][7]
    o[i] = variable_mean[i][12]
    p[i] = variable_mean[i][13]
    
t_cya_means[:] = a
t_lpp_means[:] = b
t_spp_means[:] = c
t_nh4_means[:] = d
t_no3_means[:] = e
t_po4_means[:] = f
t_zoo_means[:] = g
t_ipw_means[:] = h
t_det_means[:] = k
t_don_means[:] = l
t_sed_1_means[:] = m
t_ips_1_means[:] = n
temp_means[:] = o
salt_means[:] = p


