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

# year = input("Check sources and sinks of year:")
error_list = []
flux_sed = []
flux_3d = []
flux_surf = []
for i in range(1850,2008):
    next_year = str(i+1)
    year = str(i)
    nc_ocean_1852 = netCDF4.Dataset('../files/ocean_day3d/auswertung' + str(year)+ '.nc')
    nc_ocean_1853 = netCDF4.Dataset('../files/ocean_day3d/auswertung'+ str(next_year) +'.nc')
    nc_flux3d = netCDF4.Dataset('../files/flux3d/auswertung_flux_'+ year + '.nc')
    nc_flux_sed = netCDF4.Dataset('../files/flux_sed/auswertung_flux_' + year + '.nc')
    nc_flux_surf = netCDF4.Dataset('../files/flux_surf/auswertung_flux_surf_' + year + '.nc')
    nc_ocean_trps = netCDF4.Dataset('../files/ocean_trps/auswertung_ocean_trps_' + year + '.nc')
    nc_add_ocean_trps = netCDF4.Dataset('../files/ocean_trps/auswertung_add_ocean_trps' + year + '.nc')


    # Load OCEAN DATA
    variable_list_ocean_1852 = []
    variable_list_ocean_1853 = []
    variable_list_sed = []
    variable_list_surf = []
    variable_list_flux3d = []


    density = 1035
    sum_1852 = 0.0
    sum_1853 = 0.0
    sum_flux_3d = 0.0
    sum_flux_sed = 0.0
    sum_flux_surf = 0.0
    sum_ocean_trps = 0.0
    sum_add_ocean_trps = 0.0

    for i in nc_ocean_1852.variables:
        variable_list_ocean_1852.append(nc_ocean_1852.variables[str(i)])
        variable_list_ocean_1853.append(nc_ocean_1853.variables[str(i)])
    print '-----------------------------------------------------------'
    for i in [16,17,18,19,20,22,24,25]:
        sum_1852 += variable_list_ocean_1852[i][0]
        sum_1853 += variable_list_ocean_1853[i][0]

    print variable_list_ocean_1852[i].name,variable_list_ocean_1852[i][0],';', variable_list_ocean_1853[i].name,variable_list_ocean_1853[i][0]
    #ACHTUNG KEINE DICHTE BEI t_sed_1

    sum_1852 = sum_1852*density
    sum_1853 = sum_1853*density

    sum_1852 += variable_list_ocean_1852[26][0]
    sum_1853 += variable_list_ocean_1853[26][0]
    print '-----------------------------------------------------------'
    print year,': ',sum_1852*14/10**9,next_year,': ', sum_1853*14/10**9, 'in kT N'
    print '-----------------------------------------------------------'
    # flux3d
    p_n2_assim = nc_flux3d.variables['p_n2_assim_cya_field']
    p_det_denit_nh4 = nc_flux3d.variables['p_det_denit_nh4_field']
    p_h2s_oxno3_sul = nc_flux3d.variables['p_h2s_oxno3_sul_field']
    p_sul_oxno3_so4 = nc_flux3d.variables['p_sul_oxno3_so4_field']
    p_poc_denit = nc_flux3d.variables['p_poc_denit_field']
    # flux_sed
    p_sed_denit_nh4 = nc_flux_sed.variables['p_sed_denit_nh4_field']
    p_nh4_nitdenit_n2 = nc_flux_sed.variables['p_nh4_nitdenit_n2_field']
    p_sed_burial = nc_flux_sed.variables['p_sed_burial_field']
    # flux_surf
    runoff_flux_t_nh4 = nc_flux_surf.variables['runoff_flux_t_nh4_field']
    runoff_flux_t_no3 = nc_flux_surf.variables['runoff_flux_t_no3_field']
    runoff_flux_t_det = nc_flux_surf.variables['runoff_flux_t_det_field']
    dep_wet_t_nh4 = nc_flux_surf.variables['dep_wet_t_nh4_field']
    dep_wet_t_no3 = nc_flux_surf.variables['dep_wet_t_no3_field']
    diffusive_flux_t_no3 = nc_flux_surf.variables['diffusive_flux_t_no3_field']
    # ocean_trps
    t_no3_xflux_adv = nc_ocean_trps.variables['t_no3_xflux_adv_vyx']
    t_no3_xflux_dif = nc_ocean_trps.variables['t_no3_xflux_dif_vyx']
    t_nh4_xflux_adv = nc_ocean_trps.variables['t_nh4_xflux_adv_vyx']
    t_nh4_xflux_dif = nc_ocean_trps.variables['t_nh4_xflux_dif_vyx']
    t_det_xflux_adv = nc_ocean_trps.variables['t_det_xflux_adv_vyx']
    t_det_xflux_dif = nc_ocean_trps.variables['t_det_xflux_dif_vyx']

    # add_ocean_trps
    flux_cya = nc_add_ocean_trps.variables['flux_cya_vyx']
    flux_lpp = nc_add_ocean_trps.variables['flux_lpp_vyx']
    flux_spp = nc_add_ocean_trps.variables['flux_spp_vyx']
    flux_zoo = nc_add_ocean_trps.variables['flux_zoo_vyx']
    flux_det = nc_add_ocean_trps.variables['flux_det_vyx']

    dayspermonth = [31,28.25,31,30,31,30,31,31,30,31,30,31]
    sum_sources = numpy.ones((len(p_n2_assim[:])))
    sum_sinks = numpy.ones((len(p_n2_assim[:])))

    sum_flux_n2_assim_cya = 0.0
    year_diff = sum_1853 - sum_1852

    print 'year_diff in kT N', year_diff*14/10**9
    print '-----------------------------------------------------------'

    for i in range(len(p_n2_assim)):
        sum_flux_3d += (p_n2_assim[i]-(5.3)*p_det_denit_nh4[i]-0.4*p_h2s_oxno3_sul[i]-1.2*p_sul_oxno3_so4[i]-0.8*p_poc_denit[i])*dayspermonth[i]
        sum_flux_n2_assim_cya += (p_n2_assim[i])*dayspermonth[i]
        sum_flux_sed += ((5.3)*p_sed_denit_nh4[i]+p_nh4_nitdenit_n2[i]+p_sed_burial[i])*dayspermonth[i]
    for i in range(len(runoff_flux_t_nh4)):
        sum_flux_surf += (runoff_flux_t_nh4[i]+runoff_flux_t_no3[i]+runoff_flux_t_det[i]+dep_wet_t_nh4[i]+dep_wet_t_no3[i]+diffusive_flux_t_no3[i])*dayspermonth[i]*24*60*60
    for i in range(len(t_no3_xflux_adv)):
        sum_ocean_trps += (t_no3_xflux_adv[i]+t_no3_xflux_dif[i]+t_nh4_xflux_adv[i]+t_nh4_xflux_dif[i]+t_det_xflux_adv[i]+t_det_xflux_dif[i])*dayspermonth[i]*24*60*60
    #	print t_no3_xflux_adv[i],t_nh4_xflux_adv[i],t_det_xflux_adv[i]
    for i in range(len(flux_cya)):
        sum_add_ocean_trps += (flux_det[i]+flux_cya[i]+flux_lpp[i]+flux_spp[i]+flux_zoo[i])*dayspermonth[i]*24*60*60

    #print sum_flux_n2_assim_cya
    #print 'Ausgabe:',p_n2_assim[:]*30, p_det_denit_nh4[:]*30, p_h2s_oxno3_sul[:]*30,p_sul_oxno3_so4[:]*30

    sum_flux_3d = sum_flux_3d*density
    print 'sum_flux_3d: ',sum_flux_3d, 'in kt N', (sum_flux_3d*14/10**9)
    print 'sum_flux_sed: ',sum_flux_sed,'in kt N', (sum_flux_sed*14/10**9)
    print 'sum_flux_surf: ',sum_flux_surf,'in kt N', (sum_flux_surf*14/10**9)
    print 'sum_ocean_trps: ', sum_ocean_trps, 'in kt N', (sum_ocean_trps*14/10**9)
    print 'sum_add_ocean_trps: ', sum_add_ocean_trps, 'in kt N', (sum_add_ocean_trps*14)
    total_flux= sum_flux_3d-sum_flux_sed+sum_flux_surf+sum_ocean_trps+sum_add_ocean_trps*10**9
    print '-----------------------------------------------------------'
    print 'Total Flux: ',total_flux*14/10**9,'Year Diff: ', year_diff*14/10**9, 'total_flux-year_diff ', (total_flux-year_diff)*14/10**9, (total_flux-year_diff)/total_flux

    flux_surf.append((sum_flux_surf*14/10**9))
    flux_sed.append((sum_flux_sed*14/10**9))
    flux_3d.append((sum_flux_3d*14/10**9))
    error_list.append((total_flux-year_diff)*14/10**9)

xrang = numpy.arange(1,len(error_list)+1)
xrang = xrang + 1850
a = pl.figure()
pl.plot(xrang,error_list)
pl.xlabel('Years')
pl.ylabel('Error in [kT N]')
a.savefig('../figures/mean/check_sources_border1.png')

a = pl.figure()
pl.plot(flux_sed)
pl.title('Flux sed')
a.savefig('../figures/mean/check_sources_flux_sed.png')

a = pl.figure()
pl.plot(flux_3d)
pl.title('Flux 3d')
a.savefig('../figures/mean/check_sources_flux_3d.png')

a = pl.figure()
pl.plot(flux_surf)
pl.title('Flux Surf')
a.savefig('../figures/mean/check_sources_flux_surf.png')

