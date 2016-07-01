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


year = '1850'
nc_ocean = netCDF4.Dataset('files/ocean_day3d/auswertung' + year + '.nc')
nc_flux3d = netCDF4.Dataset('files/flux3d/auswertung_flux_'+ year + '.nc')
nc_flux_sed = netCDF4.Dataset('files/flux_sed/auswertung_flux_' + year + '.nc')
nc_flux_surf = netCDF4.Dataset('files/flux_surf/auswertung_flux_surf_' + year + '.nc')


# Load OCEAN DATA
variable_list_ocean = []
variable_list_sed = []
variable_list_surf = []
variable_list_flux3d = []

for i in nc_ocean.variables:
    variable_list_ocean.append(nc_ocean.variables[str(i)])

for i in range(16,27):
    sum =+ variable_list_ocean[i][0]
    print sum + 'im Jahr 1850'
