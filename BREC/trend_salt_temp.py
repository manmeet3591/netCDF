import sys, string
from matplotlib import rc
import numpy
import pylab as plt
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
from scipy import stats
import scipy.optimize as optimization

nc_data = netCDF4.Dataset('/silos/boergel/BREC/files/mean/mean_ocean.nc')

temp = nc_data.variables['temp_field_mea'][:]
salt = nc_data.variables['salt_field_mea'][:]

# Regression

slope, intercept, r_value, p_value, std_err = stats.linregress(temp,salt)
x0    = numpy.array([0.0, 0.0, 0.0])

def func(x, a, b, c):
    return a + b*x + c*x*x

print optimization.curve_fit(func,temp,x0)
