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

salt_list_2 = []
temp_list_2 = []
salt_list_200 = []
temp_list_200 = []
for i in range(1850,2009):
	year = str(i)
	nc_raw = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
	temp_2 = nc_raw.variables['temp'][:,2,69,118]
   	salt_2 = nc_raw.variables['salt'][:,2,69,118]
	temp = np.mean(temp_2, axis = 0)
	salt = np.mean(salt_2, axis = 0)
	salt_list_2.append(salt)
	temp_list_2.append(temp)
        temp_200 = nc_raw.variables['temp'][:,100,69,118]
        salt_200 = nc_raw.variables['salt'][:,100,69,118]
        temp_200 = np.mean(temp_200, axis = 0)
        salt_200 = np.mean(salt_200, axis = 0)
        salt_list_200.append(salt_200)
        temp_list_200.append(temp_200)

	print 'Done with' + year

print len(temp_list)

