from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
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

nc = netCDF4.Dataset('auswertung1850.nc')

lat = nc.variables['xt_ocean'][:]
lon = nc.variables['yt_ocean'][:]
time = nc.variables['time'][:]
temp = nc.variables['temp_vert'][:]

BA_LON = [8.5,25]
BA_LAT = [53,60.75]
LON,LAT = meshgrid(lon,lat)
print temp[0].shape, lon.shape, lat.shape
cp = (.75,.75,.75)

map = Basemap(projection='merc', \
            resolution= 'h', \
            lat_0=54, \
            lon_0=20, \
            lat_1=60, \
            llcrnrlon=BA_LON[0], \
            llcrnrlat=BA_LAT[0], \
            urcrnrlon=BA_LON[1], \
            area_thresh=50,\
            urcrnrlat=BA_LAT[1])

h = pl.figure()
map.drawcoastlines()
map.drawstates()
map.drawcountries()
map.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # can use HTML names or codes for colors
map.drawcounties()
mer = arange(9,34,4)
merpl = map.drawmeridians(mer,labels=[0,0,0,1],dashes=[3,1],color=cp)
par = arange(50,70,2)
parpl = map.drawparallels(par,labels=[0,1,0,0],dashes=[3,1],color=cp)
x,y=map(*meshgrid(lat,lon))
cont = map.contourf(x,y,temp[0],cmap=pl.cm.jet)
#cb = map.colorbar(temp,"bottom", size="5%", pad="2%")
axcol = h.add_axes([.78, .13, .025, .35])
colbar = pl.colorbar(cont,cax=axcol)
colbar.set_label('Temp')
pl.draw()
pl.show()
