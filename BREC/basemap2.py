from netCDF4 import Dataset
import numpy as np
my_example_nc_file = 'auswertung1850.nc'
fh = Dataset(my_example_nc_file, mode='r')
lons = fh.variables['xt_ocean'][:]
lats = fh.variables['yt_ocean'][:]
tmax = fh.variables['temp_vert'][:]
tmax_units = fh.variables['temp_vert'].units
fh.close()
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

BA_LON = [8.5,25]
BA_LAT = [53,60.75]
cp = (.75,.75,.75)

m = Basemap(projection='merc', \
            resolution= 'i', \
            lat_0=54, \
            lon_0=20, \
            lat_1=60, \
            llcrnrlon=BA_LON[0], \
            llcrnrlat=BA_LAT[0], \
            urcrnrlon=BA_LON[1], \
            area_thresh=50,\
            urcrnrlat=BA_LAT[1])

# Because our lon and lat variables are 1D,
# use meshgrid to create 2D arrays
# Not necessary if coordinates are already in 2D arrays.
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
# Plot Data
cs = m.pcolor(xi,yi,np.squeeze(tmax[0]))

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
m.fillcontinents(color='coral')

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label('Temperature in Degree')

# Add Title
plt.title('Temperature 1850')

plt.show()
