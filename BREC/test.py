import sys, string 
from matplotlib import rc
import numpy
import pylab as pl
import netCDF4
import time as t
import datetime
from pylab import load, meshgrid, title, arange, show
from netcdftime import utime
import scipy.io
import matplotlib as mpl
nc_data = netCDF4.Dataset(/silos/boergel/BREC/files/mean/flux3d_means.nc)
p_no3_assim_lpp_mea = netCDF4.Dataset['p_no3_assim_lpp_mea']
p_nh4_assim_lpp_mea = netCDF4.Dataset['p_nh4_assim_lpp_mea']
p_no3_assim_spp_mea = netCDF4.Dataset['p_no3_assim_spp_mea']
p_n2_assim_cya_mea = netCDF4.Dataset['p_n2_assim_cya_mea']
p_assim_lpp_poc_mea = netCDF4.Dataset['p_assim_lpp_poc_mea']
p_assim_spp_poc_mea = netCDF4.Dataset['p_assim_spp_poc_mea']
p_assim_cya_poc_mea = netCDF4.Dataset['p_assim_cya_poc_mea']
p_lpp_resp_nh4_mea = netCDF4.Dataset['p_lpp_resp_nh4_mea']
p_spp_resp_nh4_mea = netCDF4.Dataset['p_spp_resp_nh4_mea']
p_cya_resp_nh4_mea = netCDF4.Dataset['p_cya_resp_nh4_mea']
p_zoo_resp_nh4_mea = netCDF4.Dataset['p_zoo_resp_nh4_mea']
p_lpp_graz_zoo_mea = netCDF4.Dataset['p_lpp_graz_zoo_mea']
p_spp_graz_zoo_mea = netCDF4.Dataset['p_spp_graz_zoo_mea']
p_cya_graz_zoo_mea = netCDF4.Dataset['p_cya_graz_zoo_mea']
p_lpp_mort_det_mea = netCDF4.Dataset['p_lpp_mort_det_mea']
p_spp_mort_det_mea = netCDF4.Dataset['p_spp_mort_det_mea']
p_cya_mort_det_mea = netCDF4.Dataset['p_cya_mort_det_mea']
p_zoo_mort_det_mea = netCDF4.Dataset['p_zoo_mort_det_mea']
p_poc_resp_mea = netCDF4.Dataset['p_poc_resp_mea']
p_poc_denit_mea = netCDF4.Dataset['p_poc_denit_mea']
p_poc_sulf_mea = netCDF4.Dataset['p_poc_sulf_mea']
p_nh4_nit_no3_mea = netCDF4.Dataset['p_nh4_nit_no3_mea']
p_h2s_oxo2_sul_mea = netCDF4.Dataset['p_h2s_oxo2_sul_mea']
p_h2s_oxno3_sul_mea = netCDF4.Dataset['p_h2s_oxno3_sul_mea']
