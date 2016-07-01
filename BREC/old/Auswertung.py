for i in range(1948,2008):
    year = str(i)
    nc_r10 = netCDF4.Dataset('/silos/thomas/ModelExp/BREC/V01_R02/balt-3nm-skag-v01-r02_' + year + '/ocean_day3d.nc')
    #nc_r11 = netCDF4.Dataset('../R11/fieldsum_noyear/noyear_vert'+ year + '.nc')
    t_cya_r10 = nc_r10.variables['t_cya']
    dzt = nc_r10.variables['dzt']
    area_t = nc_r10.variables['area_t']
    #t_lpp_r10 = nc_r10.variables['t_lpp']
    #t_spp_r10 = nc_r10.variables['t_spp']
    t_cya_int = t_cya_r10*dzt*area_t
    print 'Done with round one'
    count = count + 1
