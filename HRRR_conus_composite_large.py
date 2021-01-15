#Severe Composite Map Using HRRR Data For CONUS. Can be adapted for regions/states

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from netCDF4 import num2date
import numpy as np
import xarray as xr
from siphon.catalog import TDSCatalog
from datetime import datetime
import datetime as dt
from xarray.backends import NetCDF4DataStore
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import metpy.calc as mpcalc
from metpy.units import units

#This helper function will create a new directory to store data
def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise


#Get initial dummy data using siphon to get proper runtime to create run-specific directory. Probably a stupid way of doing it, but hey it works.
best_hrrr = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/catalog.xml?dataset=grib/NCEP/HRRR/CONUS_2p5km/Best')
best_ds = best_hrrr.datasets[0]
ncss = best_ds.subset()
query = ncss.query()
query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow())
query.accept('netcdf4')
query.variables('Geopotential_height_isobaric')
data = ncss.get_data(query)

#Parse data using MetPy
ds = xr.open_dataset(NetCDF4DataStore(data))
data = ds.metpy.parse_cf()
filtered_ds = data.filter_by_attrs(standard_name='forecast_reference_time').squeeze()
coord_names = list(filtered_ds.coords)

#Pull time of run initialization
runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values

# Create new directory unique to this run
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/HRRR_CONUS')

#This loop generates the plot for each hour
for i in range(0,18):
    #Pull data specific to this hour
    best_hrrr = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/catalog.xml?dataset=grib/NCEP/HRRR/CONUS_2p5km/Best')
    best_ds = best_hrrr.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=50, south=21, east=-61, west=-130).time(datetime.utcnow()+dt.timedelta(hours=1*i))
    query.accept('netcdf4')
    query.variables('Temperature_isobaric','u-component_of_wind_isobaric','v-component_of_wind_isobaric','Geopotential_height_isobaric','Composite_reflectivity_entire_atmosphere', 'Convective_available_potential_energy_surface','Convective_inhibition_pressure_difference_layer','Pressure_reduced_to_MSL_msl')
    data = ncss.get_data(query)

    #Parse data using MetPy
    ds = xr.open_dataset(NetCDF4DataStore(data))
    data = ds.metpy.parse_cf()

    #Rename variables to useful things
    data = data.rename({
        'Temperature_isobaric': 'temperature',
        'u-component_of_wind_isobaric': 'u',
        'v-component_of_wind_isobaric': 'v',
        'Geopotential_height_isobaric': 'height',
        'Composite_reflectivity_entire_atmosphere': 'radar',
        'Convective_available_potential_energy_surface': 'cape',
        'Convective_inhibition_pressure_difference_layer': 'cin',
        'Pressure_reduced_to_MSL_msl':'mslp'
    })

    #Pull 500mb gph data
    zH5 = data['height'].metpy.sel(vertical=500 * units.hPa)
    zH5_crs = zH5.metpy.cartopy_crs
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    heights = data['height'].metpy.loc[{'time': time[0], 'vertical': 500. * units.hPa}]
    data['temperature'].metpy.convert_units('degC')

    # Create the matplotlib figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(40, 20), subplot_kw={'projection': zH5_crs})

    #Plot CAPE and include colorbar
    CAPE = data['cape'].squeeze()
    capep = ax.contourf(x, y, CAPE, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000], alpha = 0.7, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
    cb = fig.colorbar(capep, orientation='vertical', pad = 0.0000001, shrink=0.75, aspect = 25, ax = ax, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    cb.set_label('CAPE (J/kg)', size='large')

    #Plot simulated reflectivity and include colorbar
    reflectivity = data['radar'].squeeze()
    refp = ax.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Blues') #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])
    cbr = fig.colorbar(refp, orientation='vertical', pad = 0.00000001, aspect = 25,panchor = (0.999,0.5), ax = ax, extendrect=False, ticks=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], shrink=0.75)
    cbr.set_label('Composite Reflectivity (dBZ)', size = 'large')

    #Pull wind data
    wind_slice = slice(7, -7, 7)
    u_850 = data['u'].metpy.sel(vertical=850*units.hPa).squeeze()
    v_850 = data['v'].metpy.sel(vertical=850*units.hPa).squeeze()

    u_500 = data['u'].metpy.sel(vertical=500*units.hPa).squeeze()
    v_500 = data['v'].metpy.sel(vertical=500*units.hPa).squeeze()


    u1_500 = u_500.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_500 = v_500.metpy.unit_array[wind_slice, wind_slice].to('knots')

    u1_850 = u_850.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_850 = v_850.metpy.unit_array[wind_slice, wind_slice].to('knots')

    x1,y1 = np.meshgrid(x,y)
    x1 = x1[wind_slice,wind_slice]
    y1 = y1[wind_slice,wind_slice]

    #Mask so only strongest 1/3 of winds are plotted
    wspd_500 = np.sqrt(u1_500**2.+v1_500**2.)
    u_500m = np.ma.masked_where(wspd_500<2./3.*np.nanmax(wspd_500), u1_500)
    v_500m = np.ma.masked_where(wspd_500<2./3.*np.nanmax(wspd_500), v1_500)
    xm = np.ma.masked_where(wspd_500<2./3.*np.nanmax(wspd_500),x1)
    ym = np.ma.masked_where(wspd_500<2./3.*np.nanmax(wspd_500),y1)

    wspd_850 = np.sqrt(u1_850**2.+v1_850**2.)
    u_850m = np.ma.masked_where(wspd_850<2./3.*np.nanmax(wspd_850), u1_850)
    v_850m = np.ma.masked_where(wspd_850<2./3.*np.nanmax(wspd_850), v1_850)
    xm1 = np.ma.masked_where(wspd_850<2./3.*np.nanmax(wspd_850),x1)
    ym1 = np.ma.masked_where(wspd_850<2./3.*np.nanmax(wspd_850),y1)

    #Plot masked 500mb and 850mb wind barbs
    ax.barbs(xm[wind_slice,wind_slice], ym[wind_slice,wind_slice], u_500m[wind_slice,wind_slice], v_500m[wind_slice,wind_slice], length=6, color='blue', label='500-hPa Jet Core Winds (kt)')
    ax.barbs(xm1[wind_slice,wind_slice], ym1[wind_slice,wind_slice], u_850m[wind_slice,wind_slice], v_850m[wind_slice,wind_slice], length=6, color='green', label='850-hPa Jet Core Winds (kt)')

    # Plot heights and temperature as contours
    height_contour = data['height'].metpy.sel(vertical=500*units.hPa).squeeze()
    temp_contour = data['temperature'].metpy.sel(vertical=700*units.hPa).squeeze()
    mslpc = data['mslp'].metpy.unit_array.to('hPa').squeeze()
    h_contour = ax.contour(x, y, mslpc, colors='k', levels=range(940,1040,4))
    h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

    # Add geographic features
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=cfeature.COLORS['land'])
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor=cfeature.COLORS['water'])
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='#c7c783', zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor=cfeature.COLORS['water'], edgecolor='#c7c783', zorder=0)
    ax.set_extent((240, 288, 25, 48), crs = cartopy.crs.PlateCarree())

    # Set a title and show the plot
    ax.set_title('HRRR Severe Composite Valid at ' + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(), fontsize=36)
    plt.savefig(output_dir+'/HRRR_CONUS/hrrr_conus_composite_large'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png')
    fcst_hr = str(i*1)
    print('Hour '+fcst_hr+' completed!')
    plt.close()
