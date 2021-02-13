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
from scipy.ndimage import gaussian_filter
import metpy.calc as mpcalc
import numpy.ma as ma
from metpy.units import units
import scipy.ndimage as ndimage

startTime=datetime.now()

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


#Get data using siphon
best_hrrr = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/catalog.xml?dataset=grib/NCEP/HRRR/CONUS_2p5km/Best')
best_ds = best_hrrr.datasets[0]
ncss = best_ds.subset()
query = ncss.query()
query.lonlat_box(north=55, south=20, east=-60, west=-100).time(datetime.utcnow())
query.accept('netcdf4')
query.variables('Geopotential_height_isobaric')

data = ncss.get_data(query)

#Parse data using MetPy
ds = xr.open_dataset(NetCDF4DataStore(data))
data = ds.metpy.parse_cf()

filtered_ds = data.filter_by_attrs(standard_name='forecast_reference_time').squeeze()
coord_names = list(filtered_ds.coords)
runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values

# Create new directory
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/HRRR_CONUS')
mkdir_p(output_dir+'/HRRR_NE')

for i in range(0,18):
    best_hrrr = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/catalog.xml?dataset=grib/NCEP/HRRR/CONUS_2p5km/Best')
    best_ds = best_hrrr.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=50, south=21, east=-61, west=-130).time(datetime.utcnow()+dt.timedelta(hours=1*i))
    query.accept('netcdf4')
    query.variables('Dewpoint_temperature_height_above_ground','Temperature_height_above_ground','Temperature_isobaric','u-component_of_wind_height_above_ground','v-component_of_wind_height_above_ground','Geopotential_height_isobaric','Composite_reflectivity_entire_atmosphere', 'Convective_available_potential_energy_surface','Convective_inhibition_pressure_difference_layer','Pressure_reduced_to_MSL_msl')

    data = ncss.get_data(query)

    #Parse data using MetPy
    ds = xr.open_dataset(NetCDF4DataStore(data))
    data = ds.metpy.parse_cf()
    #Rename variables to useful things
    data = data.rename({
        'Temperature_isobaric': 'temperature',
        'u-component_of_wind_height_above_ground': 'u',
        'v-component_of_wind_height_above_ground': 'v',
        'Geopotential_height_isobaric': 'height',
        'Composite_reflectivity_entire_atmosphere': 'radar',
        'Convective_available_potential_energy_surface': 'cape',
        'Convective_inhibition_pressure_difference_layer': 'cin',
        'Pressure_reduced_to_MSL_msl':'mslp',
        'Temperature_height_above_ground':'sfc_temp',
        'Dewpoint_temperature_height_above_ground':'sfc_td'
    })
    zH5 = data['height'].metpy.sel(vertical=500 * units.hPa)
    zH5_crs = zH5.metpy.cartopy_crs

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)

    #f = mpcalc.coriolis_parameter(lat)
    #dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat, initstring=zH5_crs.proj4_init)
    heights = data['height'].metpy.loc[{'time': time[0], 'vertical': 500. * units.hPa}]
    #u_geo, v_geo = mpcalc.geostrophic_wind(heights, f, dx, dy)

    data['temperature'].metpy.convert_units('degC')
    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    #t2m = t2m.metpy.convert_units('degC')
    #print(t2m)
    #ax = plt.axes(projection=ccrs.LambertConformal())
    #data['height'].metpy.loc[{'time': time[0], 'vertical': 500. * units.hPa}].plot(ax=ax, transform=zH5_crs)
    #ax.coastlines()
    #plt.show()

    # Select the data for this time and level
    #data_level = data.metpy.loc[{time.name: time[0], vertical.name: 500. * units.hPa}]

    # Create the matplotlib figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(40, 20), subplot_kw={'projection': zH5_crs})
    tmp_2m = ax.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(40,110,5))
    cbr = fig.colorbar(tmp_2m, orientation = 'vertical', pad = 0.01, aspect = 25,
                        panchor = (0.999,0.5), ax = ax, extendrect=False, ticks = range(40,110,5), shrink = 0.80)
    cbr.set_label('2m Temperature (F)', fontsize = 16)

    td_2m = ax.contour(x,y,td2ms,cmap='YlGn',levels=range(30,75,5),linewidths=3)
    td_2m.clabel(fontsize=18,colors='g',inline=1,inline_spacing=4,fmt='%i',rightside_up=True,use_clabeltext=True)

    '''
    CAPE = data['cape'].squeeze()
    capep = ax.contourf(x, y, CAPE, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000], alpha = 0.7, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
    cb = fig.colorbar(capep, orientation='horizontal', pad = 0.0001, aspect = 30, ax = ax, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    cb.set_label('CAPE (J/kg)', size='large')


    reflectivity = data['radar'].squeeze()
    refp = ax.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Blues') #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])
    cbr = fig.colorbar(refp, orientation='horizontal', pad = 0.01, aspect = 30, ax = ax, extendrect=False, ticks=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65])
    cbr.set_label('Composite Reflectivity (dBZ)', size = 'large')
    '''


    mslpc = data['mslp'].metpy.unit_array.to('hPa').squeeze()
    h_contour = ax.contour(x, y, mslpc, colors='k', levels=range(940,1040,2),linewidths=2)
    h_contour.clabel(fontsize=16, colors='k', inline=1, inline_spacing=4, fmt='%i', rightside_up=True, use_clabeltext=True)

    wind_slice = slice(15,-15,15)
    u_10m = data['u'].metpy.sel(height_above_ground3=10.0).squeeze()
    v_10m = data['v'].metpy.sel(height_above_ground3=10.0).squeeze()

    u_10m = u_10m.metpy.unit_array.to('knots')
    v_10m = v_10m.metpy.unit_array.to('knots')

    ax.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=6)

    #t_contour = ax.contour(x, y, temp_contour, colors='xkcd:light red', levels=range(-50, 10, 5), alpha=0.8, linestyles='--')
    #t_contour.clabel(fontsize=8, colors='xkcd:deep blue', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
    # Add geographic features
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=cfeature.COLORS['land'])
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor=cfeature.COLORS['water'])
    ax.add_feature(cfeature.STATES.with_scale('50m'), linewidths = 3, edgecolor='#000000', zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor=cfeature.COLORS['water'], edgecolor='#c7c783', zorder=0)
    ax.set_extent((240, 288, 25, 48), crs = cartopy.crs.PlateCarree())
    # Set a title and show the plot
    ax.set_title('HRRR 2m Temperature (F), Dew Point (F), MSLP (hPa), and 10m wind (kts) at ' + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(), fontsize=36)
    plt.savefig(output_dir+'/HRRR_CONUS/hrrr_conus_temperatures'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png')
    fcst_hr = str(i*1)
    print('Hour '+fcst_hr+' completed!')
    plt.close()
    timeelapsed = datetime.now()-startTime
    print(timeelapsed)
