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
# Any import of metpy will activate the accessors
import metpy.calc as mpcalc
#from metpy.testing import get_test_data
from metpy.units import units

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
best_nam = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/NAM/CONUS_12km/catalog.xml?dataset=grib/NCEP/NAM/CONUS_12km/Best')
best_ds = best_nam.datasets[0]
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
runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values

# Create new directory
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/namsvr')

for i in range(0,25):
    best_nam = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/NAM/CONUS_12km/catalog.xml?dataset=grib/NCEP/NAM/CONUS_12km/Best')
    best_ds = best_nam.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=55, south=20, east=-60, west=-130).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Vertical_velocity_pressure_isobaric','Relative_humidity_isobaric','Temperature_isobaric','u-component_of_wind_isobaric','v-component_of_wind_isobaric','Geopotential_height_isobaric','Composite_reflectivity_entire_atmosphere_single_layer','Convective_precipitation_surface_3_Hour_Accumulation', 'Convective_available_potential_energy_surface','Convective_inhibition_pressure_difference_layer')

    data = ncss.get_data(query)

    #Parse data using MetPy
    ds = xr.open_dataset(NetCDF4DataStore(data))
    data = ds.metpy.parse_cf()

    #Rename variables to useful things
    data = data.rename({
        'Vertical_velocity_pressure_isobaric': 'omega',
        'Relative_humidity_isobaric': 'relative_humidity',
        'Temperature_isobaric': 'temperature',
        'u-component_of_wind_isobaric': 'u',
        'v-component_of_wind_isobaric': 'v',
        'Geopotential_height_isobaric': 'height',
        'Composite_reflectivity_entire_atmosphere_single_layer': 'radar',
        'Convective_precipitation_surface_3_Hour_Accumulation': 'conv_precip',
        'Convective_available_potential_energy_surface': 'cape',
        'Convective_inhibition_pressure_difference_layer': 'cin'
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

    #ax = plt.axes(projection=ccrs.LambertConformal())
    #data['height'].metpy.loc[{'time': time[0], 'vertical': 500. * units.hPa}].plot(ax=ax, transform=zH5_crs)
    #ax.coastlines()
    #plt.show()

    # Select the data for this time and level
    #data_level = data.metpy.loc[{time.name: time[0], vertical.name: 500. * units.hPa}]

    # Create the matplotlib figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': zH5_crs})

    # Plot CAPE as filled contours

    CAPE = data['cape'].squeeze()
    capep = ax.contourf(x, y, CAPE, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000], alpha = 0.7, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])

    reflectivity = data['radar'].squeeze()
    refp = ax.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Blues') #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])

    # Plot wind barbs, but not all of them
    #f = mpcalc.coriolis_parameter(lat)
    #u_geo, v_geo = mpcalc.geostrophic_wind(heights, f, x, y)
    wind_slice = slice(4, -4, 4)
    u_850 = data['u'].metpy.sel(vertical=850*units.hPa).squeeze()
    v_850 = data['v'].metpy.sel(vertical=850*units.hPa).squeeze()
    #ax.barbs(x[wind_slice], y[wind_slice], u_850.metpy.unit_array[wind_slice, wind_slice].to('knots'), v_850.metpy.unit_array[wind_slice, wind_slice].to('knots'), length=6, color = '#ff0000')

    u_500 = data['u'].metpy.sel(vertical=500*units.hPa).squeeze()
    v_500 = data['v'].metpy.sel(vertical=500*units.hPa).squeeze()
    #ax.barbs(x[wind_slice], y[wind_slice], u_500.metpy.unit_array[wind_slice, wind_slice].to('knots'), v_500.metpy.unit_array[wind_slice, wind_slice].to('knots'), length=6, color = '#0000ff')


    u1_500 = u_500.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_500 = v_500.metpy.unit_array[wind_slice, wind_slice].to('knots')

    u1_850 = u_850.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_850 = v_850.metpy.unit_array[wind_slice, wind_slice].to('knots')

    x1,y1 = np.meshgrid(x,y)
    x1 = x1[wind_slice,wind_slice]
    y1 = y1[wind_slice,wind_slice]

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

    # 500-hPa wind barbs
    ax.barbs(xm[wind_slice,wind_slice], ym[wind_slice,wind_slice],
            u_500m[wind_slice,wind_slice], v_500m[wind_slice,wind_slice],
            length=6, color='blue', label='500-hPa Jet Core Winds (kt)')
    ax.barbs(xm1[wind_slice,wind_slice], ym1[wind_slice,wind_slice],
            u_850m[wind_slice,wind_slice], v_850m[wind_slice,wind_slice],
            length=6, color='green', label='500-hPa Jet Core Winds (kt)')

    # Plot heights and temperature as contours
    height_contour = data['height'].metpy.sel(vertical=500*units.hPa).squeeze()
    temp_contour = data['temperature'].metpy.sel(vertical=700*units.hPa).squeeze()
    h_contour = ax.contour(x, y, height_contour, colors='k', levels=range(5400, 6000, 60))
    h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
    #t_contour = ax.contour(x, y, temp_contour, colors='xkcd:light red', levels=range(-50, 10, 5), alpha=0.8, linestyles='--')
    #t_contour.clabel(fontsize=8, colors='xkcd:deep blue', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

    # Add geographic features
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=cfeature.COLORS['land'])
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor=cfeature.COLORS['water'])
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='#c7c783', zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor=cfeature.COLORS['water'], edgecolor='#c7c783', zorder=0)
    ax.set_extent((235, 295, 25, 50), crs = cartopy.crs.PlateCarree())    # Set a title and show the plot
    ax.set_title('NAM 12km CAPE (J/kg), Composite Reflectivity (dBZ), 850mb and 500mb Winds (kts), and 500mb Geopotential Height (m) at ' + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item())
    plt.savefig(output_dir+'/namsvr/namconuscapewind'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png')
    fcst_hr = str(i*3)
    print('Hour '+fcst_hr+' completed!')
    plt.close()