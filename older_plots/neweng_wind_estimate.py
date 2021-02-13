#Trying to plot some basic 700mb data

#Importing relevant libraries

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

#import metpy
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
best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
best_ds = best_gfs.datasets[0]
ncss = best_ds.subset()
query = ncss.query()
query.lonlat_box(north=50, south=38, east=-65, west=-85).time(datetime.utcnow())
query.accept('netcdf4')
query.variables('Geopotential_height_isobaric')

data = ncss.get_data(query)

#Parse data using MetPy
ds = xr.open_dataset(NetCDF4DataStore(data))
data = ds.metpy.parse_cf()

#time_of_run = data['reftime'][0].dt.strftime('%Y%m%d_%H%MZ').values

filtered_ds = data.filter_by_attrs(standard_name='forecast_reference_time').squeeze()
coord_names = list(filtered_ds.coords)
runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values

# Create new directory
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/H9')

#for i in range(0,41):
#Get data using siphon
best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
best_ds = best_gfs.datasets[0]
ncss = best_ds.subset()
query = ncss.query()
query.lonlat_box(north=47, south=40, east=-65, west=-75).time(datetime.utcnow()+dt.timedelta(hours=27))
query.accept('netcdf4')
query.variables('Temperature_surface', 'Vertical_velocity_pressure_isobaric','Relative_humidity_isobaric','Temperature_isobaric','u-component_of_wind_isobaric','v-component_of_wind_isobaric','Geopotential_height_isobaric')


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
    'Temperature_surface': 'sfctemp'
})

zH9 = data['height'].metpy.sel(vertical=925 * units.hPa)
zH9_crs = zH9.metpy.cartopy_crs

#Define some coordinates
vertical, = data['temperature'].metpy.coordinates('vertical')
time = data['temperature'].metpy.time
x, y = data['height'].metpy.coordinates('x', 'y')
lat, lon = xr.broadcast(y, x)
dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat, initstring=zH9_crs.proj4_init)

#Get heights array
heights = data['height'].metpy.loc[{'time': time[0], 'vertical': 925. * units.hPa}]

#Convert temps to degrees celsius
data['temperature'].metpy.convert_units('degC')
data['sfctemp'].metpy.convert_units('degC')

# Create the matplotlib figure and axis
fig, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': zH9_crs})

ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
# Plot RH as filled contours

#rel_hum = data['relative_humidity'].metpy.sel(vertical=925*units.hPa).squeeze()
#rh = ax.contourf(x, y, rel_hum, levels=[70, 80, 90, 100], alpha = 0.7, colors=['#99ff00', '#00ff00', '#00cc00'])

#vertvel = data['omega'].metpy.sel(vertical=925*units.hPa).squeeze()
#vv = ax.contourf(x,y,vertvel, levels = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-.5], alpha = 0.6, colors = ['#ffe6f9', '#ffccf3', '#ffb3ed', '#ff99e7', '#ff80e1','#ff66db','#ff4dd5','#ff33cf','#ff1ac9','#00ffffff'])

# Plot heights and temperature as contours
height_contour = data['height'].metpy.sel(vertical=925*units.hPa).squeeze()
temp_contour = data['temperature'].metpy.sel(vertical=925*units.hPa).squeeze()
sfctemp_contour = data['sfctemp'].squeeze()
#print(height_contour.dtype)

tempdiff = sfctemp_contour-temp_contour
#heights_km = height_contour.metpy.convert_units('kilometers')
#print(tempdiff.dtype)
#print(heights_km.dtype)
lapse_rate = (tempdiff/height_contour)*1000

h_contour = ax.contour(x, y, height_contour, colors='k', levels=range(450, 990, 30))
h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
t_contour = ax.contour(x, y, temp_contour, colors='xkcd:light red', levels=range(-40, 40, 5), alpha=0.8, linestyles='--')
t_contour.clabel(fontsize=8, colors='xkcd:deep blue', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

lr_levs = np.arange(-15,15,3)
lr_cf = ax.contourf(lon, lat, lapse_rate, lr_levs, alpha = 0.7, cmap=plt.cm.RdBu, transform=zH9_crs)
cbar = plt.colorbar(lr_cf, orientation='horizontal', pad=0, aspect=50)
cbar.set_label('Surface-925mb Lapse Rate (K)')

#Plot wind barbs, but not all of them
f = mpcalc.coriolis_parameter(lat)
u_geo, v_geo = mpcalc.geostrophic_wind(heights, f, dx, dy)
wind_slice = slice(2, -2, 2)
u_wind = data['u'].metpy.sel(vertical=925*units.hPa).squeeze()
v_wind = data['v'].metpy.sel(vertical=925*units.hPa).squeeze()
ax.barbs(x[wind_slice], y[wind_slice], u_wind.metpy.unit_array[wind_slice, wind_slice].to('knots'), v_wind.metpy.unit_array[wind_slice, wind_slice].to('knots'), length=6)


# Use MetPy to calculate the wind speed for colorfill plot, change units to knots from m/s
sped_H9 = mpcalc.wind_speed(u_wind, v_wind).to('kt')
print(sped_H9)
clevs_900_sped = np.arange(30, 110, 10)
h9_wind_array = sped_H9.m
lapse_rate_array = lapse_rate.m
height_array = height_contour.m
attrs=data['heights'].attrs
sfcwind_estimate = (h9_wind_array/2)+(lapse_rate_array*0.9)-(0.05*height_array)
sfc_wind_array = xr.DataArray(sfcwind_estimate,coords=[lat,lon], dims=['lat','lon'])
sfc_wind_array.attrs=attrs

wind_contour = ax.contour(x, y, sfcwind_array, colors='xkcd:blue', levels=range(30, 100, 10), alpha=0.9, linestyles='-')
wind_contour.clabel(fontsize=8, colors='b', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
#cf = ax.contourf(lon, lat, sped_H9, clevs_900_sped, alpha = 0.7, cmap=plt.cm.BuPu, transform=zH9_crs)
#plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50)

# Add geographic features
ax.add_feature(cfeature.LAND.with_scale('50m'), edgecolor='#c7c783', zorder=0)
ax.add_feature(cfeature.OCEAN.with_scale('50m'), edgecolor='#c7c783', zorder=0)
ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='#c7c783', zorder=0)
ax.add_feature(cfeature.LAKES.with_scale('50m'), edgecolor='#c7c783', zorder=0)

# Set a title and show the plot
ax.set_title('925 hPa Heights (m), Temperature (\u00B0C), Winds (kts), and Suface-925mb Lapse Rate (C/km) at ' + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item())
#plt.savefig(output_dir+'/newng_winds'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png')
#fcst_hr = str(i*3)
#print('Hour '+fcst_hr+' completed!')
plt.show()
plt.close()
