#Trying to plot some basic 700mb data

#Importing relevant libraries

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
#plt.rcParams({'figure.max_open_warning':0})
from netCDF4 import num2date
import numpy as np
import xarray as xr
from siphon.catalog import TDSCatalog
from datetime import datetime
import datetime as dt
from xarray.backends import NetCDF4DataStore
import scipy.ndimage as ndimage

import metpy.calc as mpcalc

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
query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow())
query.accept('netcdf4')
query.variables('Geopotential_height_isobaric')

data = ncss.get_data(query)

#print(data)

#Parse data using MetPy
ds = xr.open_dataset(NetCDF4DataStore(data))
data = ds.metpy.parse_cf()

filtered_ds = data.filter_by_attrs(standard_name='forecast_reference_time').squeeze()
coord_names = list(filtered_ds.coords)
runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values

#print(runtime)

# Create new directory
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/H85')

#print(mkdir_p(output_dir+'/H85'))

for i in range(0,41):
    #Get data using siphon
    best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
    best_ds = best_gfs.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Vertical_velocity_pressure_isobaric','Relative_humidity_isobaric','Temperature_isobaric','u-component_of_wind_isobaric','v-component_of_wind_isobaric','Geopotential_height_isobaric')


    data = ncss.get_data(query)
    #print(data)

    #Parse data using MetPy
    ds = xr.open_dataset(NetCDF4DataStore(data))
    data = ds.metpy.parse_cf()

    # Set the data project (GFS is lat/lon format)
    #datacrs = ccrs.PlateCarree()

    #Rename variables to useful things
    data = data.rename({
        #'Vertical_velocity_pressure_isobaric': 'omega',
        #'Relative_humidity_isobaric': 'relative_humidity',
        'Temperature_isobaric': 'temperature',
        'u-component_of_wind_isobaric': 'u',
        'v-component_of_wind_isobaric': 'v',
        'Geopotential_height_isobaric': 'height'
    })

    zH85 = data['height'].metpy.sel(vertical=850 * units.hPa)
    zH85_crs = zH85.metpy.cartopy_crs

    #Define some coordinates
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

    #Get heights array
    heights = data['height'].metpy.loc[{'time': time[0], 'vertical': 850. * units.hPa}]
    #heights = mpcalc.smooth_n_point['heights'].metpy.sel(vertical=850*units.hPa).squeeze(), 9, 50
    heights = ndimage.gaussian_filter(heights, sigma=3.0, order=0)

    #Convert temps to degrees celsius
    data['temperature'].metpy.convert_units('degC')

    temp_ = data['temperature'].metpy.sel(vertical=850*units.hPa).squeeze()
    #temp_ = mpcalc.smooth_n_point['temperature'].metpy.sel(vertical=850*units.hPa).squeeze(), 9, 50

    u_wind = data['u'].metpy.sel(vertical=850*units.hPa).squeeze()
    v_wind = data['v'].metpy.sel(vertical=850*units.hPa).squeeze()

    uqvect, vqvect = mpcalc.q_vector(u_wind, v_wind, temp_, 850*units.hPa, dx, dy)

    # Compute the divergence of the Q-vectors calculated above
    q_div = -2*mpcalc.divergence(uqvect, vqvect, dx, dy, dim_order='yx')

    q_div = ndimage.gaussian_filter(q_div,sigma=8,order=0)


    # Create the matplotlib figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(12, 8), subplot_kw={'projection': zH85_crs})

    # Add geographic features
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=cfeature.COLORS['land'])
    #ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor=cfeature.COLORS['water'])
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='#c7c783')
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor=cfeature.COLORS['water'])#, edgecolor='#c7c783', zorder=0)

    clevs_850_tmpc = np.arange(-40, 41, 2)
    clevs_qdiv = list(range(-30, -4, 5))+list(range(5, 31, 5))
    cf = ax.contourf(lon, lat, q_div*1e18, clevs_qdiv, alpha = 0.5, cmap=plt.cm.bwr,
                 extend='both', transform=zH85_crs)
    cb = plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50, extendrect=True,
                  ticks=clevs_qdiv)
    cb.set_label('Q-Vector Div. (*10$^{18}$ m s$^{-1}$ kg$^{-1}$)')



    # Plot RH as filled contours


    #rel_hum = data['relative_humidity'].metpy.sel(vertical=700*units.hPa).squeeze()
    #rh = ax.contourf(x, y, rel_hum, levels=[70, 80, 90, 100], alpha = 0.7, colors=['#99ff00', '#00ff00', '#00cc00'])

    #vertvel = data['omega'].metpy.sel(vertical=700*units.hPa).squeeze()
    #vv = ax.contourf(x,y,vertvel, levels = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-.5], alpha = 0.6, colors = ['#ffe6f9', '#ffccf3', '#ffb3ed', '#ff99e7', '#ff80e1','#ff66db','#ff4dd5','#ff33cf','#ff1ac9','#00ffffff'])

    #Plot wind barbs, but not all of them
    #f = mpcalc.coriolis_parameter(lat)
    #print(dx)
    #u_geo, v_geo = mpcalc.geostrophic_wind(heights, f, dx, dy)
    #wind_slice = slice(10, -10, 10)
    #u_wind = data['u'].metpy.sel(vertical=700*units.hPa).squeeze()
    #v_wind = data['v'].metpy.sel(vertical=700*units.hPa).squeeze()
    #ax.barbs(x[wind_slice], y[wind_slice], u_wind.metpy.unit_array[wind_slice, wind_slice].to('knots'), v_wind.metpy.unit_array[wind_slice, wind_slice].to('knots'), length=6)

    # Plot heights and temperature as contours
    height_contour = data['height'].metpy.sel(vertical=850*units.hPa).squeeze()
    temp_contour = data['temperature'].metpy.sel(vertical=850*units.hPa).squeeze()
    h_contour = ax.contour(x, y, height_contour, colors='k', levels=range(0, 8000, 30))
    h_contour.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
    t_contour = ax.contour(x, y, temp_contour, colors='xkcd:light red', levels=range(-40, 20, 5), alpha=0.8, linestyles='--')
    t_contour.clabel(fontsize=8, colors='xkcd:deep blue', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

    # Plot 850-hPa Q-vectors, scale to get nice sized arrows

    # Plot 850-hPa Q-vectors, scale to get nice sized arrows
    # wind_slice = (slice(None, None, 5), slice(None, None, 5))
    wind_slice = slice(5,-5,5)
    ax.quiver(x[wind_slice], y[wind_slice],
        uqvect[wind_slice,wind_slice].m,
        vqvect[wind_slice,wind_slice].m,
        pivot='mid', color='black',
        scale=3e-10, scale_units='inches',
        transform=zH85_crs)


    # Set a title and show the plot
    ax.set_title('850 hPa Heights (m), Temperature (\u00B0C), and Q-Vector Divergence) at ' + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item())
    plt.savefig(output_dir+'/H85/height_temp_qvec11_'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png')
    fcst_hr = str(i*3)
    print('Hour '+fcst_hr+' completed!')
    plt.close()
