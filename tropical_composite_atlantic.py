#Trying to plot some basic 500mb data

#Importing relevant libraries
import matplotlib.lines as lines
import matplotlib.patches as mpatches
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
import scipy.ndimage as ndimage

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
query.lonlat_box(north=55, south=10, east=10, west=-110).time(datetime.utcnow())
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
mkdir_p(output_dir+'/tropical')

#Get data using siphon
for i in range(0,70):
    best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
    best_ds = best_gfs.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=45, south=5, east=0, west=-100).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Absolute_vorticity_isobaric','Convective_available_potential_energy_surface','Precipitation_rate_surface','Pressure_reduced_to_MSL_msl','Vertical_velocity_pressure_isobaric','Relative_humidity_isobaric','Temperature_isobaric','u-component_of_wind_isobaric','v-component_of_wind_isobaric','Geopotential_height_isobaric')

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
        'Pressure_reduced_to_MSL_msl': 'mslp',
        'Precipitation_rate_surface': 'prate',
        'Convective_available_potential_energy_surface': 'cape',
        'Absolute_vorticity_isobaric':'avort'})

    zH5 = data['height'].metpy.sel(vertical=500 * units.hPa)
    zH5_crs = zH5.metpy.cartopy_crs

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    #f = mpcalc.coriolis_parameter(lat)
    dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat, initstring=zH5_crs.proj4_init)
    heights = data['height'].metpy.loc[{'time': time[0], 'vertical': 500. * units.hPa}]
    mslpc = data['mslp'].metpy.unit_array.to('hPa').squeeze()
    precip = data['prate'].squeeze()
    precip = precip*3600
    CAPE = data['cape'].squeeze()
    vort = data['avort'].metpy.sel(vertical=850*units.hPa).squeeze()
    #u_geo, v_geo = mpcalc.geostrophic_wind(heights, f, dx, dy)

    data['temperature'].metpy.convert_units('degC')

    # Create the matplotlib figure and axis
    fig, ax = plt.subplots(1, 1, figsize=(40, 15), subplot_kw={'projection': zH5_crs})
    '''
    capep = ax.contourf(x, y, CAPE, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000],
    alpha = 0.5, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
    cb = fig.colorbar(capep, orientation='horizontal', fraction=0.1, pad = 0.01, shrink = 0.75, aspect = 35, ax = ax, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    cb.set_label('CAPE (J/kg)', fontsize = 24)
    '''

    vorticity = ax.contourf(x,y,vort, alpha = 0.5, cmap = 'RdPu', levels=[0.00003,0.00004,0.00005,0.00006,0.00007,0.00008,0.00009,0.0001,0.00015,0.0002,0.00025,0.0003,0.00035,0.0004,0.00045,0.0005])
    cb = fig.colorbar(vorticity, orientation='horizontal', fraction=0.1,
                    pad = 0.01, shrink = 0.75, aspect = 35, ax = ax, extendrect=False, panchor = (0.5,0.0))
                    #ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
    cb.set_label('850 hPa Absolute Vorticity (1/s)')

    pcontour = ax.contourf(x,y,precip,
                            levels = [1,2,3,4,5,10,15,20,25,30],
                            cmap='summer', alpha = 0.7)

    cbr = fig.colorbar(pcontour, orientation='vertical', pad = 0.01, aspect = 25,
                        panchor = (0.999,0.5), ax = ax, extendrect=False,
                        ticks=[1,2,3,4,5,10,15,20,25,30,35,40], shrink=0.85)
    cbr.set_label('Hourly Precipitation Rate (mm/hr)', fontsize = 24)

    #Pull wind data
    wind_slice = slice(2, -2, 2)
    u_850 = data['u'].metpy.sel(vertical=850*units.hPa).squeeze()
    v_850 = data['v'].metpy.sel(vertical=850*units.hPa).squeeze()

    u_500 = data['u'].metpy.sel(vertical=500*units.hPa).squeeze()
    v_500 = data['v'].metpy.sel(vertical=500*units.hPa).squeeze()

    u_300 = data['u'].metpy.sel(vertical=300*units.hPa).squeeze()
    v_300 = data['v'].metpy.sel(vertical=300*units.hPa).squeeze()


    u1_500 = u_500.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_500 = v_500.metpy.unit_array[wind_slice, wind_slice].to('knots')

    u1_850 = u_850.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_850 = v_850.metpy.unit_array[wind_slice, wind_slice].to('knots')

    u1_300 = u_300.metpy.unit_array[wind_slice, wind_slice].to('knots')
    v1_300 = v_300.metpy.unit_array[wind_slice, wind_slice].to('knots')

    x1,y1 = np.meshgrid(x,y)
    x1 = x1[wind_slice,wind_slice]
    y1 = y1[wind_slice,wind_slice]

    #Mask so only strongest 1/3 of winds are plotted
    wspd_300 = np.sqrt(u1_300**2.+v1_300**2.)
    wspd_300 = ndimage.gaussian_filter(wspd_300,sigma=4,order=0)
    u_300m = np.ma.masked_where(wspd_300<2./3.*np.nanmax(wspd_300), u1_300)
    v_300m = np.ma.masked_where(wspd_300<2./3.*np.nanmax(wspd_300), v1_300)
    xm2 = np.ma.masked_where(wspd_300<2./3.*np.nanmax(wspd_300),x1)
    ym2 = np.ma.masked_where(wspd_300<2./3.*np.nanmax(wspd_300),y1)

    wspd_500 = np.sqrt(u1_500**2.+v1_500**2.)
    wspd_500 = ndimage.gaussian_filter(wspd_500,sigma=4,order=0)
    u_500m = np.ma.masked_where(wspd_500<0.72*np.nanmax(wspd_500), u1_500)
    v_500m = np.ma.masked_where(wspd_500<0.72*np.nanmax(wspd_500), v1_500)
    xm = np.ma.masked_where(wspd_500<0.72*np.nanmax(wspd_500),x1)
    ym = np.ma.masked_where(wspd_500<0.72*np.nanmax(wspd_500),y1)

    wspd_850 = np.sqrt(u1_850**2.+v1_850**2.)
    wspd_850 = ndimage.gaussian_filter(wspd_850,sigma=3,order=0)
    u_850m = np.ma.masked_where(wspd_850<0.8*np.nanmax(wspd_850), u1_850)
    v_850m = np.ma.masked_where(wspd_850<0.8*np.nanmax(wspd_850), v1_850)
    xm1 = np.ma.masked_where(wspd_850<0.8*np.nanmax(wspd_850),x1)
    ym1 = np.ma.masked_where(wspd_850<0.8*np.nanmax(wspd_850),y1)

    #Plot masked 500mb and 850mb wind barbs
    jet300 = ax.barbs(xm2[wind_slice,wind_slice], ym2[wind_slice,wind_slice], u_300m[wind_slice,wind_slice], v_300m[wind_slice,wind_slice], length=8, color='b', label='300-hPa Jet Core Winds (kt)')
    jet500 = ax.barbs(xm[wind_slice,wind_slice], ym[wind_slice,wind_slice], u_500m[wind_slice,wind_slice], v_500m[wind_slice,wind_slice], length=8, color='g', label='500-hPa Jet Core Winds (kt)')
    jet850 = ax.barbs(xm1[wind_slice,wind_slice], ym1[wind_slice,wind_slice], u_850m[wind_slice,wind_slice], v_850m[wind_slice,wind_slice], length=8, color='r', label='850-hPa Jet Core Winds (kt)')


    # Plot heights and temperature as contours
    height_contour = data['height'].metpy.sel(vertical=500*units.hPa).squeeze()
    #temp_contour = data['temperature'].metpy.sel(vertical=500*units.hPa).squeeze()
    h_contour = ax.contour(x, y, height_contour, colors='k', levels=range(5400, 6000, 60), linewidths=2.5)
    h_contour.clabel(fontsize=12, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
    msl_contour = ax.contour(x, y, mslpc, colors='grey', levels=range(940,1040,4))
    msl_contour.clabel(fontsize=12, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)


    #t_contour = ax.contour(x, y, temp_contour, colors='xkcd:light red', levels=range(-50, 10, 5), alpha=0.8, linestyles='--')
    #t_contour.clabel(fontsize=8, colors='xkcd:deep blue', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

    # Add geographic features
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=cfeature.COLORS['land'])
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor=cfeature.COLORS['water'])
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='#c7c783', zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor=cfeature.COLORS['water'], edgecolor='#c7c783', zorder=0)

    # Set a title and show the plot
    ax.set_title('Tropical Composite Valid at ' + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(), fontsize=36)
    #ax.set_extent((255, 295, 10, 55), crs = zH5_crs)    # Set a title and show the plot

    # Legend
    red_line = lines.Line2D([], [], color='grey', label='MSLP (hPa)')
    black_line = lines.Line2D([], [], linestyle='solid', color='k',
                              label='500mb Geopotential Heights (m)')
    leg = plt.legend(handles=[jet300, jet500, jet850, black_line, red_line], loc=3,
                     title='Legend',framealpha=1)
    leg.set_zorder(100)

    plt.savefig(output_dir+'/tropical/tropical_composite_'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png')
    fcst_hr = str(i*3)
    print('Hour '+fcst_hr+' completed! -------------------------------')
    plt.close()
