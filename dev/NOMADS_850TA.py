import numpy as np
import xarray as xr
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import metpy.calc as mpcalc
from metpy.units import units
from datetime import datetime
import metpy.calc.tools


# make unique directory to store output
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

#grabbing data from NOMAdata
startTime=datetime.now()

year = startTime.year

if startTime.month <10:
    month = '0'+str(startTime.month)
else:
    month = str(startTime.month)

if startTime.day <10:
    day = '0'+str(startTime.day)
else:
    day = str(startTime.day)

if startTime.hour <10:
    hour = '0'+str(startTime.hour)
else:
    hour = str(startTime.hour)

mdate = str(year)+str(month)+str(day)

def get_init_hr(hour):
    if int(hour) <9:
        init_hour = '00'
    elif int(hour) <12:
        init_hour = '06'
    elif int(hour) <17:
        init_hour = '12'
    elif int(hour) <7:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

# Get dataset from NOMADS Server
ds = xr.open_dataset('http://nomads.ncep.noaa.gov:80/dods/gfs_0p25/gfs'+mdate+'/gfs_0p25_'+get_init_hr(hour)+'z')
init_hour = get_init_hr(hour)


# Create new directory to store output
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00' #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFS_para_TyS') #create subdirectory to store GFS output like this

#Now loop through the 120 forecast hours to make the plots
for i in range(0,120):
    #Get the data for the forecast hour of interest
    data = ds
    data_full = data.metpy.assign_crs(grid_mapping_name='latitude_longitude')
    # Select desired vars
    datavars = data_full[['hgtprs', 'tmpprs', 'ugrdprs', 'vgrdprs']]
    #data = data.isel(time=i)

    # Select time

    # Select level
    print(data_full)
    datalev = data_full.sel(lev=850)
    data = datalev.isel(time=i)
    time = data.time

    # Select lat/lon slice
    data = data.sel(lon=slice(220, 310), lat=slice(15, 65))

    #print(data)

	# Calcualte dx and dy
    #dx, dy  = mpcalc.lat_lon_grid_deltas(data.lon.values, data.lat.values)

	# Set calculation units
    T, u, v = data.tmpprs.values * units('K'), data.ugrdprs.values * units('m/s'), data.vgrdprs.values * units('m/s')
    dx, dy = mpcalc.lat_lon_grid_deltas(data.lon, data.lat)
    #x, y = xr.broadcast(y, x)
    #dx, dy = mpcalc.lat_lon_grid_deltas(y, x)
	# Calculate temperature advection
    #print(dx)
    Tadv = mpcalc.advection(T,u,v,dx=dx,dy=dy)#* units('K/sec')

	# Smooth data
    z     = ndimage.gaussian_filter(data.hgtprs, sigma=3, order=0)
    Tadv = ndimage.gaussian_filter(Tadv, sigma=3, order=0) * units('K / s')

	# Set plot units
    u     = u.to('kt')
    v     = v.to('kt')
    T     = T.to('degC')
    Tadv = Tadv.to('degC/hr') * 3

	# Set Projection of Data
    datacrs = ccrs.PlateCarree()

	# Set Projection of Plot
    plotcrs = ccrs.LambertConformal(central_latitude=[30, 60], central_longitude=-100)

	# Create new figure
    fig = plt.figure(figsize=(15, 12.5))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99, hspace=0.01, wspace=0.01)

	# Add the map and set the extent
    ax = plt.subplot(gs[0], projection=plotcrs)
    ax.set_extent([235, 290, 20, 55])

	# Add state/country boundaries to plot
    country_borders=cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries', scale='10m', facecolor='none')
    ax.add_feature(country_borders, edgecolor='black', linewidth=1.0)
    state_borders=cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lakes', scale='10m', facecolor='none')
    ax.add_feature(state_borders, edgecolor='black', linewidth=0.5)

	# Plot Height Contours
    clev = np.arange(900, 3000, 30)
    cs = ax.contour(data.lon, data.lat, z, clev, colors='black', linewidths=2, transform=datacrs)
    plt.clabel(cs, fontsize=10, inline=1, inline_spacing=10, fmt='%i', rightside_up=True, use_clabeltext=True)

	# Plot Temperature Contours
    clevtemp850 = np.arange(-20, 20, 2)
    cs2 = ax.contour(data.lon, data.lat, T, clevtemp850, colors='grey', linewidths=1.25, linestyles='--', transform=datacrs)
    plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=10, fmt='%i', rightside_up=True, use_clabeltext=True)

	# Plot Colorfill of Temperature Advection
    cint = np.arange(-8, 9)
    cf = ax.contourf(data.lon, data.lat, Tadv, cint[cint != 0], extend='both', cmap='bwr', transform=datacrs)
    cb = plt.colorbar(cf, ax=ax, pad=0, aspect=50, orientation='horizontal', extendrect=True, ticks=cint)
    cb.set_label(r'$^{\circ}$C 3h$^{-1}$', size='large')

	# Plot Wind Barbs
    ax.barbs(data.lon, data.lat, u.magnitude, v.magnitude, length=6, regrid_shape=20, pivot='middle', transform=datacrs)

	# Change datetiem64 to datetime
    valid = datetime.utcfromtimestamp(data.time.values.astype('O')/1e9)

	# Add plot headers
    plt.title(f'GFS 850mb Temperature Advection', loc='left')
    plt.title(f'Run: {valid.strftime("%a %Y-%m-%d %H:%M")} UTC\nValid: {valid.strftime("%a %Y-%m-%d %H:%M")} UTC', loc='right')



    # Export plot and close
    #plt.show()
    ######## Save the plot
    plt.savefig(output_dir+'/GFS_para_TyS/gfs_850ta_'+time.dt.strftime('%Y-%m-%d %H%M').item()+'.png',bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.close()
    plt.clf()
