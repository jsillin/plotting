import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import xarray as xr
import metpy
from datetime import datetime
import datetime as dt
from metpy.units import units
import scipy.ndimage as ndimage
from metpy.plots import USCOUNTIES
import cartopy
from scipy.ndimage.filters import generic_filter as gf
from metpy.plots import USCOUNTIES
from metpy.plots import SkewT
import metpy.calc as mpcalc
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import allan_soundingmaps_lcl as smap


'''
This program produces soundings derived from GFS model data obtained via
the NOMADS openDAP functionality and overlays these soundings above a
map of precipitation type and MSLP to assist in the assessment of spatial
and temporal changes in thermodynamic and moisture profiles.
This code was originally written by Jack Sillin. The NAM was added by Allan Diegan.
'''

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

#grabbing data from NOMADS
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
    if int(hour) <6:
        init_hour = '00'
    elif int(hour) <12:
        init_hour = '06'
    elif int(hour) <17:
        init_hour = '12'
    elif int(hour) <22:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

url = 'http://nomads.ncep.noaa.gov:80/dods/nam/nam'+mdate+'/nam_'+get_init_hr(hour)+'z'
init_hour = get_init_hr(hour)

# Create new directory to store output
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00' #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/NAM/TYS') #create subdirectory to store GFS output like this

#This actually opens the dataset from NOMADS and parses it with MetPy
ds = xr.open_dataset(url)
init_hr = dt.datetime(year,int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time #Pull out the time dimension
init_time = ds['time'][0]

lats = np.arange(20,55, 0.113)
lons = np.arange(240,310, 0.111)

#ds = ds.sel(lat = lats, lon = lons)

#Now loop through the 36 forecast hours to make the plots
for i in range(0,84):
    #Get the data for the forecast hour of interest
    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tmpprs': 'temperature',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'rhprs':'rh',
        'capesfc':'cape',
        'ugrd10m':'u',
        'vgrd10m':'v'
    })

    #Pull out the categorical precip type arrays
    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

    cape = data['cape'].squeeze()
    u_10m = data['u'].squeeze()
    v_10m = data['v'].squeeze()
    u_10m = u_10m*1.94384449
    v_10m = v_10m*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    #This extends each ptype one gridpoint outwards to prevent a gap between
    #different ptypes
    radius = 1
    kernel = np.zeros((2*radius+1,2*radius+1))
    y1,x1 = np.ogrid[-radius:radius+1,-radius:radius+1]
    mask=x1**2+y1**2 <=radius**2
    kernel[mask]=1

    #Make the ptype arrays nicer looking
    snowc= gf(catsnow,np.max,footprint=kernel)
    icec = gf(catice,np.max,footprint=kernel)
    sleetc = gf(catsleet,np.max,footprint=kernel)
    rainc = gf(catrain,np.max,footprint=kernel)

    #Coordinate stuff
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lats, lons = xr.broadcast(y, x)
    zH5_crs = data['temperature'].metpy.cartopy_crs
    wind_slice = slice(-1,1,1)
    #Processing surface temperature data
    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)

    #Fetch reflectivity data
    reflectivity = data['radar'].squeeze()

    #Create masked arrays for each ptype
    rain = np.ma.masked_where(rainc==0,reflectivity)
    sleet = np.ma.masked_where(sleetc==0,reflectivity)
    ice = np.ma.masked_where(icec==0,reflectivity)
    snow = np.ma.masked_where(snowc==0,reflectivity)

    #Process MSLP data
    mslp = data['mslp']/100.
    mslpc = mslp.squeeze()
    mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)

    #This creates a nice-looking datetime label
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    #twb = mpcalc.wet_bulb_temperature(mslp,t2m,td2m)

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)

    for scale, axis in zip(['500k'], [ax1]):

        ax1.coastlines(resolution='10m')
        ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
        ax1.add_feature(cfeature.STATES.with_scale('10m'), linewidth=4.0)
        axis.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='gray')

    ########## PLOTTING #######################################################
    #Plot 2m 32F isotherm
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])

    #Plot labeled MSLP contours
    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1080,4),linewidths=1,alpha=0.7)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    #Define levels and colors for plotting precip types
    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
    qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
    qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
    qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']

    capep = ax1.contourf(x, y, cape, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000], alpha = 0.6, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])


    #Plot the underlying precip type shadings
    #Use try/except so that if no ice/sleet/etc is present, things don't break
    try:
        ra = ax1.contourf(x,y,rain,colors=qr_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no rain')
    try:
        sn = ax1.contourf(x,y,snow,colors=qs_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no snow')
    try:
        ip = ax1.contourf(x,y,sleet,colors=qi_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no sleet')
    try:
        zr = ax1.contourf(x,y,ice, colors=qz_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no ice')

    ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)


    #Set a title and extent for the map
    ax1.set_title('Precipitation Type, CAPE, and Selected Soundings',fontsize=16)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n RAP Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax1.set_extent((274.1,278.1,34.375,37.625))#, crs = zH5_crs)    # Set a title and show the plot

    #################### SOUNDINGS ################################
    '''
    Ok this is the fun part. I admit this is probably the number one ugliest
    way of possibly doing this. So sorry for the frustration!
    Basically, we're going to run six for loops, one for each row of soundings.
    For each row, we're going to start a little west of the eastern plot bound
    and move west by some fixed increment that will depend on the size of your
    domain (the regional and local files will have different londelts).
    At each point of interest, we're going to interpolate the temp/rh fields
    at each pressure level and stuff that info into an array for plotting.
    Because we're using metpy, we also need to tag it with units. I'm told there
    are ways of doing this unit tagging that are 1000% better than what's shown
    here, but I always seem to find a way to break those. This way works, though
    it is very far from elegant.
    I've added a 0C isotherm with a purple dashed line since precip types are
    a pressing forecast concern at the time I wrote this (Feb 2021). I've also
    truncated the soundings at 300mb since stuff above that doesn't matter for
    precip types but you can change this with the ptop variable for severe
    applications.
    One other note of some importance is that because MetPy doesn't currently
    support data coordinates in the rect function, great care is required to
    make sure that the sounding data is actually coming from the bottom of your
    skew t plots (I set mine to come from the bottom of the purple line). I've
    put in a feature request for this, so hopefully it will be easier in the
    future.
    '''

    '''
    Local-specific stuff...
    So this sounding code is similar but a little different than the original.
    I've tried to make it a little easier to create different domains because
    manually lining up the locations of the plots with the locations from which
    data is being pulled is really time-consuming. To make your own regional
    domain without checking this, follow the following formulae:
    Step 1: decide on bounds for your area of interest. You have 4 degrees of
    longitude and 3.25 degrees of latitude to work with. Changing these numbers
    will require re-calibrating the plot locations! In set_extent above, put
    your bounds in.
    Step 2: change startlon below to be your eastern bound + 0.45 degrees. Leave
    londelt alone unless you want to recalibrate everything!
    Step 3: change startlat to be your southern bound (no addition needed)
    '''

    sound_pres = data.lev
    sound_pres = sound_pres.squeeze()
    ptop=300
    startlon=81.9+0.45
    startlat=34.38
    londelt=0.76
    sound_lons = []
    sound_lats = []
    lat_delts = [.2,.7,1.2,1.75,2.25,2.8]
    r=5
    for i in range(0,r):
        lons = -startlon-(londelt*i)
        sound_lons.append(lons)

    for i in range(0,6):
        lats = startlat+lat_delts[i]
        sound_lats.append(lats)
    print(sound_lats)

    smap.plot_soundings(fig,ax1,data['temperature'],data['rh'],34.4,81.9,'local',cape=True)

    #uncomment the two lines below (rows,cols) to plot gridlines and check that
    #your sounding bases line up with the right lats/lons used to pull the data

    #rows = ax1.gridlines(ylocs=sound_lats,linewidth=2, linestyle='--', edgecolor='dimgrey',draw_labels=True)
    #cols = ax1.gridlines(xlocs=sound_lons,linewidth=2, linestyle='--', edgecolor='dimgrey',draw_labels=True)

    ######## Save the plot
    plt.savefig(output_dir+'/NAM/TYS/NAM_hrly_pytpe_tys_sound_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1,dpi=100)
    fcst_hr = str(0)
    plt.close()
    plt.clf()
