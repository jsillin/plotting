import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
import matplotlib.patches as mpatches

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
    if int(hour) <8:
        init_hour = '00'
    elif int(hour) <13:
        init_hour = '06'
    elif int(hour) <19:
        init_hour = '12'
    elif int(hour) <24:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)
init_hour = get_init_hr(hour)
init_hr=init_hour
url='http://nomads.ncep.noaa.gov:80/dods/gefs/gefs'+mdate+'/gefs_pgrb2ap5_all_'+get_init_hr(hour)+'z'

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/GEFS')

ds = xr.open_dataset(url)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(15,80,0.5)
lons = np.arange(220,310,0.5)

for i in range(1,40):
    #fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    data = data.rename({
    #'absvprs':'avort',
    'hgtprs':'gph',
    'tmpprs':'temp',
    'ugrdprs':'u',
    'vgrdprs': 'v',
    'prmslmsl':'mslp',
    })

    vertical, = data['temp'].metpy.coordinates('vertical')
    time = data['temp'].metpy.time
    zH5_crs = data['temp'].metpy.cartopy_crs
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    z5 = data['gph'].sel(lev=500,lat=lats,lon=lons).squeeze()
    u8 = data['u'].sel(lev=850,lat=lats,lon=lons).squeeze()*1.94384449
    v8 = data['v'].sel(lev=850,lat=lats,lon=lons).squeeze()*1.94384449

    ws8 = ((u8**2)+(v8**2))**.5
    mslp = data['mslp'].sel(lat=lats,lon=lons).squeeze()/100
    x, y = z5.metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)

    ####500mb Height####
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)
    ax1.coastlines(resolution='50m')
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax1.add_feature(cfeature.STATES.with_scale('50m'))
    for j in range(1,32):
        h5e = z5.sel(ens=j)
        #print(np.max(h5e))
        #print(np.min(h5e))
        h5sp = ax1.contour(x,y,h5e,colors=['magenta','blueviolet','mediumblue','darkorange','firebrick'],levels=np.arange(4800,6000,300),linewidths=1)
    ax1.set_title('GEFS 500hPa Geopotential Height Spaghetti Plot')
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    pink = mpatches.Patch(color='magenta',label='4800m')
    purple = mpatches.Patch(color='blueviolet', label='5100m')
    blue = mpatches.Patch(color='mediumblue', label='5400m')
    orange = mpatches.Patch(color='darkorange',label='5700m')
    red = mpatches.Patch(color='firebrick',label='6000m')
    leg = plt.legend(handles=[pink,purple,blue,orange,red],loc=3,framealpha=1)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GEFS/gefs_h5_spaghetti_v3_'+dtfs+'.png')
    ax1.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_h5_spaghetti_ec_v1_'+dtfs+'_.png')
    plt.clf()
    plt.close()
    #### MSLP ####
    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111, projection = zH5_crs)
    ax2.coastlines(resolution='50m')
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax2.add_feature(cfeature.STATES.with_scale('50m'))
    for j in range(1,32):
        msle = mslp.sel(ens=j)
        #print(np.max(h5e))
        #print(np.min(h5e))
        h5sp = ax2.contour(x,y,msle,colors=['magenta','blueviolet','mediumblue','darkorange','firebrick'],levels=np.arange(960,1040,20),linewidths=1)
    ax2.set_title('GEFS MSLP Spaghetti Plot')
    ax2.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax2.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    pink = mpatches.Patch(color='magenta',label='960mb')
    purple = mpatches.Patch(color='blueviolet', label='980mb')
    blue = mpatches.Patch(color='mediumblue', label='1000mb')
    orange = mpatches.Patch(color='darkorange',label='1020mb')
    red = mpatches.Patch(color='firebrick',label='1040mb')
    leg = plt.legend(handles=[pink,purple,blue,orange,red],loc=3,framealpha=1)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GEFS/gefs_mslp_spaghetti_v2_'+dtfs+'.png')
    ax2.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_mslp_spaghetti_ec_v2_'+dtfs+'_.png')
    plt.clf()
    plt.close()
    ####850mb wind####
    fig3 = plt.figure(figsize=(15,15))
    ax3 = fig3.add_subplot(111, projection = zH5_crs)
    ax3.coastlines(resolution='50m')
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax3.add_feature(cfeature.STATES.with_scale('50m'))
    for j in range(1,32):
        llje = ws8.sel(ens=j)
        #print(np.max(h5e))
        #print(np.min(h5e))
        h5sp = ax3.contour(x,y,ndimage.gaussian_filter(llje,sigma=2,order=0),colors=['darkorange','red','magenta'],levels=[50,70,90],linewidths=1)
    ax3.set_title('GEFS 850mb Wind Spaghetti Plot')
    ax3.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax3.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    orange = mpatches.Patch(color='darkorange',label='50kt')
    red = mpatches.Patch(color='red', label='70kt')
    purple = mpatches.Patch(color='magenta', label='90kt')
    leg = plt.legend(handles=[orange,red,purple],loc=3,framealpha=1)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GEFS/gefs_llj_spaghetti_v1_'+dtfs+'.png')
    ax3.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_llj_spaghetti_ec_v1_'+dtfs+'_.png')
    plt.clf()
    plt.close()
