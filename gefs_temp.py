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
    'tmp2m':'t2m'
    })

    #vertical, = data['t2m'].metpy.coordinates('vertical')
    time = data['t2m'].metpy.time
    zH5_crs = data['t2m'].metpy.cartopy_crs
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    t = data['t2m'].sel(lat=lats,lon=lons).squeeze()-273.15

    x, y = t.metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    #### 0C ####
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)
    ax1.coastlines(resolution='50m')
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax1.add_feature(cfeature.STATES.with_scale('50m'))

    sub_freezing = t.where(t>0).count(dim='ens')/31
    sub_5 = t.where(t>5).count(dim='ens')/31
    sub_10 = t.where(t>10).count(dim='ens')/31
    sub_m5 = t.where(t>-5).count(dim='ens')/31
    sub_m10 = t.where(t>-10).count(dim='ens')/31

    ax1.set_title('GEFS Probability 2m Temperature >0C')
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    tc = ax1.contourf(x,y,sub_freezing,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    plt.savefig(output_dir+'/GEFS/gefs_tprobs_v2_'+dtfs+'.png')
    ax1.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_tprobsec_v2_'+dtfs+'_.png')
    ax1.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_tprobslocal_v2_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 5C ####
    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111, projection = zH5_crs)
    ax2.coastlines(resolution='50m')
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax2.add_feature(cfeature.STATES.with_scale('50m'))

    ax2.set_title('GEFS Probability 2m Temperature >5C')
    ax2.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax2.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax2.contourf(x,y,sub_5,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    cbar1.set_label('Probability of 2m Temperature >5C')

    plt.savefig(output_dir+'/GEFS/gefs_t5probs_v3_'+dtfs+'.png')
    ax2.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_t5probs_ec_v3_'+dtfs+'_.png')
    ax2.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_t5probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 10C ####
    fig3 = plt.figure(figsize=(15,15))
    ax3 = fig3.add_subplot(111, projection = zH5_crs)
    ax3.coastlines(resolution='50m')
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax3.add_feature(cfeature.STATES.with_scale('50m'))

    ax3.set_title('GEFS Probability 2m Temperature >10C')
    ax3.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax3.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax3.contourf(x,y,sub_10,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    cbar1.set_label('Probability of 2m Temperature >10C')

    plt.savefig(output_dir+'/GEFS/gefs_t10probs_v3_'+dtfs+'.png')
    ax3.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_t10probs_ec_v3_'+dtfs+'_.png')
    ax3.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_t10probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### -5C ####
    fig4 = plt.figure(figsize=(15,15))
    ax4 = fig4.add_subplot(111, projection = zH5_crs)
    ax4.coastlines(resolution='50m')
    ax4.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax4.add_feature(cfeature.STATES.with_scale('50m'))

    ax4.set_title('GEFS Probability 2m Temperature >-5C')
    ax4.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax4.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax4.contourf(x,y,sub_10,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    cbar1.set_label('Probability of 2m Temperature >-5C')

    plt.savefig(output_dir+'/GEFS/gefs_tm5probs_v3_'+dtfs+'.png')
    ax4.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_tm5probs_ec_v3_'+dtfs+'_.png')
    ax4.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_tm5probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### -10C ####
    fig5 = plt.figure(figsize=(15,15))
    ax5 = fig5.add_subplot(111, projection = zH5_crs)
    ax5.coastlines(resolution='50m')
    ax5.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax5.add_feature(cfeature.STATES.with_scale('50m'))

    ax5.set_title('GEFS Probability 2m Temperature >-10C')
    ax5.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax5.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax5.contourf(x,y,sub_10,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    cbar1.set_label('Probability of 2m Temperature >-10C')
    plt.savefig(output_dir+'/GEFS/gefs_tm10probs_v3_'+dtfs+'.png')
    ax5.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_tm10probs_ec_v3_'+dtfs+'_.png')
    ax5.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_tm10probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()
