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
url='http://nomads.ncep.noaa.gov:80/dods/gefs/gefs'+mdate+'/gefs_pgrb2bp5_all_'+get_init_hr(hour)+'z'

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/GEFS')

ds = xr.open_dataset(url)
times = ds['dpt2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(15,80,0.5)
lons = np.arange(220,310,0.5)

for i in range(1,40):
    #fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    data = data.rename({
    'dpt2m':'td2m'
    })

    #vertical, = data['td2m'].metpy.coordinates('vertical')
    time = data['td2m'].metpy.time
    zH5_crs = data['td2m'].metpy.cartopy_crs
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    td = data['td2m'].sel(lat=lats,lon=lons).squeeze()-273.15

    x, y = td.metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    #### Td >0C Height####
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)
    ax1.coastlines(resolution='50m')
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax1.add_feature(cfeature.STATES.with_scale('50m'))

    sub_freezing = td.where(td>0).count(dim='ens')/31
    sub_5 = td.where(td>5).count(dim='ens')/31
    sub_10 = td.where(td>10).count(dim='ens')/31

    ax1.set_title('GEFS Probability Td >0C')
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    tdc = ax1.contourf(x,y,sub_freezing,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar = fig.colorbar(tdc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    plt.savefig(output_dir+'/GEFS/gefs_tdprobs_v2_'+dtfs+'.png')
    ax1.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_tdprobsec_v2_'+dtfs+'_.png')
    ax1.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_tdprobslocal_v2_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### Td >5C ####
    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111, projection = zH5_crs)
    ax2.coastlines(resolution='50m')
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax2.add_feature(cfeature.STATES.with_scale('50m'))

    ax2.set_title('GEFS Probability Td >5C')
    ax2.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax2.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tdc1 = ax2.contourf(x,y,sub_5,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tdc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])

    plt.savefig(output_dir+'/GEFS/gefs_td5probs_v3_'+dtfs+'.png')
    ax2.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_td5probs_ec_v3_'+dtfs+'_.png')
    ax2.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_td5probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### Td >10C ####
    fig3 = plt.figure(figsize=(15,15))
    ax3 = fig3.add_subplot(111, projection = zH5_crs)
    ax3.coastlines(resolution='50m')
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax3.add_feature(cfeature.STATES.with_scale('50m'))

    ax3.set_title('GEFS Probability Td >10C')
    ax3.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax3.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tdc1 = ax3.contourf(x,y,sub_10,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tdc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])

    plt.savefig(output_dir+'/GEFS/gefs_td10probs_v3_'+dtfs+'.png')
    ax3.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_td10probs_ec_v3_'+dtfs+'_.png')
    ax3.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_td10probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    '''
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
    '''
