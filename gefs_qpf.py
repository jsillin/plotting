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
total_precip=ds['apcpsfc'].sel(lat=lats,lon=lons).isel(time=1).squeeze()*.0393700787402
for i in range(2,40):
    #fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    data = data.rename({
    'apcpsfc':'precip'
    })

    #vertical, = data['t2m'].metpy.coordinates('vertical')
    time = data['precip'].metpy.time
    zH5_crs = data['precip'].metpy.cartopy_crs
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    qpf = data['precip'].sel(lat=lats,lon=lons).squeeze()*.0393700787402
    x, y = qpf.metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    total_precip = total_precip+qpf
    print(np.max(total_precip))

    any = qpf.where(qpf>0.1).count(dim='ens')/31
    mod = qpf.where(qpf>0.5).count(dim='ens')/31
    mod2 = qpf.where(qpf>1).count(dim='ens')/31
    mod3 = qpf.where(qpf>1.5).count(dim='ens')/31
    hvy = qpf.where(qpf>2).count(dim='ens')/31
    hvy2 = qpf.where(qpf>5).count(dim='ens')/31

    tany = total_precip.where(total_precip>0.1).count(dim='ens')/31
    tmod = total_precip.where(total_precip>0.5).count(dim='ens')/31
    tmod2 = total_precip.where(total_precip>1).count(dim='ens')/31
    tmod3 = total_precip.where(total_precip>1.5).count(dim='ens')/31
    thvy = total_precip.where(total_precip>2).count(dim='ens')/31
    thvy2 = total_precip.where(total_precip>5).count(dim='ens')/31


    #### 6 0.1" ####
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)
    ax1.coastlines(resolution='50m')
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax1.add_feature(cfeature.STATES.with_scale('50m'))

    ax1.set_title('GEFS Probability 6-hourly QPF >0.1"')
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    tc = ax1.contourf(x,y,any,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar.set_label('Probability')

    plt.savefig(output_dir+'/GEFS/gefs_6q01probs_v3_'+dtfs+'.png')
    ax1.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_6q01probsec_v3_'+dtfs+'_.png')
    ax1.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_6q01probslocal_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 6 0.5" ####
    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111, projection = zH5_crs)
    ax2.coastlines(resolution='50m')
    ax2.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax2.add_feature(cfeature.STATES.with_scale('50m'))

    ax2.set_title('GEFS Probability of 6-hourly QPF >0.5"')
    ax2.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax2.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax2.contourf(x,y,mod,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig2.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability')

    plt.savefig(output_dir+'/GEFS/gefs_6q05probs_v3_'+dtfs+'.png')
    ax2.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_6q05probs_ec_v3_'+dtfs+'_.png')
    ax2.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_6q05probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 6 1" ####
    fig3 = plt.figure(figsize=(15,15))
    ax3 = fig3.add_subplot(111, projection = zH5_crs)
    ax3.coastlines(resolution='50m')
    ax3.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax3.add_feature(cfeature.STATES.with_scale('50m'))

    ax3.set_title('GEFS Probability 6-hourly QPF >1"')
    ax3.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax3.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax3.contourf(x,y,mod2,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig3.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability')

    plt.savefig(output_dir+'/GEFS/gefs_6q1probs_v3_'+dtfs+'.png')
    ax3.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_6q1probs_ec_v3_'+dtfs+'_.png')
    ax3.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_6q1probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 6 1.5" ####
    fig4 = plt.figure(figsize=(15,15))
    ax4 = fig4.add_subplot(111, projection = zH5_crs)
    ax4.coastlines(resolution='50m')
    ax4.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax4.add_feature(cfeature.STATES.with_scale('50m'))

    ax4.set_title('GEFS Probability 6-hourly QPF >1.5"')
    ax4.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax4.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax4.contourf(x,y,mod3,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig4.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability QPF >1.5"')

    plt.savefig(output_dir+'/GEFS/gefs_6q15probs_v3_'+dtfs+'.png')
    ax4.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_6q15probs_ec_v3_'+dtfs+'_.png')
    ax4.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_6q15probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 6 2" ####
    fig5 = plt.figure(figsize=(15,15))
    ax5 = fig5.add_subplot(111, projection = zH5_crs)
    ax5.coastlines(resolution='50m')
    ax5.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax5.add_feature(cfeature.STATES.with_scale('50m'))

    ax5.set_title('GEFS Probability 6-hourly QPF >2"')
    ax5.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax5.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax5.contourf(x,y,hvy,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability QPF >2"')
    plt.savefig(output_dir+'/GEFS/gefs_6q2probs_v3_'+dtfs+'.png')
    ax5.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_6q2probs_ec_v3_'+dtfs+'_.png')
    ax5.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_6q2probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 0.1" ####
    fig6 = plt.figure(figsize=(15,15))
    ax6 = fig6.add_subplot(111, projection = zH5_crs)
    ax6.coastlines(resolution='50m')
    ax6.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax6.add_feature(cfeature.STATES.with_scale('50m'))

    ax6.set_title('GEFS Probability Total QPF >0.1"')
    ax6.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax6.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    tc = ax6.contourf(x,y,tany,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar = fig6.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar.set_label('Probability')

    plt.savefig(output_dir+'/GEFS/gefs_q01probs_v3_'+dtfs+'.png')
    ax6.set_extent((265, 300, 25, 50))# crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_q01probsec_v3_'+dtfs+'_.png')
    ax6.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_q01probslocal_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 0.5" ####
    fig7 = plt.figure(figsize=(15,15))
    ax7 = fig7.add_subplot(111, projection = zH5_crs)
    ax7.coastlines(resolution='50m')
    ax7.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax7.add_feature(cfeature.STATES.with_scale('50m'))

    ax7.set_title('GEFS Probability of Total QPF >0.5"')
    ax7.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax7.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax7.contourf(x,y,tmod,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig7.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability')

    plt.savefig(output_dir+'/GEFS/gefs_q05probs_v3_'+dtfs+'.png')
    ax7.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_q05probs_ec_v3_'+dtfs+'_.png')
    ax7.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_q05probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 1" ####
    fig8 = plt.figure(figsize=(15,15))
    ax8 = fig8.add_subplot(111, projection = zH5_crs)
    ax8.coastlines(resolution='50m')
    ax8.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax8.add_feature(cfeature.STATES.with_scale('50m'))

    ax8.set_title('GEFS Probability Total QPF >1"')
    ax8.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax8.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax8.contourf(x,y,tmod2,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig8.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability')

    plt.savefig(output_dir+'/GEFS/gefs_q1probs_v3_'+dtfs+'.png')
    ax8.set_extent((265, 300, 25, 50))# crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_q1probs_ec_v3_'+dtfs+'_.png')
    ax8.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_q1probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 1.5" ####
    fig9 = plt.figure(figsize=(15,15))
    ax9 = fig9.add_subplot(111, projection = zH5_crs)
    ax9.coastlines(resolution='50m')
    ax9.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax9.add_feature(cfeature.STATES.with_scale('50m'))

    ax9.set_title('GEFS Probability Total QPF >1.5"')
    ax9.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax9.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax9.contourf(x,y,tmod3,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig9.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability QPF >1.5"')

    plt.savefig(output_dir+'/GEFS/gefs_q15probs_v3_'+dtfs+'.png')
    ax9.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_q15probs_ec_v3_'+dtfs+'_.png')
    ax9.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_q15probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 2" ####
    fig10 = plt.figure(figsize=(15,15))
    ax10 = fig10.add_subplot(111, projection = zH5_crs)
    ax10.coastlines(resolution='50m')
    ax10.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax10.add_feature(cfeature.STATES.with_scale('50m'))

    ax10.set_title('GEFS Probability Total QPF >2"')
    ax10.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax10.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax10.contourf(x,y,thvy,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig10.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    #cbar1.set_label('Probability QPF >2"')
    plt.savefig(output_dir+'/GEFS/gefs_q2probs_v3_'+dtfs+'.png')
    ax10.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_q2probs_ec_v3_'+dtfs+'_.png')
    ax10.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_q2probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()

    #### 5" ####
    fig11 = plt.figure(figsize=(15,15))
    ax11 = fig11.add_subplot(111, projection = zH5_crs)
    ax11.coastlines(resolution='50m')
    ax11.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax11.add_feature(cfeature.STATES.with_scale('50m'))

    ax11.set_title('GEFS Probability QPF >5"')
    ax11.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax11.set_title('\n GEFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    tc1 = ax11.contourf(x,y,thvy2,levels=np.linspace(0,1.05,20),cmap='RdYlBu_r')
    cbar1 = fig11.colorbar(tc,orientation='horizontal',ticks=[0,0.2,0.4,0.6,0.8,1])
    cbar1.set_label('Probability QPF >5"')
    plt.savefig(output_dir+'/GEFS/gefs_q5probs_v3_'+dtfs+'.png')
    ax11.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GEFS/gefs_q5probs_ec_v3_'+dtfs+'_.png')
    ax11.set_extent((287,292,42,46))
    plt.savefig(output_dir+'/GEFS/gefs_q5probs_local_v3_'+dtfs+'_.png')
    plt.clf()
    plt.close()
