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
from scipy.ndimage.filters import generic_filter as gf

q_levs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25]
qarlevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]
qaslevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5]
qazlevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5]
q_levs_r = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3]
qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']
qra_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00','#ffcc00','#ff9500','#ff4800','#ff2900','#ff1200','#ff0000','#cc0000','#990000','#990033','#b3003b','#ff3333','#ff6666','#ffffff']
qrs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5','#4f01f6','#7a00f5','#9e00f5','#b833ff','#d280ff','#cc00f1','#ad00cc','#820099','#4700b3']
qzr_cols = ['#ff0066','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333','#ff6666','#ff9999','#ffcccc']
qip_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc','#9933ff','#bf80ff','#e6ccff','#ffffff']

def wet_bulb(temp,dewpoint):
    tdd = temp-dewpoint
    wet_bulb = temp-((1/3)*tdd)
    return wet_bulb

def fram(ice,wet_bulb,velocity):
    ilr_p = ice
    ilr_t = (-0.0071*(wet_bulb**3))-(0.039*(wet_bulb**2))-(0.3904*wet_bulb)+0.5545
    ilr_v = (0.0014*(velocity**2))+(0.0027*velocity)+0.7574

    cond_1 = np.ma.masked_where(wet_bulb>-0.35,ice)
    cond_2 = np.ma.masked_where((wet_bulb<-0.35) & (velocity>12.),ice)
    cond_3 = np.ma.masked_where((wet_bulb<-0.35) & (velocity<=12.),ice)

    cond_1 = cond_1.filled(0)
    cond_2 = cond_2.filled(0)
    cond_3 = cond_3.filled(0)

    ilr_1 = (0.7*ilr_p)+(0.29*ilr_t)+(0.01*ilr_v)
    ilr_2 = (0.73*ilr_p)+(0.01*ilr_t)+(0.26*ilr_v)
    ilr_3 = (0.79*ilr_p)+(0.2*ilr_t)+(0.01*ilr_v)

    accretion_1 = cond_1*ilr_1
    accretion_2 = cond_2*ilr_2
    accretion_3 = cond_3*ilr_3

    total_accretion=accretion_1+accretion_2+accretion_3
    return total_accretion

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
    elif int(hour) <20:
        init_hour = '12'
    elif int(hour) <24:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

init_hour = get_init_hr(hour)
init_hr=init_hour
url='http://nomads.ncep.noaa.gov:80/dods/gefs/gefs'+mdate+'/gefs_pgrb2ap5_all_'+get_init_hr(hour)+'z'
tdurl='http://nomads.ncep.noaa.gov:80/dods/gefs/gefs'+mdate+'/gefs_pgrb2bp5_all_'+get_init_hr(hour)+'z'

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/GEFS')

ds = xr.open_dataset(url)
tds = xr.open_dataset(tdurl)

times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]
valid_times = [ds['time'][4],ds['time'][8],ds['time'][12],ds['time'][16],ds['time'][20],ds['time'][24],ds['time'][28],ds['time'][40]]

lats = np.arange(15,80,0.5)
lons = np.arange(220,310,0.5)

#Initialize ptype arrays by grabbing first hour of categorical precip
catrain = ds['crainsfc'].squeeze().sel(lat=lats,lon=lons).squeeze()
catsnow = ds['csnowsfc'].squeeze().sel(lat=lats,lon=lons).squeeze()
catsleet = ds['cicepsfc'].squeeze().sel(lat=lats,lon=lons).squeeze()
catice = ds['cfrzrsfc'].squeeze().sel(lat=lats,lon=lons).squeeze()

total_precip=ds['apcpsfc'].sel(lat=lats,lon=lons).squeeze()*.0393700787402

acc_sleet = np.ma.masked_where(catsleet==0,total_precip)
acc_ice = np.ma.masked_where(catice==0,total_precip)
acc_snow = np.ma.masked_where(catsnow==0,total_precip)
t2mi = ds['tmp2m'].sel(lat=lats,lon=lons).squeeze()-273.15
td2mi = ds['tmp2m'].sel(lat=lats,lon=lons).squeeze()-273.15
acc_sleet = acc_sleet.filled(0)+t2mi-t2mi
acc_ice = acc_ice.filled(0)
acc_snow = acc_snow.filled(0)+t2mi-t2mi

u10 = ds['ugrd10m'].sel(lat=lats,lon=lons).squeeze()*1.94384449
v10 = ds['vgrd10m'].sel(lat=lats,lon=lons).squeeze()*1.94384449
ws10 = ((u10**2)+(v10**2))**.5

acc_fram = fram(acc_ice,wet_bulb(t2mi,td2mi),ws10)
fram_sum_10 = acc_fram.isel(time=range(0,41,1)).sum(dim='time')
fram_sum_7 = acc_fram.isel(time=range(0,29,1)).sum(dim='time')
fram_sum_6 = acc_fram.isel(time=range(0,25,1)).sum(dim='time')
fram_sum_5 = acc_fram.isel(time=range(0,21,1)).sum(dim='time')
fram_sum_4 = acc_fram.isel(time=range(0,17,1)).sum(dim='time')
fram_sum_3 = acc_fram.isel(time=range(0,13,1)).sum(dim='time')
fram_sum_2 = acc_fram.isel(time=range(0,9,1)).sum(dim='time')
fram_sum_1 = acc_fram.isel(time=range(0,5,1)).sum(dim='time')

snow_sum_10 = acc_snow.isel(time=range(0,41,1)).sum(dim='time')
snow_sum_7 = acc_snow.isel(time=range(0,29,1)).sum(dim='time')
snow_sum_6 = acc_snow.isel(time=range(0,25,1)).sum(dim='time')
snow_sum_5 = acc_snow.isel(time=range(0,21,1)).sum(dim='time')
snow_sum_4 = acc_snow.isel(time=range(0,17,1)).sum(dim='time')
snow_sum_3 = acc_snow.isel(time=range(0,13,1)).sum(dim='time')
snow_sum_2 = acc_snow.isel(time=range(0,9,1)).sum(dim='time')
snow_sum_1 = acc_snow.isel(time=range(0,5,1)).sum(dim='time')

sleet_sum_10 = acc_sleet.isel(time=range(0,41,1)).sum(dim='time')
sleet_sum_7 = acc_sleet.isel(time=range(0,29,1)).sum(dim='time')
sleet_sum_6 = acc_sleet.isel(time=range(0,25,1)).sum(dim='time')
sleet_sum_5 = acc_sleet.isel(time=range(0,21,1)).sum(dim='time')
sleet_sum_4 = acc_sleet.isel(time=range(0,17,1)).sum(dim='time')
sleet_sum_3 = acc_sleet.isel(time=range(0,13,1)).sum(dim='time')
sleet_sum_2 = acc_sleet.isel(time=range(0,9,1)).sum(dim='time')
sleet_sum_1 = acc_sleet.isel(time=range(0,5,1)).sum(dim='time')

pct_over_inch_ice_1 = fram_sum_1.where(fram_sum_1>=1).count(dim='ens')/31
pct_over_inch_ice_2 = fram_sum_2.where(fram_sum_2>=1).count(dim='ens')/31
pct_over_inch_ice_3 = fram_sum_3.where(fram_sum_3>=1).count(dim='ens')/31
pct_over_inch_ice_4 = fram_sum_4.where(fram_sum_4>=1).count(dim='ens')/31
pct_over_inch_ice_5 = fram_sum_5.where(fram_sum_5>=1).count(dim='ens')/31
pct_over_inch_ice_6 = fram_sum_6.where(fram_sum_6>=1).count(dim='ens')/31
pct_over_inch_ice_7 = fram_sum_7.where(fram_sum_7>=1).count(dim='ens')/31
pct_over_inch_ice_10 = fram_sum_10.where(fram_sum_10>=1).count(dim='ens')/31

pct_over_half_inch_ice_1 = fram_sum_1.where(fram_sum_1>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_2 = fram_sum_2.where(fram_sum_2>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_3 = fram_sum_3.where(fram_sum_3>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_4 = fram_sum_4.where(fram_sum_4>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_5 = fram_sum_5.where(fram_sum_5>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_6 = fram_sum_6.where(fram_sum_6>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_7 = fram_sum_7.where(fram_sum_7>=0.5).count(dim='ens')/31
pct_over_half_inch_ice_10 = fram_sum_10.where(fram_sum_10>=0.5).count(dim='ens')/31

pct_over_quarter_inch_ice_1 = fram_sum_1.where(fram_sum_1>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_2 = fram_sum_2.where(fram_sum_2>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_3 = fram_sum_3.where(fram_sum_3>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_4 = fram_sum_4.where(fram_sum_4>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_5 = fram_sum_5.where(fram_sum_5>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_6 = fram_sum_6.where(fram_sum_6>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_7 = fram_sum_7.where(fram_sum_7>=0.25).count(dim='ens')/31
pct_over_quarter_inch_ice_10 = fram_sum_10.where(fram_sum_10>=0.25).count(dim='ens')/31

pct_over_inch_snow_1 = snow_sum_1.where(snow_sum_1>=1).count(dim='ens')/31
pct_over_inch_snow_2 = snow_sum_2.where(snow_sum_2>=1).count(dim='ens')/31
pct_over_inch_snow_3 = snow_sum_3.where(snow_sum_3>=1).count(dim='ens')/31
pct_over_inch_snow_4 = snow_sum_4.where(snow_sum_4>=1).count(dim='ens')/31
pct_over_inch_snow_5 = snow_sum_5.where(snow_sum_5>=1).count(dim='ens')/31
pct_over_inch_snow_6 = snow_sum_6.where(snow_sum_6>=1).count(dim='ens')/31
pct_over_inch_snow_7 = snow_sum_7.where(snow_sum_7>=1).count(dim='ens')/31
pct_over_inch_snow_10 = snow_sum_10.where(snow_sum_10>=1).count(dim='ens')/31

pct_over_half_inch_snow_1 = snow_sum_1.where(snow_sum_1>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_2 = snow_sum_2.where(snow_sum_2>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_3 = snow_sum_3.where(snow_sum_3>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_4 = snow_sum_4.where(snow_sum_4>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_5 = snow_sum_5.where(snow_sum_5>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_6 = snow_sum_6.where(snow_sum_6>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_7 = snow_sum_7.where(snow_sum_7>=0.5).count(dim='ens')/31
pct_over_half_inch_snow_10 = snow_sum_10.where(snow_sum_10>=0.5).count(dim='ens')/31

pct_over_quarter_inch_snow_1 = snow_sum_1.where(snow_sum_1>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_2 = snow_sum_2.where(snow_sum_2>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_3 = snow_sum_3.where(snow_sum_3>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_4 = snow_sum_4.where(snow_sum_4>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_5 = snow_sum_5.where(snow_sum_5>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_6 = snow_sum_6.where(snow_sum_6>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_7 = snow_sum_7.where(snow_sum_7>=0.25).count(dim='ens')/31
pct_over_quarter_inch_snow_10 = snow_sum_10.where(snow_sum_10>=0.25).count(dim='ens')/31

fram_sums = [fram_sum_1,fram_sum_2,fram_sum_3,fram_sum_4,fram_sum_5,fram_sum_6,fram_sum_7,fram_sum_10]
snow_sums = [snow_sum_1,snow_sum_2,snow_sum_3,snow_sum_4,snow_sum_5,snow_sum_6,snow_sum_7,snow_sum_10]
sleet_sums = [sleet_sum_1,sleet_sum_2,sleet_sum_3,sleet_sum_4,sleet_sum_5,sleet_sum_6,sleet_sum_7,sleet_sum_10]

pct_over_inchs = [pct_over_inch_ice_1,pct_over_inch_ice_2,pct_over_inch_ice_3,pct_over_inch_ice_4,pct_over_inch_ice_5,pct_over_inch_ice_6,pct_over_inch_ice_7,pct_over_inch_ice_10]
pct_over_halfs = [pct_over_half_inch_ice_1,pct_over_half_inch_ice_2,pct_over_half_inch_ice_3,pct_over_half_inch_ice_4,pct_over_half_inch_ice_5,pct_over_half_inch_ice_6,pct_over_half_inch_ice_7,pct_over_half_inch_ice_10]
pct_over_quarters = [pct_over_quarter_inch_ice_1,pct_over_quarter_inch_ice_2,pct_over_quarter_inch_ice_3,pct_over_quarter_inch_ice_4,pct_over_quarter_inch_ice_5,pct_over_quarter_inch_ice_6,pct_over_quarter_inch_ice_7,pct_over_quarter_inch_ice_10]

pct_over_inchs_snow = [pct_over_inch_snow_1,pct_over_inch_snow_2,pct_over_inch_snow_3,pct_over_inch_snow_4,pct_over_inch_snow_5,pct_over_inch_snow_6,pct_over_inch_snow_7,pct_over_inch_snow_10]
pct_over_halfs_snow = [pct_over_half_inch_snow_1,pct_over_half_inch_snow_2,pct_over_half_inch_snow_3,pct_over_half_inch_snow_4,pct_over_half_inch_snow_5,pct_over_half_inch_snow_6,pct_over_half_inch_snow_7,pct_over_half_inch_snow_10]
pct_over_quarters_snow = [pct_over_quarter_inch_snow_1,pct_over_quarter_inch_snow_2,pct_over_quarter_inch_snow_3,pct_over_quarter_inch_snow_4,pct_over_quarter_inch_snow_5,pct_over_quarter_inch_snow_6,pct_over_quarter_inch_snow_7,pct_over_quarter_inch_snow_10]

days = [1,2,3,4,5,6,7,10]
x, y = u10.metpy.coordinates('x', 'y')

for i in range(len(fram_sums)):
    fig = plt.figure(figsize=(45,15))
    axs = fig.subplots(nrows=3,ncols=10,subplot_kw={'projection': ccrs.PlateCarree()})
    for j in range(0,30):
        e_fram_sum = fram_sums[i].isel(ens=j+1)
        ax = axs.flat[j]
        ax.coastlines(resolution='50m')
        ax.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax.set_extent((265, 300, 25, 50))

        try:
            zr = ax.contourf(x,y,e_fram_sum,colors=qzr_cols,levels=qazlevs,alpha=0.7,extend='max')
        except:
            print('no ice')
        ax.set_title('Ensemble Member: '+str(j+1))
    fig.tight_layout()
    plt.suptitle('GEFS Individual Ensemble Member Total Ice Accretion (FRAM) Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/ice_accretion_stamps_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()

    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111,projection=ccrs.PlateCarree())

    for j in range(1,30):
        e_fram_sum = fram_sums[i].isel(ens=j)
        ax2.coastlines(resolution='50m')
        ax2.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax2.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax2.set_extent((265, 300, 25, 50))

        try:
            zr = ax2.contour(x,y,e_fram_sum,colors=['#ff00bf','#990073','#b30000','#ff6666'],levels=[0.1,0.5,1.0,2.0],alpha=0.7,linewidths=0.75)
        except:
            print('no ice')
    ax2.set_title('GEFS Individual Ensemble Member Total Ice Accretion (FRAM) Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/ice_accretion_spaghetti_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()

    fig3 = plt.figure(figsize=(45,15))
    axs3 = fig3.subplots(nrows=3,ncols=10,subplot_kw={'projection': ccrs.PlateCarree()})
    for j in range(0,30):
        e_snow_sum = snow_sums[i].isel(ens=j+1)
        ax = axs3.flat[j]
        ax.coastlines(resolution='50m')
        ax.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax.set_extent((265, 300, 25, 50))

        try:
            sn = ax.contourf(x,y,e_snow_sum,colors=qrs_cols,levels=qaslevs,alpha=0.7,extend='max')
        except:
            print('no snow')
        ax.set_title('Ensemble Member: '+str(j+1))
    fig3.tight_layout()
    plt.suptitle('GEFS Individual Ensemble Member Snow QPF Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/snow_accum_stamps_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()

    fig4 = plt.figure(figsize=(15,15))
    ax4 = fig4.add_subplot(111,projection=ccrs.PlateCarree())

    for j in range(1,30):
        e_snow_sum = snow_sums[i].isel(ens=j)
        ax4.coastlines(resolution='50m')
        ax4.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax4.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax4.set_extent((265, 300, 25, 50))

        try:
            sn = ax4.contour(x,y,e_snow_sum,colors=['#00adad','#0087f5','#7a00f5'],levels=[0.5,1.0,2.0],alpha=0.7,linewidths=0.75)
        except:
            print('no snow')
    ax4.set_title('GEFS Individual Ensemble Member Snow QPF Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/snow_accum_spaghetti_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()

    fig5 = plt.figure(figsize=(45,15))
    axs5 = fig5.subplots(nrows=3,ncols=10,subplot_kw={'projection': ccrs.PlateCarree()})
    for j in range(0,30):
        e_sleet_sum = sleet_sums[i].isel(ens=j+1)
        ax = axs5.flat[j]
        ax.coastlines(resolution='50m')
        ax.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax.set_extent((265, 300, 25, 50))

        try:
            ip = ax.contourf(x,y,e_sleet_sum,colors=qrs_cols,levels=qaslevs,alpha=0.7,extend='max')
        except:
            print('no sleet')
        ax.set_title('Ensemble Member: '+str(j+1))
    fig5.tight_layout()
    plt.suptitle('GEFS Individual Ensemble Member Sleet QPF Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/sleet_accum_stamps_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()

    fig6 = plt.figure(figsize=(15,15))
    ax6 = fig6.add_subplot(111,projection=ccrs.PlateCarree())

    for j in range(1,30):
        e_sleet_sum = sleet_sums[i].isel(ens=j)
        ax6.coastlines(resolution='50m')
        ax6.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax6.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax6.set_extent((265, 300, 25, 50))

        try:
            ip = ax6.contour(x,y,e_sleet_sum,colors=['#aa00ff','#660099','#9933ff'],levels=[0.5,1.0,2.0],alpha=0.7,linewidths=0.75)
        except:
            print('no sleet')
    ax6.set_title('GEFS Individual Ensemble Member Sleet QPF Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/sleet_accum_spaghetti_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()

    fig7 = plt.figure(figsize=(15,15))
    ax7 = fig7.add_subplot(111,projection=ccrs.PlateCarree())
    ax7.coastlines(resolution='50m')
    ax7.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax7.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
    ax7.set_extent((265, 300, 25, 50))
    tc = ax7.contourf(x,y,pct_over_halfs[i],levels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1],cmap='RdYlBu_r')
    cbar = fig7.colorbar(tc,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax7.set_extent((265, 300, 25, 50))
    ax7.set_title('GEFS Probability of Horizontal Ice Accretion (FRAM) Exceeding 0.5" \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/half_inch_ice_probs_d'+str(days[i])+'_.png')

    fig8 = plt.figure(figsize=(15,15))
    ax8 = fig8.add_subplot(111,projection=ccrs.PlateCarree())
    ax8.coastlines(resolution='50m')
    ax8.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax8.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
    ax8.set_extent((265, 300, 25, 50))
    tc = ax8.contourf(x,y,pct_over_quarters[i],levels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1],cmap='RdYlBu_r')
    cbar = fig8.colorbar(tc,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax8.set_extent((265, 300, 25, 50))
    ax8.set_title('GEFS Probability of Horizontal Ice Accretion (FRAM) Exceeding 0.25" \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/quarter_inch_ice_probs_d'+str(days[i])+'_.png')

    fig9 = plt.figure(figsize=(15,15))
    ax9 = fig9.add_subplot(111,projection=ccrs.PlateCarree())
    ax9.coastlines(resolution='50m')
    ax9.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax9.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
    ax9.set_extent((265, 300, 25, 50))
    tc = ax9.contourf(x,y,pct_over_inchs[i],levels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1],cmap='RdYlBu_r')
    cbar = fig9.colorbar(tc,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax9.set_extent((265, 300, 25, 50))
    ax9.set_title('GEFS Probability of Horizontal Ice Accretion (FRAM) Exceeding 1.0" \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/inch_ice_probs_d'+str(days[i])+'_.png')

    fig10 = plt.figure(figsize=(15,15))
    ax10 = fig10.add_subplot(111,projection=ccrs.PlateCarree())
    ax10.coastlines(resolution='50m')
    ax10.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax10.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
    ax10.set_extent((265, 300, 25, 50))
    tc = ax10.contourf(x,y,pct_over_quarters_snow[i],levels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1],cmap='RdYlBu_r')
    cbar = fig10.colorbar(tc,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax10.set_extent((265, 300, 25, 50))
    ax10.set_title('GEFS Probability of >0.25" QPF Falling as Snow \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/quarter_inch_snow_probs_d'+str(days[i])+'_.png')

    fig11 = plt.figure(figsize=(15,15))
    ax11 = fig11.add_subplot(111,projection=ccrs.PlateCarree())
    ax11.coastlines(resolution='50m')
    ax11.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax11.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
    ax11.set_extent((265, 300, 25, 50))
    tc = ax11.contourf(x,y,pct_over_halfs_snow[i],levels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1],cmap='RdYlBu_r')
    cbar = fig11.colorbar(tc,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax11.set_extent((265, 300, 25, 50))
    ax11.set_title('GEFS Probability of >0.5" QPF Falling as Snow \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/half_inch_snow_probs_d'+str(days[i])+'_.png')

    fig12 = plt.figure(figsize=(15,15))
    ax12 = fig12.add_subplot(111,projection=ccrs.PlateCarree())
    ax12.coastlines(resolution='50m')
    ax12.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax12.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
    ax12.set_extent((265, 300, 25, 50))
    tc = ax12.contourf(x,y,pct_over_inchs_snow[i],levels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1],cmap='RdYlBu_r')
    cbar = fig12.colorbar(tc,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax12.set_extent((265, 300, 25, 50))
    ax12.set_title('GEFS Probability of >1" QPF Falling as Snow \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/inch_snow_probs_d'+str(days[i])+'_.png')

    fig13 = plt.figure(figsize=(15,15))
    ax13 = fig13.add_subplot(111,projection=ccrs.PlateCarree())

    for j in range(1,30):
        e_sleet_sum = sleet_sums[i].isel(ens=j)
        e_snow_sum = snow_sums[i].isel(ens=j)
        e_fram_sum = fram_sums[i].isel(ens=j)

        ax13.coastlines(resolution='50m')
        ax13.add_feature(cfeature.BORDERS.with_scale('50m'))
        ax13.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.5)
        ax13.set_extent((265, 300, 25, 50))

        try:
            ip = ax13.contour(x,y,e_sleet_sum,colors=['#aa00ff','#660099','#9933ff'],levels=[0.5,1.0,2.0],alpha=0.7,linewidths=0.75)
        except:
            print('no sleet')
        try:
            sn = ax13.contour(x,y,e_snow_sum,colors=['#00adad','#0087f5','#7a00f5'],levels=[0.5,1.0,2.0],alpha=0.7,linewidths=0.75)
        except:
            print('no snow')
        try:
            zr = ax13.contour(x,y,e_fram_sum,colors=['#ff00bf','#990073','#b30000','#ff6666'],levels=[0.1,0.5,1.0,2.0],alpha=0.7,linewidths=0.75)
        except:
            print('no ice')
    ax6.set_title('GEFS Individual Ensemble Member Wintry Precipitation QPF Forecasts \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+valid_times[i].dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+valid_times[i].dt.strftime('%a %b %d %H:%MZ').item())
    plt.savefig(output_dir+'/GEFS/all_accum_spaghetti_dv1'+str(days[i])+'_.png')
    plt.clf()
    plt.close()
