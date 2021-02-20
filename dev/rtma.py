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
from math import log, exp
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import supplementary_tools as spt
import soundingmaps as smap

domainsize=input('Domain Size: ')
domainname=input('Domain Name: ')
centerlat =float(input('Center Lat: '))
centerlon =float(input('Center Lon: '))

rtma_mdate = spt.get_init_time('RTMA')[0]
rap_mdate = spt.get_init_time('RAP')[0]
rtma_init_hour = spt.get_init_time('RTMA')[1]
rap_init_hour = spt.get_init_time('RAP')[1]

rtma_url = 'http://nomads.ncep.noaa.gov:80/dods/rtma2p5/rtma2p5'+rtma_mdate+'/rtma2p5_anl_'+rtma_init_hour+'z'
rap_url = 'http://nomads.ncep.noaa.gov:80/dods/rap/rap'+rap_mdate+'/rap_'+rap_init_hour+'z'

# Create new directory
output_dir = rtma_mdate+'_'+rtma_init_hour+'00'
spt.mkdir_p(output_dir)
spt.mkdir_p(output_dir+'/RTMA')

#Parse data using MetPy
rtma_ds = xr.open_dataset(rtma_url)
rap_ds = xr.open_dataset(rap_url)
rtma_times = rtma_ds['tmp2m'].metpy.time
rap_times = rap_ds['tmp2m'].metpy.time

rtma_init_time = rtma_ds['time'][0]
rap_init_time = rap_ds['time'][0]

rtma_data = rtma_ds.metpy.parse_cf()
rap_data = rap_ds.metpy.parse_cf()

rtma_time = rtma_data['tmp2m'].metpy.time
rap_time = rap_data['tmp2m'].metpy.time
x, y = rtma_data['tmp2m'].metpy.coordinates('x', 'y')
lat, lon = xr.broadcast(y, x)
zH5_crs = rtma_data['tmp2m'].metpy.cartopy_crs

td2m = rtma_data['dpt2m'].squeeze()
t2m = rtma_data['tmp2m'].squeeze()
u10m = rtma_data['ugrd10m'].squeeze()*1.94384449
v10m = rtma_data['vgrd10m'].squeeze()*1.94384449

t2m = ((t2m - 273.15)*(9./5.))+32.
td2m = ((td2m - 273.15)*(9./5.))+32.

# RAP SOUNDING DATA #
rap_data = rap_data.isel(time=0)
prs_temps = rap_data['tmpprs']
prs_relh = rap_data['rhprs']

fig = plt.figure(figsize=(15,15))
ax1 = fig.add_subplot(111, projection = zH5_crs)
ax1.coastlines(resolution='10m')
ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
ax1.add_feature(cfeature.STATES.with_scale('10m'))
tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
tmp_con = ax1.contourf(x,y,t2m,cmap='YlGnBu_r',alpha=0.8,levels=range(-20,95,5))

tds = ax1.contour(x,y,td2m,colors=['forestgreen'],levels=[60],linewidths=2)
#wind_slice = slice(12,-12,12)

if domainsize=='regional':
    south = centerlat-6
    north = centerlat+6
    east = 360-(centerlon-7.5)
    west = 360-(centerlon+7.5)
    wind_slice = slice(12,-12,12)

elif domainsize=='local':
    south = centerlat-1.625
    north = centerlat+1.625
    east = 360-(centerlon-2)
    west = 360-(centerlon+2)
    wind_slice = slice(8,-8,8)

ax1.barbs(x[wind_slice],y[wind_slice],u10m[wind_slice,wind_slice],v10m[wind_slice,wind_slice], length=6,color='gray')

smap.plot_soundings(fig,ax1,prs_temps,prs_relh,centerlat,centerlon,domainsize,cape=True)

ax1.set_extent((west,east,south,north))

ax1.set_title('Mesoanalysis')
ax1.set_title('RTMA 10m Wind/2m Temp '+rtma_time.dt.strftime('%a %b %d %H:%MZ').item(),loc='left',fontsize=11)
ax1.set_title('RAP Soundings '+rap_init_time.dt.strftime('%a %b %d %H:%MZ').item(),loc='right',fontsize=11)

#cbr = fig.colorbar(tmp_con, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
#                    extendrect=False, ticks = range(-20,100,5), shrink=0.7)
#cbr.set_label('2m Temperature (F)', fontsize = 14)

plt.savefig(output_dir+'/RTMA/'+domainname+'_sounding_test6.png',bbox_inches='tight',pad_inches=0.1)
