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
import metpy.calc as mpcalc
import matplotlib.lines as lines
import matplotlib.patches as mpatches

### Function to make a new directory to store output files ###
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

### Generate proper NOMADS url based on the current time and date
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

### Select the proper run based on current time and when runs are available
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

url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr_parafv3/gfs'+mdate+'/gfs_0p25_1hr_parafv3_'+get_init_hr(hour)+'z'
init_hour = get_init_hr(hour)
print(url)
# Create new directory to store output files
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFSp')

#Access and parse data
ds = xr.open_dataset(url)
init_hr = dt.datetime(int(year),int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(15,70,0.25)
lons = np.arange(220,310,0.25)

forecast_hour = times[0].values
for i in range(1,121):
    data = ds.metpy.parse_cf()
    data = data.isel(time=i)
    #Rename variables to useful things
    data = data.rename({
        'absvprs':'avort',
        'hgtprs':'gph',
        'rhprs':'rh',
        'tmpprs':'temp',
        'ugrdprs':'u',
        'vgrdprs': 'v',
    })

    vertical, = data['temp'].metpy.coordinates('vertical')
    time = data['temp'].metpy.time
    zH5_crs = data['temp'].metpy.cartopy_crs
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    ### Pull dataarrays for each parameter of interest
    t7 = data['temp'].sel(lev=700.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    u7 = data['u'].sel(lev=700.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    v7 = data['v'].sel(lev=700.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    h7 = data['gph'].sel(lev=700.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    rh7 = data['rh'].sel(lev=700.0,lat=slice(15,70),lon=slice(220, 310)).squeeze()

    t8 = data['temp'].sel(lev=850.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    u8 = data['u'].sel(lev=850.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    v8 = data['v'].sel(lev=850.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    h8 = data['gph'].sel(lev=850.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    rh8 = data['rh'].sel(lev=850.0,lat=slice(15,70),lon=slice(220, 310)).squeeze()

    t9 = data['temp'].sel(lev=925.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    u9 = data['u'].sel(lev=925.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    v9 = data['v'].sel(lev=925.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    h9 = data['gph'].sel(lev=925.0,lat=slice(15, 70),lon=slice(220, 310)).squeeze()
    rh9 = data['rh'].sel(lev=925.0,lat=slice(15,70),lon=slice(220, 310)).squeeze()

    x, y = t7.metpy.coordinates('x', 'y')
    dx, dy = mpcalc.grid_deltas_from_dataarray(t7)
    lat, lon = xr.broadcast(y, x)
    wind_slice = slice(7,-7,7)

    #t7.attrs['units'] = 'K'
    #u7.attrs['units'] = v7.attrs['units'] = 'knots'
    #h7.attrs['units'] = 'm'

    #700
    u7k = u7*1.94384449
    v7k = v7*1.94384449
    t7c = t7-273.15
    t7c = ndimage.gaussian_filter(t7c,sigma=2,order=0)

    t7 = t7*units.K
    u7 = u7*units.meters/units.seconds
    v7 = v7*units.meters/units.seconds
    h7 = h7*units.m

    h7_fgen = mpcalc.frontogenesis(t7.data,u7.data,v7.data,dx,dy)
    h7_fgen = h7_fgen*1000*100*3600*3 ##convert to units of K/100km/3hrs
    h7_fgen = ndimage.gaussian_filter(h7_fgen,sigma=1,order=0)

    #850
    u8k = u8*1.94384449
    v8k = v8*1.94384449
    t8c = t8-273.15
    t8c = ndimage.gaussian_filter(t8c,sigma=2,order=0)

    t8 = t8*units.K
    u8 = u8*units.meters/units.seconds
    v8 = v8*units.meters/units.seconds
    h8 = h8*units.m

    h8_fgen = mpcalc.frontogenesis(t8.data,u8.data,v8.data,dx,dy)
    h8_fgen = h8_fgen*1000*100*3600*3 ##convert to units of K/100km/3hrs
    h8_fgen = ndimage.gaussian_filter(h8_fgen,sigma=1,order=0)

    #925
    u9k = u9*1.94384449
    v9k = v9*1.94384449
    t9c = t9-273.15
    t9c = ndimage.gaussian_filter(t9c,sigma=2,order=0)

    t9 = t9*units.K
    u9 = u9*units.meters/units.seconds
    v9 = v9*units.meters/units.seconds
    h9 = h9*units.m

    h9_fgen = mpcalc.frontogenesis(t9.data,u9.data,v9.data,dx,dy)
    h9_fgen = h9_fgen*1000*100*3600*3 ##convert to units of K/100km/3hrs
    h9_fgen = ndimage.gaussian_filter(h9_fgen,sigma=1,order=0)

    #### PLOTTING ####
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111,projection=zH5_crs)

    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax.add_feature(cfeature.STATES.with_scale('10m'))
    ax.set_extent((260, 295, 25, 50))

    rhc = ax.contourf(x,y,rh7,levels=[70,75,80,85,90,95,100],colors=['#73ea61','#5fe84a','#4be534','#38e21d','#32cb1a','#2cb517','#279e15'],alpha=0.6)
    cbar = fig.colorbar(rhc,orientation = 'horizontal', aspect = 80, ax = ax, pad = 0.01,
                        extendrect=False, ticks = [70,80,90,100])
    cbar.set_label('Relative Humidity (%)')
    ax.barbs(x[wind_slice],y[wind_slice],u7[wind_slice,wind_slice],v7[wind_slice,wind_slice], length=7)
    fc = ax.contour(x,y,h7_fgen,alpha=0.7,colors='fuchsia',levels=range(2,30,2),linewidths=3)
    ax.contour(x,y,h7,colors='dimgray',levels=range(2700,3300,30),linewidths=2)
    ax.contour(x,y,t7c,colors='b',levels=range(-40,0,5),linewidths=1.5,linestyles='dashed')
    ax.contour(x,y,t7c,colors='r',levels=range(0,20,50),linewidths=1.5,linestyles='dashed')
    ax.set_title('700mb Forecast Summary',fontsize=16)
    ax.set_title('GFS (v16) Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    pink = lines.Line2D([],[],linestyle='solid',color='fuchsia',label='Frontogenesis (K/100km/3hrs)')
    red = lines.Line2D([],[],linestyle='dashed',color='r',label='>=0C Temperatures (C)')
    blue = lines.Line2D([],[],linestyle='dashed',color='b',label='<0C Temperatures (C)')
    gray = lines.Line2D([],[],linestyle='solid',color='dimgray',label='Height (m)')
    leg = ax.legend(handles=[pink,gray,red,blue],framealpha=1,loc=4)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GFSp/ec_h7fgen_'+dtfs+'.png')
    ax.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFSp/ne_h7fgen_'+dtfs+'.png')

    #### Plot 2
    fig2 = plt.figure(figsize=(15,15))
    ax2 = fig2.add_subplot(111,projection=zH5_crs)

    ax2.coastlines(resolution='10m')
    ax2.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax2.add_feature(cfeature.STATES.with_scale('10m'))
    ax2.set_extent((260, 295, 25, 50))

    rhc = ax2.contourf(x,y,rh8,levels=[70,75,80,85,90,95,100],colors=['#73ea61','#5fe84a','#4be534','#38e21d','#32cb1a','#2cb517','#279e15'],alpha=0.6)
    cbar = fig2.colorbar(rhc,orientation = 'horizontal', aspect = 80, ax = ax2, pad = 0.01,
                        extendrect=False, ticks = [70,80,90,100])
    cbar.set_label('Relative Humidity (%)')
    ax2.barbs(x[wind_slice],y[wind_slice],u8[wind_slice,wind_slice],v8[wind_slice,wind_slice], length=7)
    fc = ax2.contour(x,y,h8_fgen,alpha=0.7,colors='fuchsia',levels=range(2,30,2),linewidths=3)
    ax2.contour(x,y,h8,colors='dimgray',levels=range(1200,1800,30),linewidths=2)
    ax2.contour(x,y,t8c,colors='b',levels=range(-40,0,5),linewidths=1.5,linestyles='dashed')
    ax2.contour(x,y,t8c,colors='r',levels=range(0,40,50),linewidths=1.5,linestyles='dashed')
    ax2.set_title('850mb Forecast Summary',fontsize=16)
    ax2.set_title('GFS (v16) Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax2.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    pink = lines.Line2D([],[],linestyle='solid',color='fuchsia',label='Frontogenesis (K/100km/3hrs)')
    red = lines.Line2D([],[],linestyle='dashed',color='r',label='>=0C Temperatures (C)')
    blue = lines.Line2D([],[],linestyle='dashed',color='b',label='<0C Temperatures (C)')
    gray = lines.Line2D([],[],linestyle='solid',color='dimgray',label='Height (m)')
    leg = ax2.legend(handles=[pink,gray,red,blue],framealpha=1,loc=4)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GFSp/ec_h8fgen_'+dtfs+'.png')
    ax2.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFSp/ne_h8fgen_'+dtfs+'.png')

    #### Plot 3
    fig3 = plt.figure(figsize=(15,15))
    ax3 = fig3.add_subplot(111,projection=zH5_crs)

    ax3.coastlines(resolution='10m')
    ax3.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax3.add_feature(cfeature.STATES.with_scale('10m'))
    ax3.set_extent((260, 295, 25, 50))

    rhc = ax3.contourf(x,y,rh9,levels=[70,75,80,85,90,95,100],colors=['#73ea61','#5fe84a','#4be534','#38e21d','#32cb1a','#2cb517','#279e15'],alpha=0.6)
    cbar = fig3.colorbar(rhc,orientation = 'horizontal', aspect = 80, ax = ax3, pad = 0.01,
                        extendrect=False, ticks = [70,80,90,100])
    cbar.set_label('Relative Humidity (%)')
    ax3.barbs(x[wind_slice],y[wind_slice],u9[wind_slice,wind_slice],v9[wind_slice,wind_slice], length=7)
    fc = ax3.contour(x,y,h9_fgen,alpha=0.7,colors='fuchsia',levels=range(2,30,2),linewidths=3)
    ax3.contour(x,y,h9,colors='dimgray',levels=range(300,1200,30),linewidths=2)
    ax3.contour(x,y,t9c,colors='b',levels=range(-40,0,5),linewidths=1.5,linestyles='dashed')
    ax3.contour(x,y,t9c,colors='r',levels=range(0,40,50),linewidths=1.5,linestyles='dashed')
    ax3.set_title('925mb Forecast Summary',fontsize=16)
    ax3.set_title('GFS (v16) Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax3.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    pink = lines.Line2D([],[],linestyle='solid',color='fuchsia',label='Frontogenesis (K/100km/3hrs)')
    red = lines.Line2D([],[],linestyle='dashed',color='r',label='>=0C Temperatures (C)')
    blue = lines.Line2D([],[],linestyle='dashed',color='b',label='<0C Temperatures (C)')
    gray = lines.Line2D([],[],linestyle='solid',color='dimgray',label='Height (m)')
    leg = ax3.legend(handles=[pink,gray,red,blue],framealpha=1,loc=4)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GFSp/ec_h9fgen_'+dtfs+'.png')
    ax3.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFSp/ne_h9fgen_'+dtfs+'.png')

    fig4 = plt.figure(figsize=(15,15))
    ax4 = fig4.add_subplot(111,projection=zH5_crs)

    ax4.coastlines(resolution='10m')
    ax4.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax4.add_feature(cfeature.STATES.with_scale('10m'))
    ax4.set_extent((260, 295, 25, 50))

    ax4.contour(x,y,h7_fgen,alpha=0.7,colors='purple',levels=range(2,30,2),linewidths=3)
    ax4.contour(x,y,h8_fgen,alpha=0.7,colors='mediumorchid',levels=range(2,30,2),linewidths=3)
    ax4.contour(x,y,h9_fgen,alpha=0.7,colors='hotpink',levels=range(2,30,2),linewidths=3)
    purple = mpatches.Patch(color='purple',label='700mb')
    orchid = mpatches.Patch(color='mediumorchid',label='850mb')
    hpink = mpatches.Patch(color='hotpink',label='925mb')
    leg = ax4.legend(handles=[purple,orchid,hpink],title='Frontogenesis',loc=4,framealpha=1)
    ax4.set_title('Multi-Level Frontogenesis Comparison',fontsize=16)
    ax4.set_title('GFS (v16) Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax4.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    plt.savefig(output_dir+'/GFSp/ec_fgencomp_'+dtfs+'.png')
    ax4.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFSp/ne_fgencomp_'+dtfs+'.png')
    print(dtfs+' done')
