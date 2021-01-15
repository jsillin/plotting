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
    if int(hour) <5:
        init_hour = '00'
    elif int(hour) <10:
        init_hour = '06'
    elif int(hour) <16:
        init_hour = '12'
    elif int(hour) <21:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)
init_hour = get_init_hr(hour)
url='http://nomads.ncep.noaa.gov:80/dods/nam/nam'+mdate+'/nam_'+get_init_hr(hour)+'z'

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/NAM')

print(url)
ds = xr.open_dataset(url)
init_hr = dt.datetime(2020,11,27,12)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(12.21990800000,61.20556254545,0.110827275)
lons = np.arange((360-152.87862300000),(360.-49.47263081081),0.11338376)

#ds = ds.sel(lat = lats, lon = lons)
total_precip=ds['apcpsfc'].isel(time=0).squeeze()*.0393700787402

sub_w1 = 281
sub_w = 284
sub_e = 295
sub_n = 49
sub_s = 39

a = 10
b = 29

for i in range(a,b):
    fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    if i>a:
        prev_t2m = t2m
        prev_u10 = u_10m
        prev_v10 = v_10m
        prev_mslp = mslp

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)
    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tcdcclm':'tcc',
        'tmpprs': 'temperature',
        'ugrd10m': 'u',
        'vgrd10m': 'v',
        'hgtprs': 'height',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf'
    })
    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

    radius = 1
    kernel = np.zeros((2*radius+1,2*radius+1))
    y1,x1 = np.ogrid[-radius:radius+1,-radius:radius+1]
    mask=x1**2+y1**2 <=radius**2
    kernel[mask]=1

    snowc= gf(catsnow,np.max,footprint=kernel)
    rainc = gf(catrain,np.max,footprint=kernel)
    #zH5 = data['height'].metpy.sel(vertical=round(250.,0) * units.hPa)

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    zH5_crs = data['temperature'].metpy.cartopy_crs

    #data['temperature'].metpy.convert_units('degC')
    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.
    if i>a:
        three_hr_deltaT = t2m-prev_t2m
        print(three_hr_deltaT)

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)

    cloudcover = data['tcc'].squeeze()
    reflectivity = data['radar'].squeeze()
    hrly_precip = data['qpf'].squeeze()*0.0393700787402
    new_precip = hrly_precip-total_precip
    total_precip = hrly_precip

    rain = np.ma.masked_where(rainc==0,reflectivity)
    sleet = np.ma.masked_where(catsleet==0,reflectivity)
    ice = np.ma.masked_where(catice==0,reflectivity)
    snow = np.ma.masked_where(snowc==0,reflectivity)

    rain = ndimage.gaussian_filter(rain,sigma=1,order=0)

    mslp = data['mslp']/100.
    mslpc = mslp.squeeze()
    mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)
    ws2 = slice(5,-5,5)
    wind_slice = slice(5,-5,5)
    u_10m = data['u'].squeeze()
    v_10m = data['v'].squeeze()

    u_10m = u_10m*1.94384449
    v_10m = v_10m*1.94384449

    if i>a:
        avg_u10 = (u_10m+prev_u10)/2
        avg_v10 = (v_10m+prev_v10)/2
        avg_mslp = (mslp+prev_mslp)/2
        avgmslpc=ndimage.gaussian_filter(avg_mslp,sigma=1,order=0)

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))

    #fig.suptitle("NAM Forecast valid at " + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=36)

    ########## PLOTTING #######################################################
    tmp_2m = ax1.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5))
    cbr.set_label('2m Temperature (F)', fontsize = 14)

    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    q_levs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25]
    qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#ffcc00','#ff9500','#ff4800']
    qs_cols = ['#b8ffff','#82ffff','#00ffff','#00adad','#007575','#0087f5','#0039f5','#1d00f5','#7a00f5']
    q_levs_r = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3]
    try:
        ra = ax1.contourf(x,y,rain,cmap='Greens',levels=ref_levs,alpha=0.7)
    except:
        print('no rain')
    try:
        sn = ax1.contourf(x,y,snow,cmap='cool',levels=ref_levs,alpha=0.7)
    except:
        print('no snow')
    try:
        ip = ax1.contourf(x,y,sleet,cmap='autumn',levels=ref_levs,alpha=0.7)
    except:
        print('no sleet')
    try:
        zr = ax1.contourf(x,y,ice, cmap='RdPu',levels=ref_levs,alpha=0.7)
    except:
        print('no ice')
    ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)
    ax1.set_title('Precip Type, 2m Temperature (F), 10m Winds (kts), and MSLP (mb)',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n NAM Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax1.set_extent((sub_w1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/NAM/menh_3hrly_pytpe_v5_'+str(i)+'.png')
    if i>a:
        fig2 = plt.figure(figsize=(15,15))
        ax2 = fig2.add_subplot(111, projection=zH5_crs)
        deltc = ax2.contourf(x,y,three_hr_deltaT,cmap='RdBu_r',levels=range(-10,10,1),alpha=0.7)
        cbar = fig.colorbar(deltc, orientation='horizontal', aspect = 80, ax = ax2, pad = 0.01,
                            extendrect=False, ticks = range(-10,11,2))
        cbar.set_label('3hr Temperature Change (F)')
        ax2.barbs(x[wind_slice],y[wind_slice],avg_u10[wind_slice,wind_slice],avg_v10[wind_slice,wind_slice], length=7)
        ax2.contour(x, y, avgmslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
        ax2.coastlines(resolution='10m')
        ax2.add_feature(cfeature.BORDERS.with_scale('10m'))
        ax2.add_feature(cfeature.STATES.with_scale('10m'))
        ax2.set_title('3hr 2m Temp Change (F), 3hr Avg 10m Winds (kts), and MSLP (mb)',fontsize=14)
        ax2.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %HZ').item(),fontsize=11,loc='right')
        ax2.set_title('\n NAM Init: '+init_time.dt.strftime('%Y-%m-%d %HZ').item(),fontsize=11,loc='left')
        ax2.set_extent((sub_w1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
        plt.savefig(output_dir+'/NAM/NAM_3hr_delta_v5_'+str(i)+'.png')
    fcst_hr = str(0)
    print('Hour '+str(i)+' completed!')
    plt.close()
    timeelapsed = datetime.now()-startTime
    print(timeelapsed)







'''
url= 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs20200903/gfs_0p25_1hr_12z'
ds = xr.open_dataset(url)
t2m_ds = ds['tmp2m']
init_hr = t2m_ds['time'][0].values
#fc_hr = t2m.ds['time'][i].values
lats = np.arange(20,50,0.25)
lons = np.arange(240,300,0.25)
t2m = t2m_ds.sel(time = init_hr, lat = lats, lon = lons)
print(t2m)

fig = plt.figure(figsize = (12,12))
fig.clf()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_extent((240,300, 20, 50), crs = ccrs.PlateCarree())
t2m_c = ax.contourf(t2m, cmap='RdPu')
plt.savefig('testingnomads6.png')
'''
