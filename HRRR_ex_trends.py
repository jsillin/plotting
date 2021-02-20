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
import matplotlib.lines as lines

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
    if int(hour) <3:
        init_hour = '00'
    elif int(hour) <9:
        init_hour = '06'
    elif int(hour) <15:
        init_hour = '12'
    elif int(hour) <21:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

def get_init_time(current_datetime):
    if current_datetime.hour<3:
        old_time = current_datetime-dt.timedelta(hours=3)
        hour = '00'
        month = current_datetime.month
        day = current_datetime.day
        year = current_datetime.year
        old_hour = '18'
        old_month = old_time.month
        old_day = old_time.day
        old_year = old_time.year
    elif current_datetime.hour<9:
        hour = '06'
        month = current_datetime.month
        day = current_datetime.day
        year = current_datetime.year
        old_month = month
        old_day = day
        old_year = year
        old_hour = '00'
    elif current_datetime.hour<15:
        hour = '12'
        month = current_datetime.month
        day = current_datetime.day
        year = current_datetime.year
        old_month = month
        old_day = day
        old_year = year
        old_hour = '06'
    elif current_datetime.hour<21:
        hour = '18'
        month = current_datetime.month
        day = current_datetime.day
        year = current_datetime.year
        old_month = month
        old_day = day
        old_year = year
        old_hour = '12'
    else:
        start_time = current_datetime+dt.timedelta(hours=3)
        month = start_time.month
        day = start_time.day
        year = start_time.year
        hour = '00'
        old_month = current_datetime.month
        old_day = current_datetime.day
        old_year = current_datetime.year
        old_hour = '18'

    if month<10:
        month = '0'+str(month)
    else:
        month = str(month)

    if day <10:
        day = '0'+str(day)
    else:
        day = str(day)

    year = str(year)

    if old_month<10:
        old_month = '0'+str(old_month)
    else:
        old_month = str(old_month)

    if old_day <10:
        old_day = '0'+str(old_day)
    else:
        old_day = str(old_day)

    old_year = str(old_year)

    mdate = year+month+day
    odate = old_year+old_month+old_day
    output=[mdate,hour,odate,old_hour]
    return output

times = get_init_time(startTime)
init_hour = times[1]
mdate = times[0]
odate = times[2]
ohour = times[3]
url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr'+mdate+'/hrrr_sfc.t'+init_hour+'z'
old_url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr'+odate+'/hrrr_sfc.t'+ohour+'z'
print(url)
print(old_url)

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/HRRR_ex')
#Parse data using MetPy
ds = xr.open_dataset(url)
ods= xr.open_dataset(old_url)
init_hr = dt.datetime(int(year),int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time
otimes = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(25,55,0.25)
lons = np.arange(260,310,0.25)

total_precip=ds['apcpsfc'].isel(time=0).squeeze()*.0393700787402
old_precip = ods['apcpsfc'].isel(time=6).squeeze()*.0393700787402
for i in range(1,43):
    fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    old_data = ods.metpy.parse_cf()
    old_data = old_data.isel(time=i+6)

    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tcdcclm':'tcc',
        'tmpprs': 'temperature',
        'ugrdprs': 'u',
        'vgrdprs': 'v',
        'mslmamsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf',
        'hgt500mb':'h5',
        'hgt700mb':'h7',
        'hgt850mb':'h8',
        'gustsfc':'sfcgust',
        'pwatclm':'pwat',
        'dptprs':'td',
        'sbt123toa':'simsat'
    })

    old_data = old_data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tcdcclm':'tcc',
        'tmpprs': 'temperature',
        'ugrdprs': 'u',
        'vgrdprs': 'v',
        'mslmamsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf',
        'hgt500mb':'h5',
        'hgt700mb':'h7',
        'hgt850mb':'h8',
        'gustsfc':'sfcgust',
        'pwatclm':'pwat',
        'dptprs':'td',
        'sbt123toa':'simsat'
    })

    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

    zH5 = data['temperature'].squeeze()
    zH5_crs = zH5.metpy.cartopy_crs

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    otime = old_data['temperature'].metpy.time

    print(time)
    print(otime)
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    dx, dy = mpcalc.lat_lon_grid_deltas(ds.lon, ds.lat)

    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    old_t2m = old_data['sfc_temp'].squeeze()
    old_t2m = ((old_t2m - 273.15)*(9./5.))+32.

    temp_diff = t2m-old_t2m

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)

    hgt5 = h5 = data['h5'].squeeze()
    hgt7 = h7 = data['h7'].squeeze()
    hgt8 = h8 = data['h8'].squeeze()

    ohgt5 = oh5 = old_data['h5'].squeeze()
    ohgt7 = oh7 = old_data['h7'].squeeze()
    ohgt8 = oh8 = old_data['h8'].squeeze()

    u2 = data['u'].sel(lev=round(250.0,0)).squeeze()*1.94384449
    v2 = data['v'].sel(lev=round(250.0,0)).squeeze()*1.94384449
    u5 = data['u'].sel(lev=round(500.0,0)).squeeze()*1.94384449
    v5 = data['v'].sel(lev=round(500.0,0)).squeeze()*1.94384449
    u7 = data['u'].sel(lev=round(700.0,0)).squeeze()*1.94384449
    v7 = data['v'].sel(lev=round(700.0,0)).squeeze()*1.94384449
    u8 = data['u'].sel(lev=round(850.0,0)).squeeze()*1.94384449
    v8 = data['v'].sel(lev=round(850.0,0)).squeeze()*1.94384449
    u9 = data['u'].sel(lev=round(925.0,0)).squeeze()*1.94384449
    v9 = data['v'].sel(lev=round(925.0,0)).squeeze()*1.94384449

    ou2 = old_data['u'].sel(lev=round(250.0,0)).squeeze()*1.94384449
    ov2 = old_data['v'].sel(lev=round(250.0,0)).squeeze()*1.94384449
    ou5 = old_data['u'].sel(lev=round(500.0,0)).squeeze()*1.94384449
    ov5 = old_data['v'].sel(lev=round(500.0,0)).squeeze()*1.94384449
    ou7 = old_data['u'].sel(lev=round(700.0,0)).squeeze()*1.94384449
    ov7 = old_data['v'].sel(lev=round(700.0,0)).squeeze()*1.94384449
    ou8 = old_data['u'].sel(lev=round(850.0,0)).squeeze()*1.94384449
    ov8 = old_data['v'].sel(lev=round(850.0,0)).squeeze()*1.94384449
    ou9 = old_data['u'].sel(lev=round(925.0,0)).squeeze()*1.94384449
    ov9 = old_data['v'].sel(lev=round(925.0,0)).squeeze()*1.94384449


    wspd2 = ((u2**2)+(v2**2))**.5
    wspd5 = ((u5**2)+(v5**2))**.5
    wspd7 = ((u7**2)+(v7**2))**.5
    wspd8 = ((u8**2)+(v8**2))**.5
    wspd9 = ((u9**2)+(v9**2))**.5

    owspd2 = ((ou2**2)+(ov2**2))**.5
    owspd5 = ((ou5**2)+(ov5**2))**.5
    owspd7 = ((ou7**2)+(ov7**2))**.5
    owspd8 = ((ou8**2)+(ov8**2))**.5
    owspd9 = ((ou9**2)+(ov9**2))**.5

    h8_wind_diff = wspd8-owspd8
    h5_hgt_diff = h5-oh5

    td7 = data['td'].sel(lev=round(700.0,0)).squeeze()

    t5 = data['temperature'].sel(lev=round(500.0,0)).squeeze()-273.15
    t7 = data['temperature'].sel(lev=round(700.0,0)).squeeze()-273.15
    t8 = data['temperature'].sel(lev=round(850.0,0)).squeeze()-273.15
    t9 = data['temperature'].sel(lev=round(925.0,0)).squeeze()-273.15

    t7u = t7.values * units.degC
    u7u = u7.values * units.knots
    v7u = v7.values * units.knots

    t8u = t8.values * units.degC
    u8u = u8.values * units.knots
    v8u = v8.values * units.knots

    td7u = td7.values * units.degK
    rh7 = mpcalc.relative_humidity_from_dewpoint(t7u,td7u)
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
    h7_fgen = ndimage.gaussian_filter(h7_fgen,sigma=2,order=0)

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
    h8_fgen = ndimage.gaussian_filter(h8_fgen,sigma=2,order=0)

    #925
    u9k = u9*1.94384449
    v9k = v9*1.94384449
    t9c = t9-273.15
    t9c = ndimage.gaussian_filter(t9c,sigma=2,order=0)

    t9 = t9*units.K
    u9 = u9*units.meters/units.seconds
    v9 = v9*units.meters/units.seconds

    h9_fgen = mpcalc.frontogenesis(t9.data,u9.data,v9.data,dx,dy)
    h9_fgen = h9_fgen*1000*100*3600*3 ##convert to units of K/100km/3hrs
    h9_fgen = ndimage.gaussian_filter(h9_fgen,sigma=1,order=0)

    simsat = data['simsat'].squeeze()
    pwat = data['pwat'].squeeze()*0.0393700787402
    cloudcover = data['tcc'].squeeze()
    reflectivity = data['radar'].squeeze()
    hrly_precip = data['qpf'].squeeze()*0.0393700787402
    old_hrly_precip = old_data['qpf'].squeeze()*0.0393700787402
    total_precip = total_precip+hrly_precip
    old_precip = old_precip+old_hrly_precip

    precip_diff = total_precip-old_precip

    rain = np.ma.masked_where(catrain==0,reflectivity)
    sleet = np.ma.masked_where(catsleet==0,reflectivity)
    ice = np.ma.masked_where(catice==0,reflectivity)
    snow = np.ma.masked_where(catsnow==0,reflectivity)

    mslpc = data['mslp'].squeeze()/100
    omslp = old_data['mslp'].squeeze()/100
    mslpc=ndimage.gaussian_filter(mslpc,sigma=3,order=0)
    omslp=ndimage.gaussian_filter(omslp,sigma=3,order=0)

    mslp_diff = mslpc-omslp

    wind_slice = slice(25,-25,25)
    wind_slice_s = slice(40,-40,40)
    wind_slice = slice(36,-36,36)
    wind_slice_ne = slice(18,-18,18)
    wind_slice_me = slice(9,-9,9)
    u_10m = data['u'].squeeze()*1.94384449
    v_10m = data['v'].squeeze()*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    swg = data['sfcgust'].squeeze()
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    ##############################################################

    fig = plt.figure(figsize=(42,15))

    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
    ax1 = fig.add_subplot(gs[:, 1], projection = zH5_crs)
    ax2 = fig.add_subplot(gs[0, 0], projection = zH5_crs)
    ax3 = fig.add_subplot(gs[1, 0], projection = zH5_crs)
    ax5 = fig.add_subplot(gs[0, 2], projection = zH5_crs)
    ax4 = fig.add_subplot(gs[1, 2], projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))

    ax2.coastlines(resolution='10m')
    ax2.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax2.add_feature(cfeature.STATES.with_scale('10m'))

    ax3.coastlines(resolution='10m')
    ax3.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax3.add_feature(cfeature.STATES.with_scale('10m'))

    ax4.coastlines(resolution='10m')
    ax4.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax4.add_feature(cfeature.STATES.with_scale('10m'))

    ax5.coastlines(resolution='10m')
    ax5.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax5.add_feature(cfeature.STATES.with_scale('10m'))


    tdiff = ax1.contourf(x,y,temp_diff,cmap='RdBu_r',levels=range(-20,21,1),transform=zH5_crs,extend='both')
    old_32 = ax1.contour(x,y,old_t2m,colors='steelblue',alpha=0.8,levels=[32],linewidths=1.5,transform=zH5_crs)
    new_32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32],linewidths=2,transform=zH5_crs)
    cbar = fig.colorbar(tdiff, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5), shrink=0.7)
    cbar.set_label('2m Temperature Difference (New Forecast-Old Forecast)', fontsize = 14)

    blue = mpatches.Patch(color='b', label='New Run 32F')
    lblue = mpatches.Patch(color='steelblue', label='Old Run 32F')
    leg = ax1.legend(handles=[blue,lblue],loc=4,framealpha=1)
    leg.set_zorder(100)

    presdiff = ax2.contourf(x,y,mslp_diff,cmap='RdYlBu_r',levels=range(-20,20,1),transform=zH5_crs,extend='both')
    cbar2 = fig.colorbar(presdiff,orientation='vertical',pad=0.01,ax=ax2,aspect=50,extendrect=False,ticks=range(-20,25,5))
    cbar2.set_label('MSLP Difference',fontsize=14)
    msl = ax2.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)

    pdiff = ax3.contourf(x,y,precip_diff,cmap='BrBG',levels=np.linspace(-2.5,2.5,60),transform=zH5_crs,extend='both')
    cbar3 = fig.colorbar(pdiff,orientation='vertical',pad=0.01,ax=ax3,aspect=50,extendrect=False,ticks=[-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5])
    cbar3.set_label('Precipitation Difference',fontsize=14)

    h8wind = ax5.contourf(x,y,h8_wind_diff,cmap='seismic',levels=range(-50,50,2),transform=zH5_crs,extend='both')
    cbar4 = fig.colorbar(h8wind, orientation='vertical', pad = 0.01, aspect = 50, ax = ax5, extendrect=False, ticks=range(-50,60,10))
    h8c = ax5.contour(x,y,ndimage.gaussian_filter(hgt8,sigma=5,order=0),colors='k',levels=range(1020,1800,30),linewidths=2)
    cbar4.set_label('850mb Wind Difference',fontsize=14)

    h5hgt = ax4.contourf(x,y,h5_hgt_diff,cmap='PuOr',levels=range(-150,150,10),transform=zH5_crs,extend='both')
    cbar5 = fig.colorbar(h5hgt, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=range(-150,180,30))
    cbar5.set_label('500mb Difference')
    h5c = ax4.contour(x,y,ndimage.gaussian_filter(hgt5,sigma=5,order=0),colors='k',levels=range(4800,6300,60),linewidths=2)

    sub_w1 = 260
    sub_w = 262
    sub_e = 295
    sub_n = 50
    sub_s = 25

    ax1.set_title('HRRR 6-hour Forecast Trends',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n HRRR Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    ax1.set_extent((sub_w1-1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    #fig.canvas.draw()
    fig.tight_layout()
    plt.savefig(output_dir+'/HRRR_ex/EC_fivepanel_dprogdt6_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    ax1.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/HRRR_ex/NE_fivepanel_dprogdt6_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.clf()
