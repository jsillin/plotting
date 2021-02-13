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
import matplotlib.lines as lines
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

total_precip=ds['apcpsfc'].isel(time=0).squeeze()*.0393700787402

for i in range(1,29):
    fc_hr = init_hr+dt.timedelta(hours=1*i)
    forecast_hour = times[0].values

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
        'ugrdprs': 'u',
        'vgrdprs': 'v',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'refcclm':'radar',
        'apcpsfc':'qpf',
        'hgtprs':'height',
        'gustsfc':'sfcgust',
        'pwatclm':'pwat',
        'dpt2m':'td',
        'absvprs':'avort',
        'rhprs':'rh',
        'ugrd10m':'u10',
        'vgrd10m':'v10'
    })

    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

    zH5 = data['temperature'].squeeze()
    zH5_crs = zH5.metpy.cartopy_crs

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    dx, dy = mpcalc.lat_lon_grid_deltas(ds.lon, ds.lat)

    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    h2 = data['height'].sel(lev=round(250.0,0)).squeeze()
    hgt5 = data['height'].sel(lev=round(500.0,0)).squeeze()
    hgt7 = data['height'].sel(lev=round(700.0,0)).squeeze()
    hgt8 = data['height'].sel(lev=round(850.0,0)).squeeze()

    h7=hgt7
    h8=hgt8

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
    u_10m = data['u10'].squeeze()*1.94384449
    v_10m = data['v10'].squeeze()*1.94384449

    wspd2 = ((u2**2)+(v2**2))**.5
    wspd5 = ((u5**2)+(v5**2))**.5
    wspd7 = ((u7**2)+(v7**2))**.5
    wspd8 = ((u8**2)+(v8**2))**.5
    wspd9 = ((u9**2)+(v9**2))**.5

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

    rh7 = data['rh'].sel(lev=round(700.0,0)).squeeze()

    av5 = data['avort'].sel(lev=round(500.0,0)).squeeze()*1e5
    #h7_fgen = mpcalc.frontogenesis(t7u,u7u,v7u,dx,dy,dim_order='xy')
    #h7_fgen = h7_fgen*1000*100*3600*3 ##convert to units of K/100km/3hrs
    #print(h7_fgen)

    pwat = data['pwat'].squeeze()*0.0393700787402
    cloudcover = data['tcc'].squeeze()
    reflectivity = data['radar'].squeeze()
    hrly_precip = data['qpf'].squeeze()*0.0393700787402
    total_precip = total_precip+hrly_precip

    rain = np.ma.masked_where(catrain==0,reflectivity)
    sleet = np.ma.masked_where(catsleet==0,reflectivity)
    ice = np.ma.masked_where(catice==0,reflectivity)
    snow = np.ma.masked_where(catsnow==0,reflectivity)

    mslpc = data['mslp'].squeeze()/100
    mslpc=ndimage.gaussian_filter(mslpc,sigma=3,order=0)
    wind_slice = slice(10,-10,10)
    wind_slice_s = slice(20,-20,20)
    #u_10m = data['u'].squeeze()*1.94384449
    #v_10m = data['v'].squeeze()*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    swg = data['sfcgust'].squeeze()

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(43,15))

    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
    gs.update(wspace=0.01, hspace=0.01)
    ax1 = fig.add_subplot(gs[:, 1], projection = zH5_crs)
    ax2 = fig.add_subplot(gs[0, 0], projection = zH5_crs)
    ax3 = fig.add_subplot(gs[1, 0], projection = zH5_crs)
    ax4 = fig.add_subplot(gs[0, 2], projection = zH5_crs)
    ax5 = fig.add_subplot(gs[1, 2], projection = zH5_crs)

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

    #fig.suptitle("HRRR Forecast valid at " + time[0].dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=36)
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    print(dtfs)
    ########## PLOTTING #######################################################
    tmp_2m = ax1.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32],linewidths=2)
    h9_0c = ax1.contour(x,y,ndimage.gaussian_filter(t9,sigma=5,order=0),colors='royalblue',alpha=0.8,levels=[0],linewidths=3)
    h8_0c = ax1.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='cornflowerblue',alpha=0.8,levels=[0],linewidths=3)
    h7_0c = ax1.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='lightsteelblue',alpha=0.8,levels=[0],linewidths=3)
    h8_wsc = ax1.contour(x,y,gaussian_filter(wspd8,sigma=5,order=0),cmap='YlOrRd',alpha=0.7,levels=range(40,90,10),linewidths=3)
    cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01, shrink = .7,
                        extendrect=False, ticks = range(-20,100,5))
    cbr.set_label('2m Temperature (F)', fontsize = 14)

    m_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=3)
    m_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=1, fmt='%i mb', rightside_up=True, use_clabeltext=True)
    #h5_contour = ax1.contour(x,y,hgt5,colors='gray',levels=range(4800,6000,60),linewidths=3)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]

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

    ax1.barbs(x[wind_slice],y[wind_slice],u8[wind_slice,wind_slice],v8[wind_slice,wind_slice], length=7)
    ax1.set_title('Precip Type, 2m Temperature (F), 850mb Winds (kts), and MSLP (mb)',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n NAM Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11)

    ### AX2 = top left = 200mb
    w2c = ax2.contourf(x,y,wspd2, alpha = 0.7,levels=range(50,200,5), cmap='Blues')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(w2c,orientation='vertical',pad=0.01,ax=ax2,shrink=.7,aspect=50,extendrect=False, ticks=np.arange(50,200,20))
    cbar2.set_label('250mb Wind Speed (kts)',fontsize=14)
    ax2.barbs(x[wind_slice_s],y[wind_slice_s],u2[wind_slice_s,wind_slice_s],v2[wind_slice_s,wind_slice_s],length=7)
    h2c = ax2.contour(x,y,h2,colors='dimgray', levels = range(9000,11000,60),linewidths=1.5)

    rhlevs = []
    for j in range(1,35):
        lev = 0.69+0.01*j
        rhlevs.append(lev)

    ### AX3 = bottom left = 700mb
    rhc = ax3.contourf(x,y,rh7, alpha = 0.7, cmap = 'winter_r',transform=zH5_crs, levels=range(68,102,2))
    h7c = ax3.contour(x,y,hgt7,colors='k',levels=range(2700,3300,30),linewidths=2)
    t7c = ax3.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='r',levels=range(0,18,3),linestyles='dashed',linewidths=2)
    t7c2 = ax3.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    cbar3 = fig.colorbar(rhc,orientation='vertical',pad=0.01,shrink=.8,ax=ax3,aspect=50,extendrect=False,ticks=[70,80,90,100])
    cbar3.set_label('700mb RH',fontsize=14)
    ax3.contour(x,y,h7_fgen,alpha=0.7,colors='fuchsia',levels=range(4,50,2),linewidths=3)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax3.legend(handles=[purple,dashed_blue_line,dashed_red_line,black_line],loc=4,framealpha=1)


    ### AX4 = top right = 500mb
    a5c = ax4.contourf(x,y,av5,cmap='autumn_r',levels=range(10,65,2),alpha=0.8)
    a5cb = fig.colorbar(a5c, orientation = 'vertical', aspect = 80, ax = ax4, pad = 0.01, shrink=.8,
                        extendrect=False, ticks = range(10,60,5))
    a5cb.set_label('500mb Absolute Vorticity ($s^{-1}$)', fontsize = 12)
    h5c = ax4.contour(x,y,hgt5,colors='k',levels=range(4800,6000,60),linewidths=2)
    #cb = fig.colorbar(capep, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000])
    #cb.set_label('CAPE (J/kg)', size='large')
    ax4.barbs(x[wind_slice_s],y[wind_slice_s],u5[wind_slice_s,wind_slice_s],v5[wind_slice_s,wind_slice_s],length=7)

    #ax4.barbs(x[ws2],y[ws2],s5u[ws2,ws2],s5v[ws2,ws2], length = 7)

    ### AX5 = bottom right = 850mb
    pwlevs = []
    for j in range(1,25):
        lev = 0.1*j
        pwlevs.append(lev)

    pwatc = ax5.contourf(x,y,pwat, levels=pwlevs, cmap = 'BrBG',alpha=0.5)
    h8c = ax5.contour(x,y,ndimage.gaussian_filter(hgt8,sigma=5,order=0),colors='k',levels=range(1020,1800,30),linewidths=2)
    t8c = ax5.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='r',levels=range(0,30,3),linestyles='dashed',linewidth=2)
    t8c2 = ax5.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    cbar4 = fig.colorbar(pwatc, orientation='vertical',shrink=0.8,pad=0.01,ax=ax5, aspect = 50, extendrect=False, ticks = [0.5,1,1.5,2,2.5])
    cbar4.set_label('Precipitable Water (in)')
    ax5.contour(x,y,h8_fgen,alpha=0.7,colors='fuchsia',levels=range(4,50,2),linewidths=3)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax5.legend(handles=[purple,dashed_blue_line,dashed_red_line,black_line],loc=4,framealpha=1)


    #refp3 = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs)
    #cbar4 = fig.colorbar(tccp,orientation='horizontal',pad=0.01,ax=ax4,aspect=50,extendrect=False)

    #ax1.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax2.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax3.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax4.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax5.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #north=50, south=15, east=-70, west=-115

    sub_w1 = 269
    sub_w = 270
    sub_e = 294
    sub_n = 50
    sub_s = 32

    ax1.set_extent((sub_w1-2, sub_e+2, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/NAM/ec_fivepanelwinter_v9_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.01)
    plt.clf()
    plt.close()

    ##########################################################
    ws2r = slice(15,-15,15)

    ### SECOND PLOT ###
    fig2 = plt.figure(figsize=(15,15))
    ax6 = fig2.add_subplot(111, projection = zH5_crs)

    ax6.coastlines(resolution='10m')
    ax6.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax6.add_feature(cfeature.STATES.with_scale('10m'))

    tmp_2m = ax6.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax6.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr6 = fig2.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax6, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5))
    cbr6.set_label('2m Temperature (F)', fontsize = 14)

    h_contour = ax6.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]

    try:
        ra = ax6.contourf(x,y,rain,cmap='Greens',levels=ref_levs,alpha=0.7)
    except:
        print('no rain')
    try:
        sn = ax6.contourf(x,y,snow,cmap='cool',levels=ref_levs,alpha=0.7)
    except:
        print('no snow')
    try:
        ip = ax6.contourf(x,y,sleet,cmap='autumn',levels=ref_levs,alpha=0.7)
    except:
        print('no sleet')
    try:
        zr = ax6.contourf(x,y,ice, cmap='RdPu',levels=ref_levs,alpha=0.7)
    except:
        print('no ice')

    ax6.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)

    ax6.set_extent((sub_w-1, sub_e+1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax6.set_title('NAM Composite Forecast valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/NAM/NE_ptype_composite_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.01)
    plt.clf()
    plt.close()
    ### END SECOND PLOT ###

    ### THIRD PLOT ###
    fig3 = plt.figure(figsize=(15,15))
    ax7 = fig3.add_subplot(111,projection=zH5_crs)

    ax7.coastlines(resolution='10m')
    ax7.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax7.add_feature(cfeature.STATES.with_scale('10m'))

    w2c = ax7.contourf(x,y,wspd2, alpha = 0.7,levels=range(50,200,5), cmap='Blues')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(w2c,orientation='horizontal',pad=0.01,ax=ax7,shrink=.8,aspect=50,extendrect=False, ticks=np.arange(50,200,20))
    cbar2.set_label('250mb Wind Speed (kts)',fontsize=14)
    ax7.barbs(x[ws2r],y[ws2r],u2[ws2r,ws2r],v2[ws2r,ws2r],length=7)
    h2c = ax7.contour(x,y,h2,colors='dimgray', levels = range(9000,11000,60),linewidths=1.5)

    ax7.set_extent((sub_w-1, sub_e+1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax7.set_title('NAM 250mb Forecast Valid ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/NAM/NE_250_2_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.01)
    plt.clf()
    plt.close()

    ### FOURTH PLOT ###
    fig4 = plt.figure(figsize=(15,15))
    ax8 = fig4.add_subplot(111,projection=zH5_crs)

    ax8.coastlines(resolution='10m')
    ax8.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax8.add_feature(cfeature.STATES.with_scale('10m'))

    rhc = ax8.contourf(x,y,rh7, alpha = 0.7, cmap = 'winter_r',transform=zH5_crs, levels=range(68,102,2))
    h7c = ax8.contour(x,y,hgt7,colors='k',levels=range(2700,3300,30),linewidths=2)
    t7c = ax8.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='r',levels=range(0,18,3),linestyles='dashed',linewidths=2)
    t7c2 = ax8.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    cbar3 = fig.colorbar(rhc,orientation='horizontal',pad=0.01,shrink=.8,ax=ax8,aspect=50,extendrect=False,ticks=[70,80,90,100])
    cbar3.set_label('700mb RH (%)',fontsize=14)

    ax8.set_title('NAM 700mb Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    ax8.set_extent((sub_w-1, sub_e+1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot

    h_contouqr = ax8.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contouqr.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ax8.contour(x,y,h7_fgen,alpha=0.7,colors='fuchsia',levels=range(2,30,2),linewidths=3)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax8.legend(handles=[purple,dashed_blue_line,dashed_red_line,black_line],loc=4,framealpha=1)


    plt.savefig(output_dir+'/NAM/NE_700_2_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.01)
    plt.clf()
    plt.close()

    ### FIFTH PLOT ###
    fig5 = plt.figure(figsize=(15,15))
    ax9 = fig5.add_subplot(111,projection=zH5_crs)

    ax9.coastlines(resolution='10m')
    ax9.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax9.add_feature(cfeature.STATES.with_scale('10m'))

    a5c = ax9.contourf(x,y,av5,cmap='autumn_r',levels=range(10,65,2),alpha=0.8)
    a5cb = fig.colorbar(a5c, orientation = 'horizontal', aspect = 80, ax = ax4, pad = 0.01,
                        extendrect=False, ticks = range(10,60,5))
    a5cb.set_label('500mb Absolute Vorticity ($s^{-1}$)', fontsize = 12)
    h5c = ax9.contour(x,y,hgt5,colors='k',levels=range(4800,6000,60),linewidths=2)
    #cb = fig.colorbar(capep, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000])
    #cb.set_label('CAPE (J/kg)', size='large')
    ax9.barbs(x[ws2r],y[ws2r],u5[ws2r,ws2r],v5[ws2r,ws2r],length=7)

    ax9.set_title('NAM 500mb Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    ax9.set_extent((sub_w-1, sub_e+1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/NAM/NE_500_2_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.01)
    plt.clf()
    plt.close()

    #### SIXTH PLOT ####
    fig6 = plt.figure(figsize=(15,15))
    ax10 = fig6.add_subplot(111,projection=zH5_crs)

    ax10.coastlines(resolution='10m')
    ax10.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax10.add_feature(cfeature.STATES.with_scale('10m'))

    pwlevs = []
    for j in range(1,25):
        lev = 0.1*j
        pwlevs.append(lev)

    pwatc = ax10.contourf(x,y,pwat, levels=pwlevs, cmap = 'BrBG',alpha=0.5)
    h8c = ax10.contour(x,y,ndimage.gaussian_filter(hgt8,sigma=5,order=0),colors='k',levels=range(1020,1800,30),linewidths=2)
    t8c = ax10.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='r',levels=range(0,30,3),linestyles='dashed',linewidth=2)
    t8c2 = ax10.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    ax10.contour(x,y,h8_fgen,alpha=0.7,colors='fuchsia',levels=range(2,30,2),linewidths=3)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax10.legend(handles=[purple,dashed_blue_line,dashed_red_line,black_line],loc=4,framealpha=1)
    cbar4 = fig.colorbar(pwatc, orientation='horizontal',pad=0.01,ax=ax10, aspect = 50, extendrect=False, ticks = [0.5,1,1.5,2,2.5])
    cbar4.set_label('Precipitable Water (in)')
    ax10.barbs(x[ws2r],y[ws2r],u8[ws2r,ws2r],v8[ws2r,ws2r], length=7)
    ax10.set_title('NAM 850mb Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    ax10.set_extent((sub_w-1, sub_e+1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/NAM/NE_850_2_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.01)
    plt.clf()
    plt.close()
    print(str(i)+'_Done!')
