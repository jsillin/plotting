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
import supplementary_tools as spt

mdate = spt.get_init_time('HRRR')[0]
init_hour = spt.get_init_time('HRRR')[1]
url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr'+mdate+'/hrrr_sfc.t'+init_hour+'z'
#url='http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr20201231/hrrr_sfc.t00z'
print(url)

# Create new directory
output_dir = mdate+'_'+init_hour+'00'
#output_dir = '20201231_0000'
spt.mkdir_p(output_dir)
spt.mkdir_p(output_dir+'/HRRR_ex')
#Parse data using MetPy
ds = xr.open_dataset(url)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(25,55,0.25)
lons = np.arange(260,310,0.25)

total_precip=ds['apcpsfc'].isel(time=0).squeeze()*.0393700787402

for i in range(1,49):

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
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    dx, dy = mpcalc.lat_lon_grid_deltas(ds.lon, ds.lat)

    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)

    hgt5 = h5 = data['h5'].squeeze()
    hgt7 = h7 = data['h7'].squeeze()
    hgt8 = h8 = data['h8'].squeeze()

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

    wspd2 = ((u2**2)+(v2**2))**.5
    wspd5 = ((u5**2)+(v5**2))**.5
    wspd7 = ((u7**2)+(v7**2))**.5
    wspd8 = ((u8**2)+(v8**2))**.5
    wspd9 = ((u9**2)+(v9**2))**.5

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
    total_precip = total_precip+hrly_precip

    rain = np.ma.masked_where(catrain==0,reflectivity)
    sleet = np.ma.masked_where(catsleet==0,reflectivity)
    ice = np.ma.masked_where(catice==0,reflectivity)
    snow = np.ma.masked_where(catsnow==0,reflectivity)

    mslpc = data['mslp'].squeeze()/100
    mslpc=ndimage.gaussian_filter(mslpc,sigma=3,order=0)
    wind_slice = slice(25,-25,25)
    wind_slice_s = slice(40,-40,40)
    wind_slice = slice(36,-36,36)
    wind_slice_ne = slice(18,-18,18)
    wind_slice_me = slice(9,-9,9)
    u_10m = data['u'].squeeze()*1.94384449
    v_10m = data['v'].squeeze()*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    swg = data['sfcgust'].squeeze()

    ########## SET UP FIGURE ##################################################
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
    cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5))
    cbr.set_label('2m Temperature (F)', fontsize = 14)

    m_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=3)
    m_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

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
    ax1.set_title('Precip Type, 2m Temperature (F), 10m Winds (kts), and MSLP (mb)',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    ### AX2 = top left = 200mb
    w2c = ax2.contourf(x,y,wspd2, alpha = 0.7,levels=range(50,200,5), cmap='Blues')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(w2c,orientation='vertical',pad=0.01,ax=ax2,shrink=.8,aspect=50,extendrect=False, ticks=np.arange(50,200,20))
    cbar2.set_label('250mb Wind Speed (kts)',fontsize=14)
    ax2.barbs(x[wind_slice_s],y[wind_slice_s],u2[wind_slice_s,wind_slice_s],v2[wind_slice_s,wind_slice_s],length=7)

    rhlevs = []
    for j in range(1,30):
        lev = 0.69+0.01*j
        rhlevs.append(lev)

    ### AX3 = bottom left = 700mb
    fc = ax3.contour(x,y,h7_fgen,alpha=0.7,colors='fuchsia',levels=range(4,60,4),linewidths=3)

    rhc = ax3.contourf(x,y,rh7, alpha = 0.7, cmap = 'winter_r',transform=zH5_crs, levels=rhlevs,extend='max')
    h7c = ax3.contour(x,y,hgt7,colors='k',levels=range(2700,3300,30),linewidths=2)
    t7c = ax3.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='r',levels=range(0,18,3),linestyles='dashed',linewidths=2)
    t80c = ax3.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='purple',levels=[0],linestyles='dashed',linewidths=2)
    t7c2 = ax3.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    cbar3 = fig.colorbar(rhc,orientation='vertical',pad=0.01,shrink=.8,ax=ax3,aspect=50,extendrect=False,ticks=[0.7,0.8,0.9,1])
    cbar3.set_label('700mb RH',fontsize=14)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_purple_line = lines.Line2D([],[],linestyle='dashed',color='purple',label='Temperature (0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax3.legend(handles=[purple,dashed_blue_line,dashed_purple_line,dashed_red_line,black_line],loc=4,framealpha=1)

    ### AX4 = top right = 500mb
    #refp = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs) #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])
    ss = ax4.contourf(x,y,simsat,cmap='bone',levels=range(170,330,5),alpha=0.7)
    h9c = ax4.contour(x,y,hgt5,colors='k',levels=range(300,1200,30),linewidths=2)
    f9c = ax4.contour(x,y,h9_fgen,alpha=0.7,colors='fuchsia',levels=range(6,60,4),linewidths=3)
    t9c = ax4.contour(x,y,ndimage.gaussian_filter(t9,sigma=5,order=0),colors='r',levels=range(0,18,3),linestyles='dashed',linewidths=2)
    t90c = ax4.contour(x,y,ndimage.gaussian_filter(t9,sigma=5,order=0),colors='purple',levels=[0],linestyles='dashed',linewidths=2)
    t9c2 = ax4.contour(x,y,ndimage.gaussian_filter(t9,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    cb = fig.colorbar(ss, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000])
    cb.set_label('Simulated Satellite')
    ax4.barbs(x[wind_slice_s],y[wind_slice_s],u9[wind_slice_s,wind_slice_s],v9[wind_slice_s,wind_slice_s],length=7)

    #ax4.barbs(x[ws2],y[ws2],s5u[ws2,ws2],s5v[ws2,ws2], length = 7)

    ### AX5 = bottom right = 850mb
    pwlevs = []
    for j in range(1,25):
        lev = 0.1*j
        pwlevs.append(lev)

    pwatc = ax5.contourf(x,y,pwat, levels=pwlevs, cmap = 'BrBG',alpha=0.5)
    h8c = ax5.contour(x,y,ndimage.gaussian_filter(hgt8,sigma=5,order=0),colors='k',levels=range(1020,1800,30),linewidths=2)
    t8c = ax5.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='r',levels=range(3,30,3),linestyles='dashed',linewidth=2)
    t80c = ax5.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='purple',levels=[0],linestyles='dashed',linewidths=2)
    t8c2 = ax5.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidth=2)
    cbar4 = fig.colorbar(pwatc, orientation='vertical',pad=0.01,ax=ax5, aspect = 50, shrink=0.8, extendrect=False, ticks = [0.5,1,1.5,2,2.5])
    cbar4.set_label('Precipitable Water (in)')
    ax5.contour(x,y,h8_fgen,alpha=0.7,colors='fuchsia',levels=range(4,60,4),linewidths=3)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_purple_line = lines.Line2D([],[],linestyle='dashed',color='purple',label='Temperature (0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax5.legend(handles=[purple,dashed_blue_line,dashed_purple_line,dashed_red_line,black_line],loc=4,framealpha=1)

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
    sub_n = 49
    sub_s = 32

    ax1.set_extent((sub_w1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    #fig.canvas.draw()
    fig.tight_layout()
    plt.savefig(output_dir+'/HRRR_ex/EC_fivepanelwinter_pres_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    #ax1.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=7)
    ax1.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/HRRR_ex/NE_fivepanelwinter_pres_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.clf()

    ### SECOND PLOT ###
    fig2 = plt.figure(figsize=(15,15))
    ax6 = fig2.add_subplot(111, projection = zH5_crs)

    ax6.coastlines(resolution='10m')
    ax6.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax6.add_feature(cfeature.STATES.with_scale('10m'))

    fc = ax6.contour(x,y,h7_fgen,alpha=0.7,colors='fuchsia',levels=range(4,50,2),linewidths=3)

    rhc = ax6.contourf(x,y,rh7, alpha = 0.7, cmap = 'winter_r',transform=zH5_crs, levels=rhlevs,extend='max')
    h7c = ax6.contour(x,y,hgt7,colors='k',levels=range(2700,3300,30),linewidths=2)
    t7c = ax6.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='r',levels=range(0,18,3),linestyles='dashed',linewidths=2)
    t7c2 = ax6.contour(x,y,ndimage.gaussian_filter(t7,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidths=2)
    cbar3 = fig.colorbar(rhc,orientation='vertical',pad=0.01,shrink=.8,ax=ax6,aspect=50,extendrect=False,ticks=[0.7,0.8,0.9,1])
    cbar3.set_label('700mb RH',fontsize=14)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax6.legend(handles=[purple,dashed_blue_line,dashed_red_line,black_line],loc=4,framealpha=1)


    ax6.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax6.set_title('HRRR 700mb Forecast valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/HRRR_ex/EC_700mb_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    ax6.set_extent((283, 295, 39, 49))
    plt.savefig(output_dir+'/HRRR_ex/NE_700mb1_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.clf()
    ### END SECOND PLOT ###

    ### THIRD PLOT ###
    fig3 = plt.figure(figsize=(15,15))
    ax7 = fig3.add_subplot(111,projection=zH5_crs)

    ax7.coastlines(resolution='10m')
    ax7.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax7.add_feature(cfeature.STATES.with_scale('10m'))

    pwatc = ax7.contourf(x,y,pwat, levels=pwlevs, cmap = 'BrBG',alpha=0.5)
    h8c = ax7.contour(x,y,ndimage.gaussian_filter(hgt8,sigma=5,order=0),colors='k',levels=range(1020,1800,30),linewidths=2)
    t8c = ax7.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='r',levels=range(0,30,3),linestyles='dashed',linewidth=2)
    t8c2 = ax7.contour(x,y,ndimage.gaussian_filter(t8,sigma=5,order=0),colors='b',levels=range(-54,-3,3),linestyles='dashed',linewidth=2)
    cbar4 = fig.colorbar(pwatc, orientation='vertical',pad=0.01,ax=ax7, aspect = 50, extendrect=False, ticks = [0.5,1,1.5,2,2.5])
    cbar4.set_label('Precipitable Water (in)')
    ax7.contour(x,y,h8_fgen,alpha=0.7,colors='fuchsia',levels=range(4,50,2),linewidths=3)
    purple = mpatches.Patch(color='fuchsia',label='Frontogenesis')
    dashed_blue_line = lines.Line2D([], [], linestyle='dashed', color='b', label='Temperature (<0C)')
    dashed_red_line = lines.Line2D([], [], linestyle='dashed', color='r', label='Temperature (>0C)')
    black_line = lines.Line2D([], [], color='k', label='Geopotential Height')
    leg = ax7.legend(handles=[purple,dashed_blue_line,dashed_red_line,black_line],loc=4,framealpha=1)


    ax7.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax7.set_title('HRRR 850mb Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/HRRR_ex/EC_850mb_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    ax7.set_extent((283, 295, 39, 49))
    plt.savefig(output_dir+'/HRRR_ex/NE_850mb_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.clf()
