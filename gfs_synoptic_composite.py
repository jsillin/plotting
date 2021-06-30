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
import supplementary_tools as spt
import metpy.calc as mpcalc
import matplotlib.lines as lines
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def addcolorbar(var,clevs):
    axes_bbox = ax.get_position()
    left = axes_bbox.x1 + 0.015
    bottom = axes_bbox.y0
    width = 0.015
    height = axes_bbox.y1 - bottom
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(var, cax=cax, ticks=clevs, orientation='vertical')
    cbar.ax.tick_params(labelsize=10)
    #cbar.set_label('500-hPa temperature [K]', size=11)  # MODIFY THIS for other fields!!

startTime=datetime.now()

mdate = spt.get_init_time('GFS')[0]
init_hour = spt.get_init_time('GFS')[1]
url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25/gfs'+mdate+'/gfs_0p25_'+init_hour+'z'

# Create new directory
output_dir = mdate+'_'+init_hour+'00'
#output_dir = '20201231_0000'
spt.mkdir_p(output_dir)
spt.mkdir_p(output_dir+'/GFS')

#Parse data using MetPy
ds = xr.open_dataset(url)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(15,70,0.25)
lons = np.arange(220,310,0.25)

for i in range(1,64):
    data = ds.metpy.parse_cf()
    data = data.isel(time=i)
    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'absvprs':'avort',
        'hgtprs':'gph',
        'rhprs':'rh',
        'tmpprs':'temp',
        'ugrdprs':'u',
        'vgrdprs': 'v',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf'
    })

    #This extends each ptype one gridpoint outwards to prevent a gap between
    #different ptypes
    radius = 1
    kernel = np.zeros((2*radius+1,2*radius+1))
    y1,x1 = np.ogrid[-radius:radius+1,-radius:radius+1]
    mask=x1**2+y1**2 <=radius**2
    kernel[mask]=1

    ########### GRABBING AND PARSING DATA #################
    catrain = data['catrain'].sel(lat=lats,lon=lons).squeeze()
    catsnow = data['catsnow'].sel(lat=lats,lon=lons).squeeze()
    catsleet = data['catsleet'].sel(lat=lats,lon=lons).squeeze()
    catice = data['catice'].sel(lat=lats,lon=lons).squeeze()

    snowc= gf(catsnow,np.max,footprint=kernel)
    icec = gf(catice,np.max,footprint=kernel)
    sleetc = gf(catsleet,np.max,footprint=kernel)
    rainc = gf(catrain,np.max,footprint=kernel)

    vertical, = data['temp'].metpy.coordinates('vertical')
    time = data['temp'].metpy.time
    zH5_crs = data['temp'].metpy.cartopy_crs

    t2m = data['sfc_temp'].sel(lat=lats,lon=lons).squeeze()
    t2mc = t2m-273.15
    t2m = ((t2m - 273.15)*(9./5.))+32.

    t2 = data['temp'].sel(lev=250.0,lat=lats,lon=lons)-273.15
    u2 = data['u'].sel(lev=250.0,lat=lats,lon=lons).squeeze()*1.94384449
    v2 = data['v'].sel(lev=250.0,lat=lats,lon=lons).squeeze()*1.94384449
    av2 = data['avort'].sel(lev=250.0,lat=lats,lon=lons).squeeze()*1e5
    h2 = data['gph'].sel(lev=250.0,lat=lats,lon=lons).squeeze()

    t5 = data['temp'].sel(lev=500.0,lat=lats,lon=lons)-273.15
    u5 = data['u'].sel(lev=500.0,lat=lats,lon=lons).squeeze()*1.94384449
    v5 = data['v'].sel(lev=500.0,lat=lats,lon=lons).squeeze()*1.94384449
    av5 = data['avort'].sel(lev=500.0,lat=lats,lon=lons).squeeze()*1e5
    rh5 = data['rh'].sel(lev=500.0,lat=lats,lon=lons).squeeze()
    h5 = data['gph'].sel(lev=500.0,lat=lats,lon=lons).squeeze()

    t7 = data['temp'].sel(lev=700.0,lat=lats,lon=lons).squeeze()
    u7 = data['u'].sel(lev=700.0,lat=lats,lon=lons).squeeze()
    v7 = data['v'].sel(lev=700.0,lat=lats,lon=lons).squeeze()
    h7 = data['gph'].sel(lev=700.0,lat=lats,lon=lons).squeeze()
    rh7 = data['rh'].sel(lev=700.0,lat=lats,lon=lons).squeeze()

    t8 = data['temp'].sel(lev=850.0,lat=lats,lon=lons).squeeze()
    u8 = data['u'].sel(lev=850.0,lat=lats,lon=lons).squeeze()
    v8 = data['v'].sel(lev=850.0,lat=lats,lon=lons).squeeze()
    h8 = data['gph'].sel(lev=850.0,lat=lats,lon=lons).squeeze()
    rh8 = data['rh'].sel(lev=850.0,lat=lats,lon=lons).squeeze()

    t9 = data['temp'].sel(lev=925.0,lat=lats,lon=lons).squeeze()
    u9 = data['u'].sel(lev=925.0,lat=lats,lon=lons).squeeze()
    v9 = data['v'].sel(lev=925.0,lat=lats,lon=lons).squeeze()
    h9 = data['gph'].sel(lev=925.0,lat=lats,lon=lons).squeeze()
    rh9 = data['rh'].sel(lev=925.0,lat=lats,lon=lons).squeeze()

    reflectivity = data['radar'].sel(lat=lats,lon=lons).squeeze()
    rain = np.ma.masked_where(rainc==0,reflectivity)
    sleet = np.ma.masked_where(sleetc==0,reflectivity)
    ice = np.ma.masked_where(icec==0,reflectivity)
    snow = np.ma.masked_where(snowc==0,reflectivity)

    mslp = data['mslp']/100.
    mslpc = mslp.sel(lat=lats,lon=lons).squeeze()
    mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)

    #av5 = ndimage.gaussian_filter(av5,sigma=1,order=0)

    x, y = t7.metpy.coordinates('x', 'y')
    dx, dy = mpcalc.grid_deltas_from_dataarray(t7)
    lat, lon = xr.broadcast(y, x)
    wind_slice = slice(7,-7,7)
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

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
    #h7_adv = mpcalc.advection(t7, [u7, v7], (dx, dy), dim_order='yx') * units('K/sec')

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
    #h8_adv = mpcalc.advection(t8, [u8, v8], (dx, dy), dim_order='yx') * units('K/sec')

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
    #h9_adv = mpcalc.advection(t9, [u9, v9], (dx, dy), dim_order='yx') * units('K/sec')

    x, y = t2.metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    wind_slice = slice(5,-5,5)

    h2_wspd = ((u2**2)+(v2**2))**.5

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
    qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
    qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
    qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']
    qra_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00','#ffcc00','#ff9500','#ff4800','#ff2900','#ff1200','#ff0000','#cc0000','#990000','#990033','#b3003b','#ff3333','#ff6666','#ffffff']
    qrs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5','#4f01f6','#7a00f5','#9e00f5','#b833ff','#d280ff','#cc00f1','#ad00cc','#820099','#4700b3']
    qzr_cols = ['#ff0066','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333','#ff6666','#ff9999','#ffcccc','#ffffff']
    qip_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc','#9933ff','#bf80ff','#e6ccff','#ffffff']


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
    ax.set_title('GFS Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    pink = lines.Line2D([],[],linestyle='solid',color='fuchsia',label='Frontogenesis (K/100km/3hrs)')
    red = lines.Line2D([],[],linestyle='dashed',color='r',label='>=0C Temperatures (C)')
    blue = lines.Line2D([],[],linestyle='dashed',color='b',label='<0C Temperatures (C)')
    gray = lines.Line2D([],[],linestyle='solid',color='dimgray',label='Height (m)')
    leg = ax.legend(handles=[pink,gray,red,blue],framealpha=1,loc=4)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GFS/ec_h7fgen_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    ax.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFS/ne_h7fgen_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)

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
    ax2.set_title('GFS Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax2.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    pink = lines.Line2D([],[],linestyle='solid',color='fuchsia',label='Frontogenesis (K/100km/3hrs)')
    red = lines.Line2D([],[],linestyle='dashed',color='r',label='>=0C Temperatures (C)')
    blue = lines.Line2D([],[],linestyle='dashed',color='b',label='<0C Temperatures (C)')
    gray = lines.Line2D([],[],linestyle='solid',color='dimgray',label='Height (m)')
    leg = ax2.legend(handles=[pink,gray,red,blue],framealpha=1,loc=4)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GFS/ec_h8fgen_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    ax2.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFS/ne_h8fgen_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)

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
    ax3.contour(x,y,t9c,colors='b',levels=range(-40,0,2),linewidths=1.5,linestyles='dashed')
    ax3.contour(x,y,t9c,colors='r',levels=range(0,40,2),linewidths=1.5,linestyles='dashed')
    ax3.set_title('925mb Forecast Summary',fontsize=16)
    ax3.set_title('GFS Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax3.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    pink = lines.Line2D([],[],linestyle='solid',color='fuchsia',label='Frontogenesis (K/100km/3hrs)')
    red = lines.Line2D([],[],linestyle='dashed',color='r',label='>=0C Temperatures (C)')
    blue = lines.Line2D([],[],linestyle='dashed',color='b',label='<0C Temperatures (C)')
    gray = lines.Line2D([],[],linestyle='solid',color='dimgray',label='Height (m)')
    leg = ax3.legend(handles=[pink,gray,red,blue],framealpha=1,loc=4)
    leg.set_zorder(100)
    plt.savefig(output_dir+'/GFS/ec_h9fgen_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    ax3.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFS/ne_h9fgen_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)

    fig4 = plt.figure(figsize=(15,15))

    ax4 = fig4.add_subplot(111,projection=zH5_crs)
    ax4.coastlines(resolution='10m')
    ax4.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax4.add_feature(cfeature.STATES.with_scale('10m'))

    h5c = ax4.contour(x,y,h5,colors='whitesmoke', levels = range(4800,6200,60),linewidths=1.5)
    a5c = ax4.contourf(x,y,av5,cmap='autumn_r',levels=range(15,50,2),alpha=0.3,extend='max',antialiased=True)
    w2c = ax4.contourf(x,y,h2_wspd,cmap='PuRd',levels=range(80,200,20),alpha=0.3)
    msl = ax4.contour(x, y, mslpc, colors='lightgray', levels=range(940,1040,4),linewidths=1,alpha=0.6)
    tmp_2m32 = ax4.contour(x,y,t2m,colors='deeppink', alpha = 0.8, levels = [32])
    t8w = ax4.contour(x,y,t8c,colors='lightcoral',linestyle='dashed',levels=range(5,50,5),alpha=0.6,linewidths=1.5)
    t8c0 = ax4.contour(x,y,t8c,colors='magenta',linestyle='dashed',levels=[0],alpha=0.6,linewidths=1.5)
    t8c = ax4.contour(x,y,t8c,colors='turquoise',linestyle='dashed',levels=range(-50,0,5),alpha=0.6,linewidths=1.5)

    try:
        ra = ax4.contourf(x,y,rain,colors=qr_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no rain')
    try:
        sn = ax4.contourf(x,y,snow,colors=qs_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no snow')
    try:
        ip = ax4.contourf(x,y,sleet,colors=qi_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no sleet')
    try:
        zr = ax4.contourf(x,y,ice, colors=qz_cols,levels=ref_levs,alpha=0.4,extend='max')
    except:
        print('no ice')

    ax4.set_extent((230, 295, 20, 60))
    ax4.set_title('Synoptic Forecast Summary',fontsize=16)
    ax4.set_title('GFS Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=12,loc='left')
    ax4.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=12,loc='right')
    ax4.set_facecolor('dimgray')

    #cbax1 = fig4.add_axes([0.65,0.05,0.3,0.1])
    #sncb = fig4.colorbar(sn,cax=cbax1,orientation='horizontal',ticks=ref_levs)
    #addcolorbar(sn,ref_levs)
    plt.savefig(output_dir+'/GFS/gfs_synoptic_composite_ec_v13_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)

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
