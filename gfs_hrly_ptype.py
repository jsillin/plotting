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

# make unique directory to store output
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

#grabbing data from NOMADS
startTime=datetime.now()

m_date='20200903'
m_hour='12'

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

url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs'+mdate+'/gfs_0p25_1hr_'+get_init_hr(hour)+'z'
init_hour = get_init_hr(hour)

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFS')

#Parse data using MetPy
ds = xr.open_dataset(url)
init_hr = dt.datetime(year,int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(25,55,0.25)
lons = np.arange(260,310,0.25)

ds = ds.sel(lat = lats, lon = lons)

#Initialize ptype arrays by grabbing first hour of categorical precip
catrain = ds['crainsfc'].squeeze().isel(time=1).squeeze()
catsnow = ds['csnowsfc'].squeeze().isel(time=1).squeeze()
catsleet = ds['cicepsfc'].squeeze().isel(time=1).squeeze()
catice = ds['cfrzrsfc'].squeeze().isel(time=1).squeeze()

total_precip=ds['apcpsfc'].isel(time=1).squeeze()*.0393700787402

acc_rain = np.ma.masked_where(catrain==0,total_precip)
acc_sleet = np.ma.masked_where(catsleet==0,total_precip)
acc_ice = np.ma.masked_where(catice==0,total_precip)
acc_snow = np.ma.masked_where(catsnow==0,total_precip)

acc_rain = acc_rain.filled(0)
acc_sleet = acc_sleet.filled(0)
acc_ice = acc_ice.filled(0)
acc_snow = acc_snow.filled(0)

t2mi = ds['tmp2m'].isel(time=1).squeeze()-273.15
td2mi = ds['tmp2m'].isel(time=1).squeeze()-273.15

u10 = ds['ugrd10m'].isel(time=1).squeeze()*1.94384449
v10 = ds['vgrd10m'].isel(time=1).squeeze()*1.94384449
ws10 = ((u10**2)+(v10**2))**.5
acc_fram = fram(acc_ice,wet_bulb(t2mi,td2mi),ws10)
print("INITIALIZATION SUCCESSFUL")
#Now loop through the rest to come up with a ptype map each hour and
#ideally a four-panel accumulation map
for i in range(2,120):
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

    #This extends each ptype one gridpoint outwards to prevent a gap between
    #different ptypes
    radius = 1
    kernel = np.zeros((2*radius+1,2*radius+1))
    y1,x1 = np.ogrid[-radius:radius+1,-radius:radius+1]
    mask=x1**2+y1**2 <=radius**2
    kernel[mask]=1

    snowc= gf(catsnow,np.max,footprint=kernel)
    rainc = gf(catrain,np.max,footprint=kernel)

    #Coordinate stuff
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['height'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    zH5_crs = data['temperature'].metpy.cartopy_crs

    #data['temperature'].metpy.convert_units('degC')
    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    wb2m = wet_bulb(t2m,td2m)
    wb2mc = (wb2m-32.)*(5./9.)

    cloudcover = data['tcc'].squeeze()
    reflectivity = data['radar'].squeeze()
    hrly_precip = data['qpf'].squeeze()*0.0393700787402
    new_precip = hrly_precip-total_precip
    total_precip = hrly_precip

    rain = np.ma.masked_where(rainc==0,new_precip)
    sleet = np.ma.masked_where(catsleet==0,new_precip)
    ice = np.ma.masked_where(catice==0,new_precip)
    snow = np.ma.masked_where(snowc==0,new_precip)

    #Generate running accumulation total arrays for each ptype
    acc_snow = acc_snow+snow.filled(0)
    acc_sleet = acc_sleet+sleet.filled(0)
    acc_ice = acc_ice+ice.filled(0)
    acc_rain = acc_rain+rain.filled(0)


    #Smooth rain
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
    ws10m = ((u_10m**2)+(v_10m**2))**.5

    ###COMPUTE ICE ACCRETION###

    fram_accretion=fram(ice,wb2mc,ws10m)
    fram_accretion=fram_accretion.filled(0)
    acc_fram = acc_fram+fram_accretion

    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))


    ########## PLOTTING #######################################################
    tmp_2m = ax1.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5))
    cbr.set_label('2m Temperature (F)', fontsize = 14)

    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    q_levs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25]
    qarlevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10]
    qaslevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5]
    qazlevs = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3]
    qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
    qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
    qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
    qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']
    qra_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00','#ffcc00','#ff9500','#ff4800','#ff2900','#ff1200','#ff0000','#cc0000','#990000','#990033','#b3003b','#ff3333','#ff6666','#ffffff']
    qrs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5','#4f01f6','#7a00f5','#9e00f5','#b833ff','#d280ff','#cc00f1','#ad00cc','#820099','#4700b3']
    qzr_cols = ['#ff0066','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333','#ff6666','#ff9999','#ffcccc','#ffffff']
    qip_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc','#9933ff','#bf80ff','#e6ccff','#ffffff']
    q_levs_r = [0.01,0.05,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3]
    try:
        ra = ax1.contourf(x,y,rain,colors=qr_cols,levels=q_levs,alpha=0.7)
    except:
        print('no rain')
    try:
        sn = ax1.contourf(x,y,snow,colors=qs_cols,levels=q_levs,alpha=0.7)
    except:
        print('no snow')
    try:
        ip = ax1.contourf(x,y,sleet,colors=qi_cols,levels=q_levs,alpha=0.7)
    except:
        print('no sleet')
    try:
        zr = ax1.contourf(x,y,ice, colors=qz_cols,levels=q_levs,alpha=0.7)
    except:
        print('no ice')
    ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)
    ax1.set_title('Precip Type, 2m Temperature (F), 10m Winds (kts), and MSLP (mb)',fontsize=16)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax1.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GFS/gfs_hrly_pytpe_v7_'+dtfs+'.png')
    fcst_hr = str(0)
    plt.close()
    plt.clf()

    fig2 = plt.figure(figsize=(25,15))
    ax2 = fig2.add_subplot(221,projection=zH5_crs)
    ax3 = fig2.add_subplot(222,projection=zH5_crs)
    ax4 = fig2.add_subplot(223,projection=zH5_crs)
    ax5 = fig2.add_subplot(224,projection=zH5_crs)

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

    fig2.suptitle('GFS Accumulated Precipitation By Precipitation Type \n Initialized: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' \n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item()+' | '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14)

    try:
        ra = ax5.contourf(x,y,acc_rain,colors=qra_cols,levels=qarlevs,alpha=0.7)
        rac = fig2.colorbar(ra,orientation = 'vertical', aspect = 20, ax = ax5, pad = 0.01,
                            extendrect=False, ticks = [0.5,1,1.5,2,3,4,5,7,9],shrink=0.7)
    except:
        print('no rain')
    try:
        sn = ax2.contourf(x,y,acc_snow,colors=qrs_cols,levels=qaslevs,alpha=0.7)
        snc = fig2.colorbar(sn,orientation = 'vertical', aspect = 20, ax = ax2, pad = 0.01,
                            extendrect=False, ticks = range(0,5,1),shrink=0.7)
    except:
        print('no snow')
    try:
        ip = ax3.contourf(x,y,acc_sleet,colors=qip_cols,levels=qazlevs,alpha=0.7)
        ipc = fig2.colorbar(ip,orientation = 'vertical', aspect = 20, ax = ax3, pad = 0.01,
                            extendrect=False, ticks = [0.01,0.1,0.5,1,1.5,2,3],shrink=0.7)
    except:
        print('no sleet')
    try:
        zr = ax4.contourf(x,y,acc_ice, colors=qzr_cols,levels=qazlevs,alpha=0.7)
        zrc = fig2.colorbar(zr,orientation = 'vertical', aspect = 20, ax = ax4, pad = 0.01,
                            extendrect=False, ticks = [0.01,0.1,0.5,1,1.5,2,3],shrink=0.7)
    except:
        print('no ice')

    ax2.set_title('Accumulated Liquid Equivalent Snow')
    ax3.set_title('Accumulated Liquid Equivalent Sleet')
    ax4.set_title('Accumulated Liquid Equivalent Freezing Rain')
    ax5.set_title('Accumulated Rain')
    ax2.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    fig2.tight_layout()
    plt.savefig(output_dir+'/GFS/ec_accum_ptype_'+dtfs+'.png')
    ax2.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax3.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax4.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax5.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    plt.savefig(output_dir+'/GFS/ne_accum_ptype_'+dtfs+'.png')
    plt.close()
    plt.clf()

    fig3 = plt.figure(figsize=(15,15))
    ax6 = fig3.add_subplot(111,projection=zH5_crs)
    try:
        zr_acc = ax6.contourf(x,y,acc_fram,colors=qzr_cols,levels=qazlevs,alpha=0.8)
        zrc = fig3.colorbar(zr_acc,orientation = 'horizontal', aspect = 80, ax = ax6, pad = 0.01,
                                extendrect=False, ticks = [0.01,0.1,0.5,1,1.5,2,3],shrink=0.9)
        zrc.set_label('Inches of Horizontal Accretion')
    except:
        print('no ice')

    ax6.coastlines(resolution='10m')
    ax6.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax6.add_feature(cfeature.STATES.with_scale('10m'))
    ax6.set_title('GFS Ice Accretion Forecast (FRAM)',fontsize=14)
    ax6.set_title(' \n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax6.set_title(' \n Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax6.set_extent((265, 300, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/GFS/ec_ice_accretion_fram_'+dtfs+'.png')
    ax6.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='gray',linewidths=0.5)
    ax6.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/GFS/ne_ice_accretion_fram_'+dtfs+'.png')
    ax6.set_extent((260,271,36,42))#, crs = zH5_crs)    # Set a title and show the plot
    #plt.savefig(output_dir+'/GFS/kc_outer_ice_accretion_fram_'+dtfs+'.png')
    #ax6.set_extent((263,268,37.5,40.5))#, crs = zH5_crs)    # Set a title and show the plot
    #plt.savefig(output_dir+'/GFS/kc_ice_accretion_fram_'+dtfs+'.png')
    plt.close()
    plt.clf()

    print('Hour '+str(i)+' completed!')
    timeelapsed = datetime.now()-startTime
    print(timeelapsed)
