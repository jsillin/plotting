import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from netCDF4 import num2date
import numpy as np
import xarray as xr
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
from metpy.plots import USCOUNTIES
import matplotlib.patches as mpatches


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
init_hour = get_init_hr(hour)
url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr'+mdate+'/hrrr_sfc.t'+get_init_hr(hour)+'z'
#url='http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr20201231/hrrr_sfc.t00z'
print(url)

# Create new directory
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00'
#output_dir = '20201231_0000'
mkdir_p(output_dir)
mkdir_p(output_dir+'/HRRR_ex')
#Parse data using MetPy
ds = xr.open_dataset(url)
init_hr = dt.datetime(int(year),int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(25,55,0.25)
lons = np.arange(260,310,0.25)
total_precip=ds['apcpsfc'].isel(time=0).squeeze()*.0393700787402

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

for i in range(1,49):
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
        'mslmamsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'apcpsfc':'qpf',
        'capesfc':'cape',
        'gustsfc':'sfcgust',
        'hcdchcll':'high_cloud',
        'mcdcmcll':'mid_cloud',
        'lcdclcll':'low_cloud',
        'vissfc':'sfcvis',
        'hgt263_k':'hgt_m10c',
        'hgt253_k':'hgt_m20c',
        'ltngclm':'lightning',
        'sbt124toa':'simsat',
        'hgt0c':'0chgt'
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

    t2m = data['sfc_temp'].squeeze()
    t2mc = t2m-273.15
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2mc = td2m-273.15
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    wb2mc = wet_bulb(t2mc,td2mc)
    #wb2mc = (wb2m-32.)*(5./9.)

    swg = data['sfcgust'].squeeze()
    cloudcover = data['tcc'].squeeze()
    high_cloud = data['high_cloud'].squeeze()
    mid_cloud = data['mid_cloud'].squeeze()
    low_cloud = data['low_cloud'].squeeze()
    vis = data['sfcvis'].squeeze()*0.000621371
    reflectivity = data['radar'].squeeze()
    cape = data['cape'].squeeze()
    lightning=data['lightning'].squeeze()
    dgz_depth = data['hgt_m20c'].squeeze()-data['hgt_m10c'].squeeze()
    simsat = data['simsat'].squeeze()
    hgt0c = data['0chgt'].squeeze()*3.28084
    hrly_precip = data['qpf'].squeeze()*0.0393700787402
    total_precip = total_precip+hrly_precip

    rain = np.ma.masked_where(catrain==0,reflectivity)
    sleet = np.ma.masked_where(catsleet==0,reflectivity)
    ice = np.ma.masked_where(catice==0,reflectivity)
    snow = np.ma.masked_where(catsnow==0,reflectivity)

    qrain = np.ma.masked_where(catrain==0,hrly_precip)
    qsleet = np.ma.masked_where(catsleet==0,hrly_precip)
    qice = np.ma.masked_where(catice==0,hrly_precip)
    qsnow = np.ma.masked_where(catsnow==0,hrly_precip)

    #Generate running accumulation total arrays for each ptype
    acc_snow = acc_snow+qsnow.filled(0)
    acc_sleet = acc_sleet+qsleet.filled(0)
    acc_ice = acc_ice+qice.filled(0)
    acc_rain = acc_rain+qrain.filled(0)

    mslpc = data['mslp'].squeeze()/100
    mslpc=ndimage.gaussian_filter(mslpc,sigma=3,order=0)
    wind_slice = slice(35,-35,35)
    wind_slice_ne = slice(15,-15,15)
    wind_slice_me = slice(10,-10,10)
    u_10m = data['u'].squeeze()
    v_10m = data['v'].squeeze()
    u_10m = u_10m*1.94384449
    v_10m = v_10m*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    ###COMPUTE ICE ACCRETION###

    fram_accretion=fram(qice,wb2mc,wspd)
    fram_accretion=fram_accretion.filled(0)
    acc_fram = acc_fram+fram_accretion

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(44,15))

    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
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
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5), shrink=0.7)
    cbr.set_label('2m Temperature (F)', fontsize = 14)

    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=True, fmt='%i mb', rightside_up=True, use_clabeltext=True)

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

    #cax1 = fig.add_axes([.35,.5,.02,.25])
    #cax2 = fig.add_axes([.35,.75,.02,.25])

    refp = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs) #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])
    capep = ax4.contourf(x, y, cape, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000], alpha = 0.7, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
    lgt = ax4.contour(x,y,lightning,levels=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
    cb = fig.colorbar(capep, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000])
    cb.set_label('CAPE (J/kg)', size='large')

    #ax4.barbs(x[ws2],y[ws2],s5u[ws2,ws2],s5v[ws2,ws2], length = 7)
    blue = mpatches.Patch(color='#6ec6fd', label='Snow')
    orange = mpatches.Patch(color='#e3aa52', label='Sleet')
    pink = mpatches.Patch(color='#ff94aa', label='Freezing Rain')
    green = mpatches.Patch(color='#9cc78d', label='Rain')
    leg = ax1.legend(handles=[blue,orange,pink,green],loc=4,title='Precipitation Type',framealpha=1)
    leg.set_zorder(100)

    swgc = ax5.contourf(x,y,swg, levels=range(0,60,5), cmap = 'hot')
    cbar4 = fig.colorbar(swgc, orientation='vertical',pad=0.01,ax=ax5, aspect = 50, extendrect=False, ticks = np.arange(0,70,10))
    cbar4.set_label('Surface Wind Gust (kts)')

    ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)
    ax1.set_title('Precip Type, 2m Temperature (F), 10m Winds (kts), and MSLP (mb)',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n HRRR Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

    cloud_levs = [40,45,50,55,60,65,70,75,80,85,90,95]
    #high_cols = ['#33c7f4f4','#40b6f0f1','#4da5edee','#5990e8ea','#6678ebed','#7359e6e8','#8030e5e8','#8c18d5d8','#9916c2c5','#9914b5b8','#9912a5a5','#99119579','#990e7a7c','#990d7173','#990c686a','#990b6365']
    #mid_cols = ['#4dbcfac1','#59aef9b4','#66a0f8a7','#7391f79a','#8083f68d','#8c75f57f','#9966f472','#a658f365','#b34af258','#b33bf14b','#b32df03d','#b325e935','#b31ce32d','#b31ad52a','#b319c827','#b317ba25']
    #low_cols = ['#66ea99f4','#73e78cf3','#80e47ef1','#8ce170f0','#99db54ed','#a6d846ec','#b3d539ea','#bfd32be9','#ccd01de7','#ccc617de','#ccba16d0','#ccae14c2','#cca113b4','#cc9511a6','#cc881098','#cc7c0f8a']

    high_cols = ['#c7f4f4','#b6f0f1','#a5edee','#90e8ea','#78ebed','#59e6e8','#30e5e8','#18d5d8','#16c2c5','#14b5b8','#12a5a5','#119579']
    mid_cols = ['#bcfac1','#aef9b4','#a0f8a7','#91f79a','#83f68d','#75f57f','#66f472','#58f365','#4af258','#3bf14b','#2df03d','#25e935']
    low_cols = ['#ea99f4','#e78cf3','#e47ef1','#e170f0','#db54ed','#d846ec','#d539ea','#d32be9','#d01de7','#c617de','#ba16d0','#ae14c2']

    blue = mpatches.Patch(color='#119579', label='High Clouds')
    green = mpatches.Patch(color='#25e935', label='Mid-Level Clouds')
    purple = mpatches.Patch(color='#ae14c2',label='Low-Level Clouds')
    leg = ax2.legend(handles=[blue,green,purple],loc=3,framealpha=1)
    leg.set_zorder(100)
    tccp = ax2.contourf(x,y,cloudcover,cmap='Greys',levels=cloud_levs,alpha=0,extend='max')
    lccp = ax2.contourf(x,y,low_cloud, colors=low_cols,levels=cloud_levs,alpha = 0.35,extend='max')
    mccp = ax2.contourf(x,y,mid_cloud, colors=mid_cols,levels=cloud_levs,alpha = 0.25,extend='max')
    hccp = ax2.contourf(x,y,high_cloud, colors=high_cols,levels=cloud_levs,alpha = 0.15,extend='max')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(tccp,orientation='vertical',pad=0.01,ax=ax2,shrink=.8,aspect=50,extendrect=False, ticks=np.arange(10,100,10))
    cbar2.set_label('Cloud Cover (%)',fontsize=14)

    tprecip = ax3.contourf(x,y,total_precip, alpha = 0.7, cmap = 'cool',transform=zH5_crs, levels=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    tcprecip = ax3.contour(x,y,total_precip,colors=['b','darkblue','darkviolet'],levels=[0.5,1,1.5],linewidths=2)
    cbar3 = fig.colorbar(tprecip,orientation='vertical',pad=0.01,shrink=.8,ax=ax3,aspect=50,extendrect=False,ticks=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    cbar3.set_label('Total Precipitation (inches)',fontsize=14)

    #refp3 = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs)
    #cbar4 = fig.colorbar(tccp,orientation='horizontal',pad=0.01,ax=ax4,aspect=50,extendrect=False)

    #ax1.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax2.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax3.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax4.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #ax5.set_extent((255, 290, 25, 45))#, crs = zH5_crs)    # Set a title and show the plot
    #north=50, south=15, east=-70, west=-115

    sub_w1 = 260
    sub_w = 262
    sub_e = 295
    sub_n = 50
    sub_s = 25

    ax1.set_extent((sub_w1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    #fig.canvas.draw()
    fig.tight_layout()
    plt.savefig(output_dir+'/HRRR_ex/EC_fivepanelwinter_'+dtfs+'_.png')
    ax1.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=7)
    ax1.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/HRRR_ex/NE_fivepanelwinter_'+dtfs+'_.png')
    plt.clf()

    plt.clf()

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
    h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

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
    ax6.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax6.set_title('HRRR Composite Forecast valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/HRRR_ex/EC_ptype_composite_'+dtfs+'_.png')
    ax6.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax6.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=7)
    plt.savefig(output_dir+'/HRRR_ex/NE_ptype_composite_'+dtfs+'_.png')
    ax6.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    wsl = slice(5,-5,5)
    ax6.barbs(x[wsl],y[wsl],u_10m[wsl,wsl],v_10m[wsl,wsl], length=7)
    ax6.set_extent((289,291,43,45))
    plt.savefig(output_dir+'/HRRR_ex/local_ptype_composite_'+dtfs+'_.png')
    plt.close()
    plt.clf()
    plt.close()
    plt.clf()
    ### END SECOND PLOT ###

    ### THIRD PLOT ###
    fig3 = plt.figure(figsize=(15,15))
    ax7 = fig3.add_subplot(111,projection=zH5_crs)

    ax7.coastlines(resolution='10m')
    ax7.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax7.add_feature(cfeature.STATES.with_scale('10m'))

    tprecip = ax7.contourf(x,y,total_precip, alpha = 0.7, cmap = 'cool',transform=zH5_crs, levels=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    tcprecip = ax7.contour(x,y,total_precip,colors=['b','darkblue','darkviolet'],levels=[0.5,1,1.5],linewidths=2)
    cbar3 = fig3.colorbar(tprecip,orientation='vertical',pad=0.01,shrink=.6,ax=ax7,aspect=50,extendrect=False,ticks=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    cbar3.set_label('Total Precipitation (inches)',fontsize=14)

    ax7.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax7.set_title('HRRR Precipitation Forecast Valid Through ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    plt.savefig(output_dir+'/HRRR_ex/EC_total_precip_'+dtfs+'_.png')
    plt.close()
    plt.clf()
    ### FOURTH PLOT ###
    fig4 = plt.figure(figsize=(15,15))
    ax8 = fig4.add_subplot(111,projection=zH5_crs)

    ax8.coastlines(resolution='10m')
    ax8.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax8.add_feature(cfeature.STATES.with_scale('10m'))

    sfcwinds = ax8.contourf(x,y,wspd, cmap='PuRd', levels=range(5,60,5), alpha=0.75)
    cbr = fig4.colorbar(sfcwinds, orientation = 'vertical', pad = 0.01, aspect = 25,
                        panchor = (0.999,0.5), ax = ax8, extendrect=False, ticks = range(5,75,5), shrink = 0.80)
    cbr.set_label('10m Wind Speed (kts)')
    ax8.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=7)
    ax8.set_title('HRRR 10m Wind Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
    ax8.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot

    h_contouqr = ax8.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contouqr.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    plt.savefig(output_dir+'/HRRR_ex/EC_wind_'+dtfs+'_.png')
    ax8.set_extent((281, 295, 39, 49))
    plt.savefig(output_dir+'/HRRR_ex/NE_wind_'+dtfs+'_.png')
    ax8.set_extent((289,291,43,45))
    ax8.barbs(x[wsl],y[wsl],u_10m[wsl,wsl],v_10m[wsl,wsl], length=7)
    ax8.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    plt.savefig(output_dir+'/HRRR_ex/local_wind_'+dtfs+'_.png')
    plt.close()
    plt.clf()
    plt.close()
    plt.clf()

    ### FIFTH PLOT ###
    fig5 = plt.figure(figsize=(15,15))
    ax9 = fig5.add_subplot(111,projection=zH5_crs)

    ax9.coastlines(resolution='10m')
    ax9.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax9.add_feature(cfeature.STATES.with_scale('10m'))
    ax9.set_title('HRRR Cloud Cover Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)

    blue = mpatches.Patch(color='#119579', label='High Clouds')
    green = mpatches.Patch(color='#25e935', label='Mid-Level Clouds')
    purple = mpatches.Patch(color='#ae14c2',label='Low-Level Clouds')
    leg = ax9.legend(handles=[blue,green,purple],loc=3,framealpha=1)
    leg.set_zorder(100)
    tccp = ax9.contourf(x,y,cloudcover,cmap='Greys',levels=cloud_levs,alpha=0,extend='max')
    lccp = ax9.contourf(x,y,low_cloud, colors=low_cols,levels=cloud_levs,alpha = 0.35,extend='max')
    mccp = ax9.contourf(x,y,mid_cloud, colors=mid_cols,levels=cloud_levs,alpha = 0.25,extend='max')
    hccp = ax9.contourf(x,y,high_cloud, colors=high_cols,levels=cloud_levs,alpha = 0.15,extend='max')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(tccp,orientation='vertical',pad=0.01,ax=ax9,shrink=.8,aspect=50,extendrect=False, ticks=np.arange(10,100,10))
    cbar2.set_label('Cloud Cover (%)',fontsize=14)
    ax9.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/HRRR_ex/EC_clouds_'+dtfs+'_.png')
    ax9.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/HRRR_ex/NE_clouds_'+dtfs+'_.png')
    ax9.set_extent((289,291,43,45))
    ax9.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    plt.savefig(output_dir+'/HRRR_ex/local_clouds_'+dtfs+'_.png')
    plt.clf()
    plt.close()
    plt.clf()
    plt.close()

    ###PLOT ACCUMS#####
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

    fig6 = plt.figure(figsize=(25,15))
    ax10 = fig6.add_subplot(221,projection=zH5_crs)
    ax11 = fig6.add_subplot(222,projection=zH5_crs)
    ax12 = fig6.add_subplot(223,projection=zH5_crs)
    ax13 = fig6.add_subplot(224,projection=zH5_crs)

    ax10.coastlines(resolution='10m')
    ax10.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax10.add_feature(cfeature.STATES.with_scale('10m'))

    ax11.coastlines(resolution='10m')
    ax11.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax11.add_feature(cfeature.STATES.with_scale('10m'))

    ax12.coastlines(resolution='10m')
    ax12.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax12.add_feature(cfeature.STATES.with_scale('10m'))

    ax13.coastlines(resolution='10m')
    ax13.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax13.add_feature(cfeature.STATES.with_scale('10m'))


    try:
        ra = ax13.contourf(x,y,acc_rain,colors=qra_cols,levels=qarlevs,alpha=0.7)
        rac = fig6.colorbar(ra,orientation = 'vertical', aspect = 20, ax = ax13, pad = 0.01,
                            extendrect=False, ticks = qarlevs,shrink=0.7)
    except:
        print('no rain')
    try:
        sn = ax10.contourf(x,y,acc_snow,colors=qrs_cols,levels=qaslevs,alpha=0.7)
        snc = fig6.colorbar(sn,orientation = 'vertical', aspect = 20, ax = ax10, pad = 0.01,
                            extendrect=False, ticks = qaslevs,shrink=0.7)
    except:
        print('no snow')
    try:
        ip = ax11.contourf(x,y,acc_sleet,colors=qip_cols,levels=qazlevs,alpha=0.7)
        ipc = fig6.colorbar(ip,orientation = 'vertical', aspect = 20, ax = ax11, pad = 0.01,
                            extendrect=False, ticks = qazlevs,shrink=0.7)
    except:
        print('no sleet')
    try:
        zr = ax12.contourf(x,y,acc_ice, colors=qzr_cols,levels=qazlevs,alpha=0.7)
        zrc = fig6.colorbar(zr,orientation = 'vertical', aspect = 20, ax = ax12, pad = 0.01,
                            extendrect=False, ticks = qazlevs,shrink=0.7)
    except:
        print('no ice')

    ax10.set_title('Accumulated Liquid Equivalent Snow (in)')
    ax10.set_title('HRRR Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=14,loc='left')
    ax11.set_title('Accumulated Liquid Equivalent Sleet (in)')
    ax11.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14,loc='right')
    ax12.set_title('Accumulated Liquid Equivalent Freezing Rain (in)')
    ax13.set_title('Accumulated Rain (in)')
    ax10.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax11.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax12.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    ax13.set_extent((260, 295, 25, 50))#, crs = zH5_crs)    # Set a title and show the plot
    fig6.tight_layout()
    plt.savefig(output_dir+'/HRRR_ex/ec_accum_ptype_'+dtfs+'.png')
    ax10.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax11.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax12.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    ax13.set_extent((281,295,39,49))#, crs = zH5_crs)    # Set a title and show the plot
    plt.savefig(output_dir+'/HRRR_ex/ne_accum_ptype_'+dtfs+'.png')
    ax10.set_extent((289,291,43,45))
    ax10.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    ax11.set_extent((289,291,43,45))
    ax11.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    ax12.set_extent((289,291,43,45))
    ax12.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    ax13.set_extent((289,291,43,45))
    ax13.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    plt.savefig(output_dir+'/HRRR_ex/local_accum_ptype_'+dtfs+'.png')
    plt.close()
    plt.clf()

    fig7 = plt.figure(figsize=(15,15))
    ax14 = fig7.add_subplot(111,projection=zH5_crs)

    ax14.coastlines(resolution='10m')
    ax14.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax14.add_feature(cfeature.STATES.with_scale('10m'))

    try:
        zr = ax14.contourf(x,y,acc_fram,colors=qzr_cols,levels=qazlevs,alpha=0.7,extend='max')
    except:
        print('no ice')
    ax14.set_title('Total Ice Accretion (FRAM)')
    ax14.set_title('HRRR Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=14,loc='left')
    ax14.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14,loc='right')

    ax14.set_extent((260, 295, 25, 50))
    plt.savefig(output_dir+'/HRRR_ex/ec_accum_ice_'+dtfs+'.png')
    ax14.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/HRRR_ex/ne_accum_ice_'+dtfs+'.png')
    ax14.set_extent((289,291,43,45))
    ax14.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    plt.savefig(output_dir+'/HRRR_ex/local_accum_ice_'+dtfs+'.png')
    plt.close()
    plt.clf()

    fig8 = plt.figure(figsize=(15,15))
    ax15 = fig8.add_subplot(111,projection=zH5_crs)

    ax15.coastlines(resolution='10m')
    ax15.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax15.add_feature(cfeature.STATES.with_scale('10m'))

    vis = ax15.contourf(x,y,vis,levels=[0.00625,0.125,0.25,0.5,0.75,1,2,3,4,5],extend='min',colors=['#ffb3ff','magenta','deeppink','hotpink','mediumvioletred','crimson','orangered','darkorange','orange','gold'],alpha=0.7)
    cbr = fig.colorbar(vis,orientation = 'horizontal', aspect = 80, ax = ax15, pad = 0.01,
                        extendrect=False, ticks = [0.00625,0.125,0.25,0.5,0.75,1,2,3,4,5], shrink=0.7)
    cbr.set_label('Surface Visibility (mi)')
    ax15.set_title('Surface Visibility')
    ax15.set_title('HRRR Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=14,loc='left')
    ax15.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14,loc='right')

    ax15.set_extent((260, 295, 25, 50))
    plt.savefig(output_dir+'/HRRR_ex/ec_vis_'+dtfs+'.png')
    ax15.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/HRRR_ex/ne_vis_'+dtfs+'.png')
    ax15.set_extent((289,291,43,45))
    ax15.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    plt.savefig(output_dir+'/HRRR_ex/local_vis_'+dtfs+'.png')
    plt.close()
    plt.clf()

    fig9 = plt.figure(figsize=(15,15))
    ax16 = fig9.add_subplot(111,projection=zH5_crs)

    ax16.coastlines(resolution='10m')
    ax16.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax16.add_feature(cfeature.STATES.with_scale('10m'))

    hgt_frz = ax16.contourf(x,y,hgt0c,cmap='rainbow',levels=range(0,10000,500),extend='max')
    cbar2 = fig9.colorbar(hgt_frz,orientation='vertical', pad = 0.01, aspect = 50, ax = ax16, extendrect=False,shrink=0.7)
    cbr.set_label('Freezing Level Height (ft)')

    ax16.set_title('Freezing Level Height')
    ax16.set_title('HRRR Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=14,loc='left')
    ax16.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14,loc='right')

    ax16.set_extent((260, 295, 25, 50))
    plt.savefig(output_dir+'/HRRR_ex/ec_frz_height_'+dtfs+'.png')
    ax16.set_extent((281,295,39,49))
    plt.savefig(output_dir+'/HRRR_ex/ne_frz_height_'+dtfs+'.png')
    ax16.set_extent((289,291,43,45))
    ax16.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
    plt.savefig(output_dir+'/HRRR_ex/local_frz_height_'+dtfs+'.png')
    plt.close()
    plt.clf()

    print(str(i)+'_Done!')
