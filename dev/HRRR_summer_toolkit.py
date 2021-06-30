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
import supplementary_tools as spt
import matplotlib.colors as col
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import hrrr_dash_cbars as cbars

######### CONFIGURE PLOTTING SEPARATE FIELDS #############
plot_precip = False
plot_winds = False
plot_vis = False
plot_clouds = False
plot_overview = True

######## CONFIGURE EXTRA DOMAIN OF INTEREST ###############
extradoms = True
centerlat = [33,40,37,32.7]
centerlon = [99,101,87,88.3]
state_abbr = ['tx','ks','ky','al']

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

total_precip=ds['apcpsfc'].isel(time=1).squeeze()*.0393700787402

t2mi = ds['tmp2m'].isel(time=1).squeeze()-273.15
td2mi = ds['tmp2m'].isel(time=1).squeeze()-273.15

u10 = ds['ugrd10m'].isel(time=1).squeeze()*1.94384449
v10 = ds['vgrd10m'].isel(time=1).squeeze()*1.94384449
ws10 = ((u10**2)+(v10**2))**.5
print("INITIALIZATION SUCCESSFUL")

for i in range(0,49):

    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    #Rename variables to useful things
    data = data.rename({
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
        'hgt0c':'0chgt',
        'hgt500mb':'height5',
        'hgt700mb':'height7',
        'pwatclm':'pwat',
        'dptprs':'dpt',
        'ugrdprs':'uprs',
        'vgrdprs':'vprs',
        'cape0_3000m':'3cape',
        'hgtl5':'lclh',
        'hlcy1000_0m':'1kmhelicity',
        'hlcy3000_0m':'3kmhelicity',
        'ustm0_6000m':'u_storm',
        'vstm0_6000m':'v_storm',
        'pot2m':'2mthetae',
        'vucsh0_6000m':'ushr6km',
        'vvcsh0_6000m':'vshr6km',
        'vucsh0_1000m':'ushr1km',
        'vvcsh0_1000m':'vshr1km'
    })

    zH5 = data['temperature'].squeeze()
    zH5_crs = zH5.metpy.cartopy_crs

    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    dx, dy = mpcalc.lat_lon_grid_deltas(lat,lon)

    t2m = data['sfc_temp'].squeeze()
    t2mc = t2m-273.15
    t2m = ((t2m - 273.15)*(9./5.))+32.

    td2m = data['sfc_td'].squeeze()
    td2mc = td2m-273.15
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    wb2mc = spt.wet_bulb(t2mc,td2mc)

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

    cape3 = data['3cape'].squeeze()
    helicity_1km = data['1kmhelicity'].squeeze()
    helicity_3km = data['3kmhelicity'].squeeze()

    u_storm = data['u_storm'].squeeze()*1.94384449
    v_storm = data['v_storm'].squeeze()*1.94384449

    u1km_shr = data['ushr1km'].squeeze()*1.94384449
    v1km_shr = data['vshr1km'].squeeze()*1.94384449

    u6km_shr = data['ushr6km'].squeeze()*1.94384449
    v6km_shr = data['vshr6km'].squeeze()*1.94384449

    capesmoothed = ndimage.gaussian_filter(cape,sigma=2,order=0)

    thetae_2m = data['2mthetae'].squeeze()
    lcls = data['lclh'].squeeze()
    lcl_height = np.ma.masked_where(cape<25,lcls)
    lclsmoothed = ndimage.gaussian_filter(lcls,sigma=2,order=0)
    mslpc = data['mslp'].squeeze()/100
    mslpc=ndimage.gaussian_filter(mslpc,sigma=3,order=0)

    x_2d, y_2d = np.meshgrid(x, y)
    quiver_slices = (slice(None,None,50),slice(None,None,50))
    quiver_slices_reg = (slice(None,None,30),slice(None,None,30))
    wind_slice_zo = slice(60,-60,60)
    wind_slice = slice(36,-36,36)
    wind_slice_ne = slice(18,-18,18)
    wind_slice_me = slice(9,-9,9)

    u_10m = data['u'].squeeze()
    v_10m = data['v'].squeeze()
    u_10m = u_10m*1.94384449
    v_10m = v_10m*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5
    print(u_10m)
    print(v_10m)
    print(u_10m.values)
    print(v_10m.values)
    print(dx)
    print(dy)
    sfcvort = mpcalc.vorticity(u_10m.values*units.knots,v_10m.values*units.knots,dx,dy)
    print(sfcvort)
    print(np.max(sfcvort))
    sfcvort = sfcvort*1e5
    print(sfcvort)
    print(np.max(sfcvort))

    t500 = data['temperature'].sel(lev=500).squeeze()
    t700 = data['temperature'].sel(lev=700).squeeze()
    h500 = data['height5'].squeeze()
    h700 = data['height7'].squeeze()

    delt = t700-t500
    delp = h700-h500
    mllapse = (delt/delp)*1000
    mllapse = ndimage.gaussian_filter(mllapse,sigma=1,order=0)

    u3 = data['uprs'].sel(lev=300).squeeze()*1.94384449
    v3 = data['vprs'].sel(lev=300).squeeze()*1.94384449
    u5 = data['uprs'].sel(lev=500).squeeze()*1.94384449
    v5 = data['vprs'].sel(lev=500).squeeze()*1.94384449
    u7 = data['uprs'].sel(lev=700).squeeze()*1.94384449
    v7 = data['vprs'].sel(lev=700).squeeze()*1.94384449
    u8 = data['uprs'].sel(lev=850).squeeze()*1.94384449
    v8 = data['vprs'].sel(lev=850).squeeze()*1.94384449
    u9 = data['uprs'].sel(lev=925).squeeze()*1.94384449
    v9 = data['vprs'].sel(lev=925).squeeze()*1.94384449
    u10 = data['uprs'].sel(lev=1000).squeeze()*1.94384449
    v10 = data['vprs'].sel(lev=1000).squeeze()*1.94384449

    ws3 = ((u3**2)+(v3**2))**.5
    ws5 = ((u5**2)+(v5**2))**.5
    ws8 = ((u8**2)+(v8**2))**.5

    u3m = np.ma.masked_where(ws3<35,u3)
    v3m = np.ma.masked_where(ws3<35,v3)

    u5m = np.ma.masked_where(ws5<35,u5)
    v5m = np.ma.masked_where(ws5<35,v5)

    u8m = np.ma.masked_where(ws8<25,u8)
    v8m = np.ma.masked_where(ws8<25,v8)

    umean = (u5+u7+u8+u9+u10)/5
    vmean = (v5+v7+v8+v9+v10)/5

    pwat = data['pwat'].squeeze()*0.0393700787402

    uivt = umean*(pwat/0.0393700787402)
    vivt = vmean*(pwat/0.0393700787402)

    ###GRAB LOCAL DATA###
    stations=['PWM','AUG','RKD','SFM','IZG','LEW','EEN','MHT','BML','CVA','1P1','FAR','BGR','GNR','CAR','MAC',
            'RUT','BTV','CDA','ORH','TAN','BAF','BOS','ALB','SYR','ROC','ITH','ART','MSV','HFD','BDL','OKX',
            'PVD','YSC']
    coords=[[43.644940, -70.309360],[44.317911, -69.796462],[44.062355, -69.095368],[43.398890, -70.711034],
            [43.989267, -70.946733],[44.047577, -70.284566],[42.901405, -72.269613],[42.929656, -71.434166],
            [44.577748, -71.177514],[45.084798, -70.217251],[43.778606, -71.753213],[44.670005, -70.151436],
            [44.807750, -68.818231],[45.463921, -69.553345],[45.463921, -69.553345],[44.727026, -67.473533],
            [43.527006, -72.949843],[44.466761, -73.141541],[44.570582, -72.017706],[44.466761, -73.141541],
            [42.266704, -71.873740],[41.873395, -71.015191],[42.159667, -72.715977],[42.358814, -71.057869],
            [42.749352, -73.804664],[43.112376, -76.110738],[43.122772, -77.672485],[42.490250, -76.457773],
            [43.992712, -76.022970],[41.702975, -74.798589],[41.737535, -72.649142],[41.478649, -73.135379],
            [41.724628, -71.427676],[45.438474, -71.691100]]
    pwlevs = []
    for j in range(3,13):
        lev = 0.2*j
        pwlevs.append(lev)

    station_qpf = []
    station_temp = []

    for i in range(len(stations)):
        station = stations[i]
        lat = coords[i][0]
        lon = coords[i][1]
        qpf = total_precip.interp(lat=lat,lon=lon).values
        temp = t2m.interp(lat=lat,lon=lon).values
        station_qpf.append(np.round(qpf,1))
        station_temp.append(np.round(temp,1))

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(44,15))

    gs = fig.add_gridspec(ncols=3,nrows= 2, width_ratios=[1,2,1])
    gs.update(hspace=0.01,wspace=0.01)
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

    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    print(dtfs)
    ########## PLOTTING #######################################################
    wsl = slice(5,-5,5)
    quiver_kwargs = {'headlength': 4, 'headwidth': 3, 'angles': 'uv', 'scale_units': 'xy',
                 'scale': 10,'linewidth':2}
    tmp_2m = ax1.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])


    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=True, fmt='%i mb', rightside_up=True, use_clabeltext=True)
    c_contour = ax1.contour(x,y,capesmoothed,colors=['plum','violet','deeppink','magenta'],levels=[500,1000,1500,2000],linewidths=1.5)
    ref_levs = range(5,75,1)

    norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
    newcmap = ListedColormap(cmap_ref(range(40, 194)))
    newnorm = Normalize(0,77)
    refc = ax1.contourf(x,y,reflectivity,norm=newnorm,cmap=newcmap,levels=range(5,75,1),alpha=0.7)
    cbars.addt2mcolorbar(ax1,fig,tmp_2m,range(-20,100,5))
    cbars.addrefcolorbar(ax1,fig,refc,ref_levs)

    tdc = ax1.contour(x,y,td2m,colors=['springgreen','lawngreen','limegreen','green'],levels=[55,60,65,70],linewidths=1.5)
    tdf = ax1.contourf(x,y,td2m,colors=['springgreen','lawngreen','limegreen','green'],levels=[55,60,65,70,100],alpha=0.15)

    ax1.quiver(x_2d[quiver_slices],y_2d[quiver_slices],u_storm[quiver_slices],v_storm[quiver_slices],color="white")

    ax1.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=6)
    ax1.set_title('Reflectivity (dBZ), 2m Temp (F), 2m Dewpoint, 10m Wind (kts), MSLP (mb), SBCAPE (J/kg), and Mean Storm Motion (kts)',fontsize=14)
    ax1.set_title('\n Valid: '+time.dt.strftime('%a %b %d %HZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n HRRR Init: '+init_time.dt.strftime('%Y-%m-%d %HZ').item(),fontsize=11,loc='left')

    blue = mpatches.Patch(color='green', label='>70F')
    orange = mpatches.Patch(color='limegreen', label='65-70F')
    pink = mpatches.Patch(color='lawngreen', label='60-65F')
    green = mpatches.Patch(color='springgreen', label='55-60F')
    leg = ax1.legend(handles=[blue,orange,pink,green],loc=4,title='2m Dew Point',framealpha=1)
    leg.set_zorder(100)

    #### TOP LEFT PANEL #########
    cloud_levs = [40,45,50,55,60,65,70,75,80,85,90,95]
    high_cols = ['#c7f4f4','#b6f0f1','#a5edee','#90e8ea','#78ebed','#59e6e8','#30e5e8','#18d5d8','#16c2c5','#14b5b8','#12a5a5','#119579']
    mid_cols = ['#bcfac1','#aef9b4','#a0f8a7','#91f79a','#83f68d','#75f57f','#66f472','#58f365','#4af258','#3bf14b','#2df03d','#25e935']
    low_cols = ['#ea99f4','#e78cf3','#e47ef1','#e170f0','#db54ed','#d846ec','#d539ea','#d32be9','#d01de7','#c617de','#ba16d0','#ae14c2']

    tccp = ax2.contourf(x,y,cloudcover,cmap='Greys',levels=cloud_levs,alpha=0,extend='max')
    lccp = ax2.contourf(x,y,low_cloud, colors=low_cols,levels=cloud_levs,alpha = 0.35,extend='max')
    mccp = ax2.contourf(x,y,mid_cloud, colors=mid_cols,levels=cloud_levs,alpha = 0.25,extend='max')
    hccp = ax2.contourf(x,y,high_cloud, colors=high_cols,levels=cloud_levs,alpha = 0.15,extend='max')#colors=['dimgray','gray','darkgray','slategrey','silver','lightgray'])
    cbar2 = fig.colorbar(tccp,orientation='vertical',pad=0.01,ax=ax2,aspect=50,extendrect=False, ticks=np.arange(10,100,10),shrink=0.9)
    cbar2.set_label('Cloud Cover (%)',fontsize=14)

    thec = ax2.contour(x,y,ndimage.gaussian_filter(thetae_2m,sigma=4,order=0),cmap='Blues_r',levels=range(270,370,4),extend='both')
    ax2.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], color='lightgray', length=6)
    ax2.quiver(x_2d[quiver_slices],y_2d[quiver_slices],u_storm[quiver_slices],v_storm[quiver_slices],color="white",**quiver_kwargs)
    ax2.quiver(x_2d[quiver_slices],y_2d[quiver_slices],u6km_shr[quiver_slices],v6km_shr[quiver_slices],color="yellow",**quiver_kwargs)


    blue = mpatches.Patch(color='#119579', label='High Clouds')
    green = mpatches.Patch(color='#25e935', label='Mid-Level Clouds')
    purple = mpatches.Patch(color='#ae14c2',label='Low-Level Clouds')
    leg = ax2.legend(handles=[blue,green,purple],loc=4,framealpha=1)
    leg.set_zorder(100)

    #### BOTTOM LEFT PANEL ########
    tprecip = ax3.contourf(x,y,total_precip, alpha = 0.7, cmap = 'cool',transform=zH5_crs,
                            levels=[0.01,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3,3.5,4,4.5,5])
    tcprecip = ax3.contour(x,y,total_precip,colors=['b','darkblue','darkviolet'],levels=[1,3,5],linewidths=1.5)
    pwc = ax3.contourf(x,y,pwat,cmap='Greens',levels=pwlevs,alpha=0.5)

    ax3.quiver(x_2d[quiver_slices],y_2d[quiver_slices],uivt[quiver_slices],vivt[quiver_slices],color="lime")
    cbar3 = fig.colorbar(pwc,orientation='vertical',pad=0.01,ax=ax3,aspect=50,extendrect=False,
                        ticks=pwlevs,shrink=0.9)
    cbar3.set_label('Precipitable Water (inches)',fontsize=14)

    #### TOP RIGHT PANEL ########
    mlccf = ax4.contourf(x,y,mllapse,colors=['darkblue','blue','royalblue','cornflowerblue'],levels=[-10,-8,-7,-6],alpha=0.1)
    refp = ax4.contourf(x,y,reflectivity, levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Greens',transform=zH5_crs) #colors=['#0099ff00', '#4D8080ff', '#666666ff', '#804d4dff','#993333ff','#B33333ff','#CC1a1aff','#E60000ff','#0000e6','#0000cc','#0000b3','#2d00b3','#5900b3','#8600b3','#b300b3','#b30086'])
    capep = ax4.contourf(x, y, cape, levels=[100,250,500,1000,1500,2000,2500,3000], extend='max', alpha = 0.7, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
    lgt = ax4.contour(x,y,lightning,levels=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
    mlcc = ax4.contour(x,y,mllapse,colors=['darkblue','blue','royalblue','cornflowerblue'],levels=[-10,-8,-7,-6],linewidths=1.5)

    cb = fig.colorbar(capep, orientation='vertical', pad = 0.01, aspect = 50, ax = ax4, extendrect=False, ticks=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000],shrink=0.9)
    cb.set_label('CAPE (J/kg)', size='large')
    ax4.barbs(x[wind_slice_zo],y[wind_slice_zo],u5m[wind_slice_zo,wind_slice_zo],v5m[wind_slice_zo,wind_slice_zo],length=6,color='blue')
    ax4.barbs(x[wind_slice_zo],y[wind_slice_zo],u8m[wind_slice_zo,wind_slice_zo],v8m[wind_slice_zo,wind_slice_zo],length=6,color='yellow')

    lblue = mpatches.Patch(color='cornflowerblue', label='6-7 C/km')
    blue = mpatches.Patch(color='royalblue', label='7-8 C/km')
    dblue = mpatches.Patch(color='blue', label='8-9 C/km')
    ddblue = mpatches.Patch(color='darkblue', label='>9 C/km')
    leg = ax4.legend(handles=[ddblue,dblue,blue,lblue],loc=4,title='H5-7 Lapse Rate',framealpha=1)
    leg.set_zorder(100)


    #ax4.barbs(x[ws2],y[ws2],s5u[ws2,ws2],s5v[ws2,ws2], length = 7)

    #### BOTTOM RIGHT PANEL ########
    cape3c = ax5.contourf(x,y,cape3,cmap='Reds',extend='max',levels=range(25,200,25))
    lclcc = ax5.contour(x,y,lclsmoothed,levels=[500,1000,1500,2000],
            colors=['darkgreen','limegreen','lightgreen','greenyellow','yellowgreen'],linewidths=1.0)
    lclcf = ax5.contourf(x,y,lclsmoothed,levels=[500,1000,1500,2000],
            colors=['darkgreen','limegreen','lightgreen','greenyellow','yellowgreen'],extend='both',alpha=0.15)
    cape3_100c = ax5.contour(x,y,cape3,color='4',levels=[100],linewidths=1.5)
    print(np.max(sfcvort))
    print(np.min(sfcvort))
    sfcvortc = ax5.contour(x,y,sfcvort,levels=range(1,50,1),linewidths=0.5,cmap='Blues_r')
    h_contour = ax5.contour(x, y, mslpc, colors='lightgray', levels=range(940,1040,2),linewidths=0.5)
    ax5.barbs(x[wind_slice_zo],y[wind_slice_zo],u1km_shr[wind_slice_zo,wind_slice_zo],v1km_shr[wind_slice_zo,wind_slice_zo], length=6, color='white')
    cbar5 = fig.colorbar(cape3c,orientation='vertical',pad=0.01,ax=ax5,aspect=50,extendrect=False,ticks=range(25,200,25),shrink=0.9)
    cbar5.set_label('0-3km CAPE (J/kg)',fontsize=14)

    darkgreen = mpatches.Patch(color='darkgreen', label='<0.5km')
    green = mpatches.Patch(color='limegreen', label='.5-1km')
    lightgreen = mpatches.Patch(color='lightgreen', label='1-1.5km')
    yellowgreen = mpatches.Patch(color='greenyellow', label='1.5-2km')
    yellowish = mpatches.Patch(color='yellowgreen', label='>2km')
    leg = ax5.legend(handles=[darkgreen,green,lightgreen,yellowgreen,yellowish],loc=4,title='LCL Height',framealpha=1)
    leg.set_zorder(100)

    ax2.set_title('Cloud Cover/Storm Mode Toolkit')
    ax3.set_title('Flash Flood Toolkit')
    ax4.set_title('General Severe Toolkit')
    ax5.set_title('Tornado Toolkit')

    ax3.set_facecolor('dimgray')
    ax2.set_facecolor('dimgray')
    ax4.set_facecolor('dimgray')
    ax5.set_facecolor('dimgray')
    sub_w1 = 260
    sub_w = 262
    sub_e = 295
    sub_n = 50
    sub_s = 25

    ax1.set_extent((sub_w1-1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((sub_w, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
    #fig.canvas.draw()
    fig.tight_layout()
    plt.savefig(output_dir+'/HRRR_ex/EC_fivepanel_summer_v16_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    if extradoms == True:
        for j in range(len(centerlon)):
            sub_w1 = 360-(float(centerlon[j])-7)
            sub_e1 = 360-(float(centerlon[j])+7)
            sub_n1 = float(centerlat[j])+5
            sub_s1 = float(centerlat[j])-5

            ax1.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=6)
            #ax1.quiver(x_2d[quiver_slices_reg],y_2d[quiver_slices_reg],u_storm[quiver_slices_reg],v_storm[quiver_slices_reg],color="white")

            ax1.set_extent((sub_w1, sub_e1, sub_s1, sub_n1))#, crs = zH5_crs)    # Set a title and show the plot
            ax2.set_extent((sub_w1, sub_e1, sub_s1, sub_n1))#, crs = zH5_crs)    # Set a title and show the plot
            ax3.set_extent((sub_w1, sub_e1, sub_s1, sub_n1))#, crs = zH5_crs)    # Set a title and show the plot
            ax4.set_extent((sub_w1, sub_e1, sub_s1, sub_n1))#, crs = zH5_crs)    # Set a title and show the plot
            ax5.set_extent((sub_w1, sub_e1, sub_s1, sub_n1))#, crs = zH5_crs)    # Set a title and show the plot
            #fig.canvas.draw()
            fig.tight_layout()
            plt.savefig(output_dir+'/HRRR_ex/'+state_abbr[j]+'_svr_dash_v22_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)

    ax1.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=6)
    ax1.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax2.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax3.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax4.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
    ax5.set_extent((283, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot

    for i in range(len(stations)):
        station = stations[i]
        lat = coords[i][0]
        lon = coords[i][1]
        ax3.text(360+lon,lat,str(station_qpf[i]),ha='center',transform=zH5_crs,color='white')
        ax1.text(360+lon,lat,str(station_temp[i]),ha='center',transform=zH5_crs)

    plt.savefig(output_dir+'/HRRR_ex/NE_fivepanel_summer_v16_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
    plt.clf()

    if plot_overview == True:
        fig11 = plt.figure(figsize=(15,15))
        ax11 = fig11.add_subplot(111,projection=zH5_crs)
        ax11.coastlines(resolution='10m')
        ax11.add_feature(cfeature.BORDERS.with_scale('10m'))
        ax11.add_feature(cfeature.STATES.with_scale('10m'))

        tmp_2m = ax11.contourf(x,y,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=zH5_crs)
        tmp_2m32 = ax11.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])


        h_contour = ax11.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
        h_contour.clabel(fontsize=14, colors='dimgray', inline=True, fmt='%i mb', rightside_up=True, use_clabeltext=True)

        ref_levs = range(5,75,1)

        norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
        newcmap = ListedColormap(cmap_ref(range(40, 194)))
        newnorm = Normalize(0,77)
        refc = ax11.contourf(x,y,reflectivity,norm=newnorm,cmap=newcmap,levels=range(5,75,1),alpha=0.7)
        cbars.addt2mcolorbar(ax11,fig,tmp_2m,range(-20,100,5))
        cbars.addrefcolorbar(ax11,fig,refc,ref_levs)

        tdc = ax11.contour(x,y,td2m,colors=['springgreen','lawngreen','limegreen','green'],levels=[55,60,65,70],linewidths=1.5)
        tdf = ax11.contourf(x,y,td2m,colors=['springgreen','lawngreen','limegreen','green'],levels=[55,60,65,70,100],alpha=0.15)

        blue = mpatches.Patch(color='green', label='>70F')
        orange = mpatches.Patch(color='limegreen', label='65-70F')
        pink = mpatches.Patch(color='lawngreen', label='60-65F')
        green = mpatches.Patch(color='springgreen', label='55-60F')
        leg = ax11.legend(handles=[blue,orange,pink,green],loc=4,title='2m Dew Point',framealpha=1)
        leg.set_zorder(100)
        ax11.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=6)
        ax11.set_title('Reflectivity (dBZ), 2m Temperature (F), 10m Winds (kts), and MSLP (mb)',fontsize=14)
        ax11.set_title('\n Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=11,loc='right')
        ax11.set_title('\n HRRR Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
        ax11.set_extent((sub_w1-1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
        #fig.canvas.draw()
        fig11.tight_layout()
        plt.savefig(output_dir+'/HRRR_ex/EC_overview_summer_v1_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax11.barbs(x[wind_slice_ne],y[wind_slice_ne],u_10m[wind_slice_ne,wind_slice_ne],v_10m[wind_slice_ne,wind_slice_ne], length=6)
        ax11.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
        for i in range(len(stations)):
            station = stations[i]
            lat = coords[i][0]
            lon = coords[i][1]
            ax11.text(360+lon,lat,str(station_temp[i]),ha='center',transform=zH5_crs)
        plt.savefig(output_dir+'/HRRR_ex/NE_overview_summer_v1_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        plt.clf()

    if plot_precip == True:
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
        plt.savefig(output_dir+'/HRRR_ex/EC_total_precip_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax7.set_extent((281, 295, 39, 49))
        plt.savefig(output_dir+'/HRRR_ex/NE_total_precip_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax7.set_extent((289,291,43,45))
        ax7.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
        plt.savefig(output_dir+'/HRRR_ex/local_total_precip_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        plt.close()
        plt.clf()

    if plot_winds == True:
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
        ax8.barbs(x[wind_slice],y[wind_slice],u_10m[wind_slice,wind_slice],v_10m[wind_slice,wind_slice], length=6)
        ax8.set_title('HRRR 10m Wind Forecast Valid at ' + time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=24)
        ax8.set_extent((sub_w-1, sub_e-1, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot

        h_contouqr = ax8.contour(x, y, mslpc, colors='dimgray', levels=range(940,1040,4),linewidths=2)
        h_contouqr.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

        plt.savefig(output_dir+'/HRRR_ex/EC_wind_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax8.set_extent((281, 295, 39, 49))
        plt.savefig(output_dir+'/HRRR_ex/NE_wind_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax8.set_extent((289,291,43,45))
        ax8.barbs(x[wsl],y[wsl],u_10m[wsl,wsl],v_10m[wsl,wsl], length=6)
        ax8.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
        plt.savefig(output_dir+'/HRRR_ex/local_wind_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        plt.close()
        plt.clf()
        plt.close()
        plt.clf()

    if plot_clouds ==True:
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
        plt.savefig(output_dir+'/HRRR_ex/EC_clouds_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax9.set_extent((281, 295, 39, 49))#, crs = zH5_crs)    # Set a title and show the plot
        plt.savefig(output_dir+'/HRRR_ex/NE_clouds_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        ax9.set_extent((289,291,43,45))
        ax9.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
        plt.savefig(output_dir+'/HRRR_ex/local_clouds_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
        plt.clf()
        plt.close()
        plt.clf()
        plt.close()

    if plot_vis == True:
        fig8 = plt.figure(figsize=(15,15))
        ax15 = fig8.add_subplot(111,projection=zH5_crs)

        ax15.coastlines(resolution='10m')
        ax15.add_feature(cfeature.BORDERS.with_scale('10m'))
        ax15.add_feature(cfeature.STATES.with_scale('10m'))

        visc1 = ax15.contourf(x,y,vis,levels=[0.00625,0.125,0.25,0.5,0.75,1,2,3,4,5],extend='min',colors=['#ffb3ff','magenta','deeppink','hotpink','mediumvioletred','crimson','orangered','darkorange','orange','gold'],alpha=0.7)
        cbr = fig.colorbar(visc1,orientation = 'horizontal', aspect = 80, ax = ax15, pad = 0.01,
                            extendrect=False, ticks = [0.00625,0.125,0.25,0.5,0.75,1,2,3,4,5], shrink=0.7)
        cbr.set_label('Surface Visibility (mi)')
        ax15.set_title('Surface Visibility')
        ax15.set_title('HRRR Init: '+init_time.dt.strftime('%m-%d %H:%MZ').item(),fontsize=14,loc='left')
        ax15.set_title('Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=14,loc='right')

        ax15.set_extent((260, 295, 25, 50))
        plt.savefig(output_dir+'/HRRR_ex/ec_vis_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
        ax15.set_extent((281,295,39,49))
        plt.savefig(output_dir+'/HRRR_ex/ne_vis_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
        ax15.set_extent((289,291,43,45))
        ax15.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')
        plt.savefig(output_dir+'/HRRR_ex/local_vis_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
        plt.close()
        plt.clf()

    print(str(i)+'_Done!')
