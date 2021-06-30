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
from metpy.plots import USCOUNTIES
import matplotlib.patches as mpatches
from siphon.simplewebservice.wyoming import WyomingUpperAir

#Fumctoion to make directory to store output
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

#Initialize some empty lists and select all RAOB stations of interest
coords = []
obs_z5 = []
ideal_ua_station_ids = ['UIL','OTX','TFX','GGW','BIS','ABR','MPX','RAP','RIW','BOI','MFR','SLE','INL','GRB','DVN','ILX','ILN','WMW','YMO','WPL','YQD','WSE','ZXS','CWVK',
                    'PANT','PAVA','YZT','OAK','REV','LKN','SLC','DNR','DDC','LFB','OAX','TOP','SGF','AMA','OUN','LZK','BNA','RNK','WAL','ALB','CHH','OKX','GYX','GSC','MHX',
                    'YSM','YYE','YAH','YZU','CAR','BUF','PIT','APX','DTX','YYQ','CHS','FFC','BMX','JAN','LIX','LCH','CRP','BRO','DRT','MAF','EPZ','ABQ','TUS','FGZ',
                    'VEF','NKX','VBG','FWD','YGI','PADQ','PACD','PANC','PAKN','PAFA','PABE','PAMC','PASN','PAOM','YEV','YVQ','YCB','YRB','YUX','YBK','YBF','YYR','YJT','AVT','PHLI','PHTO',
                    'TLH','MFL','XMR','JAX','TBW','EYW','BRO','AVT']
ua_station_ids = []

#Get current time and pick proper model run and verification time accordingly
current_utc = datetime.utcnow()
print(current_utc)

curr_utc_yr = current_utc.year
curr_utc_mo = current_utc.month
curr_utc_dy = current_utc.day

def get_verif_time(current_utc):
    if current_utc.hour <1:
        prev_utc = current_utc-dt.timedelta(days=1)
        verif_yr = prev_utc.year
        verif_mo = prev_utc.month
        verif_dy = prev_utc.day
        verif_time = datetime(verif_yr,verif_mo,verif_dy,12)
    elif current_utc.hour<13:
        verif_time = datetime(curr_utc_yr,curr_utc_mo,curr_utc_dy,0)
    else:
        verif_time = datetime(curr_utc_yr,curr_utc_mo,curr_utc_dy,12)
    return verif_time

def model_init_time(current_utc):
    if current_utc.hour <5:
        prev_utc = current_utc-dt.timedelta(days=1)
        init_yr = prev_utc.year
        init_mo = prev_utc.month
        init_dy = prev_utc.day
        init_time = datetime(init_yr,init_mo,init_dy,12)
    elif current_utc.hour<17:
        init_yr = current_utc.year
        init_mo = current_utc.month
        init_dy = current_utc.day
        init_time = datetime(init_yr,init_mo,init_dy,0)
    else:
        init_time = datetime(curr_utc_yr,curr_utc_mo,curr_utc_dy,12)
    return init_time

def model_init_string(init_time):
    year = init_time.year
    if init_time.month <10:
        month = '0'+str(init_time.month)
    else:
        month = str(init_time.month)

    if init_time.day <10:
        day = '0'+str(init_time.day)
    else:
        day = str(init_time.day)

    if init_time.hour <10:
        hour = '0'+str(init_time.hour)
    else:
        hour = str(init_time.hour)
    minute = '00'

    mdate = str(year)+str(month)+str(day)+'_'+hour+minute
    return mdate

verification_time = get_verif_time(current_utc)
model_init_time = model_init_time(current_utc)
mdate = model_init_string(model_init_time)
print(verification_time)
print(model_init_time)
print(mdate)

#Now that we have the model and verification time, get the RAOB data
for i in range(len(ideal_ua_station_ids)):
    try:
        df = WyomingUpperAir.request_data(verification_time,ideal_ua_station_ids[i])
        lat = df['latitude'][0]
        lon = df['longitude'][0]
        p = df['pressure']
        h5 = df[(df['pressure']==500)]
        z5 = h5['height'].values[0]
        obs_z5.append(z5)
        coord = (lat,lon)
        coords.append(coord)
        print(coord)
        ua_station_ids.append(ideal_ua_station_ids[i])
    except:
        'Data Unavailable For '+ideal_ua_station_ids[i]

#and the model data
url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs'+mdate[0:8]+'/gfs_0p25_1hr_'+mdate[9:11]+'z'

# Create new directory
output_dir = str(mdate)
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFS'
)
#Parse data using MetPy
ds = xr.open_dataset(url)

times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(15,70,0.25)
lons = np.arange(220,330,0.25)

data = ds.metpy.parse_cf()
data = data.sel(time=verification_time)

gph = data['hgtprs'].sel(lev=500).squeeze()
zH5_crs = data['hgtprs'].metpy.cartopy_crs
x, y = gph.metpy.coordinates('x', 'y')

predicted_z5 = []

#Grab predicted heights for each point
for i in range(len(ua_station_ids)):
    station_lat = coords[i][0]
    station_lon = coords[i][1]
    pz5 = gph.interp(lat=station_lat,lon=360+station_lon).values
    predicted_z5.append(round(float(pz5),1))

#Make a figure
fig = plt.figure(figsize=(15,15))
ax1 = fig.add_subplot(111, projection = zH5_crs)

ax1.coastlines(resolution='50m',edgecolor='gray',linewidths=0.5)
ax1.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray',linewidth=0.5)
ax1.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray',linewidth=0.5)

h5c = ax1.contour(x,y,gph,colors='dimgray', levels = range(4800,6200,60),linewidths=.75)
#Compute and plot the model error
for i in range(len(ua_station_ids)):
    diff = obs_z5[i]-predicted_z5[i]
    station = ua_station_ids[i]
    lat = coords[i][0]
    lon = coords[i][1]
    if diff <-10:
        ax1.plot(360+lon,lat,color='blue',linewidth=2,marker='o',transform=zH5_crs)
    elif diff >10:
        ax1.plot(360+lon,lat,color='red',linewidth=2,marker='o',transform=zH5_crs)
    else:
        ax1.plot(360+lon,lat,color='black',linewidth=2,marker='o',transform=zH5_crs)
    ax1.text(360+coords[i][1],coords[i][0]+.75,ua_station_ids[i],ha='center',transform=zH5_crs)
    ax1.text(360+coords[i][1],coords[i][0]-1.25,str(round(diff,1)),ha='center',transform=zH5_crs)

blue = mpatches.Patch(color='blue', label='Observed Heights >10m Lower Than GFS Forecast')
black = mpatches.Patch(color='black', label='Observed Heights Within 10m Of GFS Forecast')
red = mpatches.Patch(color='red', label='Observed Heights >10m Higher Than GFS Forecast')
leg = ax1.legend(handles=[red,black,blue],loc=4,title='500mb Height Verification',framealpha=1)
leg.set_zorder(100)


ax1.set_title(str(mdate)+' GFS vs '+str(verification_time)[0:16]+' RAOBs')
ax1.set_extent((360-135,360-65,20,60))
plt.savefig(output_dir+'/GFS_H5_CONUS_Verif_'+mdate+'.png')
ax1.set_extent((360-175,360-140,35,70))
plt.savefig(output_dir+'/GFS_H5_AK_Verif_'+mdate+'.png')
ax1.set_extent((360-174,360-57,10,70))
plt.savefig(output_dir+'/GFS_H5_NAM_Verif_'+mdate+'.png')

###########################################################

#Same as fig1 but with a slightly different color table
fig2 = plt.figure(figsize=(15,15))
ax2 = fig2.add_subplot(111, projection = zH5_crs)

ax2.coastlines(resolution='50m',edgecolor='gray',linewidths=0.5)
ax2.add_feature(cfeature.BORDERS.with_scale('50m'),edgecolor='gray',linewidth=0.5)
ax2.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='gray',linewidth=0.5)

h5c = ax2.contour(x,y,gph,colors='dimgray', levels = range(4800,6200,60),linewidths=.75)
for i in range(len(ua_station_ids)):
    diff = obs_z5[i]-predicted_z5[i]
    station = ua_station_ids[i]
    lat = coords[i][0]
    lon = coords[i][1]
    if diff <-15:
        ax2.plot(360+lon,lat,color='blue',linewidth=2,marker='o',transform=zH5_crs)
    elif (diff <-10)&(diff>-15):
        ax2.plot(360+lon,lat,color='lightskyblue',linewidth=2,marker='o',transform=zH5_crs)
    elif (diff >10)&(diff<15):
        ax2.plot(360+lon,lat,color='orange',linewidth=2,marker='o',transform=zH5_crs)
    elif diff >15:
        ax2.plot(360+lon,lat,color='red',linewidth=2,marker='o',transform=zH5_crs)
    else:
        ax2.plot(360+lon,lat,color='black',linewidth=2,marker='o',transform=zH5_crs)
    ax2.text(360+coords[i][1],coords[i][0]+.75,ua_station_ids[i],ha='center',transform=zH5_crs)
    ax2.text(360+coords[i][1],coords[i][0]-1.25,str(round(diff,1)),ha='center',transform=zH5_crs)

blue = mpatches.Patch(color='blue', label='Observed Heights >15m Lower Than GFS Forecast')
lblue = mpatches.Patch(color='lightskyblue', label='Observed Heights 10-15m Lower Than GFS Forecast')
black = mpatches.Patch(color='black', label='Observed Heights Within 10m Of GFS Forecast')
orange = mpatches.Patch(color='orange', label='Observed Heights 10-15m Higher Than GFS Forecast')
red = mpatches.Patch(color='red', label='Observed Heights >15m Higher Than GFS Forecast')
leg = ax2.legend(handles=[red,orange,black,lblue,blue],loc=4,title='500mb Height Verification',framealpha=1)
leg.set_zorder(100)


ax2.set_title(str(mdate)+' GFS vs '+str(verification_time)[0:16]+' RAOBs')
ax2.set_extent((360-135,360-65,18,60))
plt.savefig(output_dir+'/GFS_H5_CONUS_Verif__v1'+mdate+'.png')
ax2.set_extent((360-175,360-140,35,70))
plt.savefig(output_dir+'/GFS_H5_AK_Verif__v1'+mdate+'.png')
ax2.set_extent((360-176,360-55,8,70))
plt.savefig(output_dir+'/GFS_H5_NAM_Verif__v1'+mdate+'.png')
