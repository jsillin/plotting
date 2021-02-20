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
from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.plots import SkewT
import matplotlib.lines as lines

def get_mdate(year,month,day,hour):
    if month <10:
        month = '0'+str(month)
    else:
        month = str(month)

    if day <10:
        day = '0'+str(day)
    else:
        day = str(day)

    if hour <10:
        hour = '0'+str(hour)
    else:
        hour = str(hour)

    mdate = str(year)+month+day
    init_hr = hour

    return [mdate,init_hr]

df = WyomingUpperAir.request_data(datetime(2021,2,16,12),'LCH')

obs_temp = df['temperature']
obs_pres = df['pressure']

print(obs_temp)

lat = df['latitude'][0]
lon = df['longitude'][0]

mdates = []
init_hrs = [0,6,12,18]
lats = np.arange(20,90,0.25)
lons = np.arange(190,350,0.25)

fig = plt.figure(figsize=(15, 15))
skew = SkewT(fig, rotation=45)
plt.title('Lake Charles Sounding Verification valid 12z 2-16-21')

lightred_line = lines.Line2D([], [], color='indianred', label='GFS Forecast Soundings')
red_line = lines.Line2D([], [], color='r',label='Observed Sounding')
leg = skew.ax.legend(handles=[lightred_line,red_line],title='Legend',loc=4,framealpha=1)
leg.set_zorder(100)

for i in range(7,16):
    for j in range(0,3):
        mdate = get_mdate(2021,2,i,init_hrs[j])[0]
        init_hr = get_mdate(2021,2,i,init_hrs[j])[1]
        print(mdate)
        print(init_hr)
        url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25/gfs'+mdate+'/gfs_0p25_'+init_hr+'z'
        ds = xr.open_dataset(url)
        data = ds.metpy.parse_cf()
        print(data.time)
        data = data.sel(time=datetime(2021,2,16,12))
        temp = data['tmpprs']
        #relh = data['rhprs']
        pres = data.lev
        sta_temp = temp.interp(lat = lat, lon = 360+lon)-273.15
        #sta_rh = relh.interp(lat = lat, lon = 360+lon)-273.15
        #sound_dp = mpcalc.dewpoint_from_relative_humidity(sta_temp.data*units.degC,sta_rh.data*units.percent)

        skew.plot(pres,sta_temp,color='indianred',linewidth=0.5)
        #skew.plot(pres,sta_rh,'midnightblue',linewidth=0.5)

        plt.savefig('lch_sounding_verif_'+str(i)+'_'+init_hr+'.png')

skew.plot(obs_pres,obs_temp,'r',linewidth=3)
plt.savefig('lch_final_verif.png')
