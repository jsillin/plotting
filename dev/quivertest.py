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
import matplotlib.colors as col
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize

#Data URL
url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr20210426/hrrr_sfc.t18z'

#Set vector spacing
wind_slice_zo = slice(60,-60,60)

#Set map bounds
sub_w1 = 260
sub_w = 262
sub_e = 295
sub_n = 50
sub_s = 25

#Parse data using MetPy
ds = xr.open_dataset(url)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

lats = np.arange(20,55,0.25)
lons = np.arange(240,300,0.25)

data = ds.metpy.assign_crs(grid_mapping_name='latitude_longitude')
data = data.isel(time=42)
print(data['hgt700mb'])

#Rename variables to useful things
data = data.rename({
        'ustm0_6000m':'u_storm',
        'vstm0_6000m':'v_storm',
})

zH5 = data['u_storm'].squeeze()
zH5_crs = zH5.metpy.cartopy_crs

time = data['u_storm'].metpy.time
x, y = data['u_storm'].metpy.coordinates('x', 'y')
lat, lon = xr.broadcast(y, x)
x_2d, y_2d = np.meshgrid(x, y)

dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

u_storm = data['u_storm'].squeeze()*1.94384449
v_storm = data['v_storm'].squeeze()*1.94384449


fig11 = plt.figure(figsize=(15,15))
ax11 = fig11.add_subplot(111,projection=zH5_crs)
ax11.coastlines(resolution='10m')
ax11.add_feature(cfeature.BORDERS.with_scale('10m'))
ax11.add_feature(cfeature.STATES.with_scale('10m'))

#Maybe this is how you're supposed to set vector spacing?
quiver_slices = (slice(None,None,50),slice(None,None,50))
ax11.quiver(x_2d[quiver_slices],y_2d[quiver_slices],u_storm[quiver_slices],v_storm[quiver_slices],color="lime")
ax11.set_title('Storm Motion',fontsize=14)
ax11.set_title('\n Valid: '+time.dt.strftime('%a %b %d %H:%MZ').item(),fontsize=11,loc='right')
ax11.set_title('\n HRRR Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
ax11.set_extent((sub_w1-1, sub_e, sub_s, sub_n))#, crs = zH5_crs)    # Set a title and show the plot
plt.savefig('NE_smv_v1_'+dtfs+'_.png',bbox_inches='tight',pad_inches=0.1)
