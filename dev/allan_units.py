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
from metpy.plots import SkewT
import metpy.calc as mpcalc
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import soundingmaps as smap
import supplementary_tools as spt
from metpy.calc import wind_components, lcl, lfc, el, dry_lapse, parcel_profile, dewpoint, wet_bulb_temperature, mean_pressure_weighted
from metpy.calc import most_unstable_parcel, parcel_profile_with_lcl, bulk_shear, storm_relative_helicity
from metpy.calc import wind_speed, wind_direction, thermo, vapor_pressure, bunkers_storm_motion, pressure_to_height_std
from metpy.plots import SkewT, Hodograph
from metpy.units import units, concatenate
import sharppy.sharptab.profile as profile
from sharppy.sharptab import utils, winds, params, interp, thermo, watch_type, fire
import sharppy.sharptab as tab
from sharppy.sharptab.constants import *
from sharppy.sharptab.constants import MISSING, TOL

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
#grabbing data from NOMADS
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

# Create new directory to store output
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00' #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFS/STL') #create subdirectory to store GFS output like this

#This actually opens the dataset from NOMADS and parses it with MetPy
ds = xr.open_dataset(url)
init_hr = dt.datetime(year,int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time #Pull out the time dimension
init_time = ds['time'][0]

#Subset the data to only work with certain lats and lons of interest
lats = np.arange(20,55,0.25)
lons = np.arange(240,310,0.25)

ds = ds.sel(lat = lats, lon = lons)

#Now loop through the 36 forecast hours to make the plots
for i in range(0,84):
    #Get the data for the forecast hour of interest
    data = ds.metpy.assign_crs(grid_mapping_name='latitude_longitude')
    data = data.isel(time=i)

    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tmpprs': 'temperature',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'rhprs':'rh',
        'capesfc':'cape',
        'ugrd10m':'u',
        'vgrd10m':'v',
        'hlcy3000_0m':'srh',
        'ugrdprs':'uwnd',
        'vgrdprs':'vwnd',
        'hlcy3000_0m':'heli',
        'pressfc':'spres',
        'hgtsfc':'surface_hgt',
        'hgtprs':'z',
        'hgt0c':'hgt0c'
    })

    cape = data['cape'].squeeze()
    u_10m = data['u'].squeeze()
    v_10m = data['v'].squeeze()
    u_10m = u_10m*1.94384449
    v_10m = v_10m*1.94384449
    wspd = ((u_10m**2)+(v_10m**2))**.5

    #Coordinate stuff
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lats, lons = xr.broadcast(y, x)
    zH5_crs = data['temperature'].metpy.cartopy_crs
    wind_slice = slice(-1,1,1)
    #Processing surface temperature data
    t2m = data['sfc_temp'].squeeze()
    t2m = ((t2m - 273.15)*(9./5.))+32.

    #This creates a nice-looking datetime label
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    td2m = data['sfc_td'].squeeze()
    td2m = ((td2m - 273.15)*(9./5.))+32.
    td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)
    sfc_pressure = data['spres'].squeeze()/100
    sound_pres = data.lev
    ptop=100
    soundlat = 38.6270
    soundlon = 90.1994

    soundlat = soundlat
    soundlon = 360-(soundlon)
    sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

    sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
    sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
    #sound_wb = mpcalc.wet_bulb_temperature(sound_pres,sound_temps.data*units.degC,sound_dp)
    u_wind_s = data['uwnd'].interp(lat=soundlat,lon=soundlon)
    v_wind_s = data['vwnd'].interp(lat=soundlat,lon=soundlon)
    u7 = u_wind_s*units.knots
    v7 = v_wind_s*units.knots
    spd = ((u_wind_s**2)+(v_wind_s**2))**.5
    w_dir = wind_direction(u7.data, v7.data) * units.deg
    u,v = u_wind_s*1.94384,v_wind_s*1.94384 # m/s to knots
    spres = sfc_pressure.interp(lat=soundlat,lon=soundlon)
    surface_hgt = data['surface_hgt'].interp(lat=soundlat,lon=soundlon)
    z = data['z'].interp(lat=soundlat,lon=soundlon)
    z7 = z*units.meter
    hgt0c = data['hgt0c'].interp(lat=soundlat,lon=soundlon)

    plt.rcParams['figure.figsize'] = (12, 14)
    fig = plt.figure(figsize=(24, 14))

    gs = fig.add_gridspec(ncols=2,nrows=1)

	# identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
    #ax3 = fig.add_subplot(gs[1, 0])
    #ax4 = fig.add_subplot(gs[1, 1])

    # Grid for plots
    #skew = SkewT(fig, rotation=45, subplot=gs[0, 0])

    # Grid for plots
    skew = SkewT(fig, rotation=45, subplot=gs[0, 0])


	# Plot the data
    skew.plot(sound_pres, sound_temps, 'r', linewidth=2)
    #skew.plot(sound_pres, sound_dp, 'b', linewidth=2)
    #skew.plot(p, td2, 'y')
    skew.plot(sound_pres, sound_dp, 'g', linewidth=2)
    skew.plot_moist_adiabats(color='grey',alpha=0.2)
    skew.plot_mixing_lines(color='grey',alpha=0.2)
    skew.plot_dry_adiabats(color='grey',alpha=0.2)
    #skew.plot(p, wetbulb, 'b', linewidth=1)

    #Only want data above the ground
    abv_sfc_temp = spt.mask_below_terrain(spres,sound_temps,sound_pres)[0]
    abv_sfc_dewp = spt.mask_below_terrain(spres,sound_dp,sound_pres)[0]
    #abv_sfc_wetb = spt.mask_below_terrain(spres,sound_wb,sound_pres)[0]
    pres_abv_ground = spt.mask_below_terrain(spres,sound_temps,sound_pres)[1]
    abv_sfc_u = spt.mask_below_terrain(spres,u_wind_s,sound_pres)[0]
    abv_sfc_v = spt.mask_below_terrain(spres,v_wind_s,sound_pres)[0]

    # Calculate LCL height and plot as black dot
    lcl_pressure, lcl_temperature = lcl(sound_pres, sound_temps.data*units.degC, sound_dp[0])

    # Calculate LFC height and plot as yellow dash
    lfc_pressure, lfc_temperature = lfc(sound_pres, sound_temps.data*units.degC, sound_dp)

    el_pressure, el_temperature = el(sound_pres, sound_temps.data*units.degC, sound_dp)

    sb_cape, sb_cin = mpcalc.surface_based_cape_cin(sound_pres, sound_temps.data*units.degC, sound_dp)
    ml_cape, ml_cin = mpcalc.mixed_layer_cape_cin(sound_pres, sound_temps.data*units.degC, sound_dp)
    mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(sound_pres, sound_temps.data*units.degC, sound_dp)

    mixed_0_3 = mpcalc.mixed_parcel(sound_pres, sound_temps.data*units.degC, sound_dp, depth=3000 * units.meter)
    print(mixed_0_3)

    sbcape = np.round(sb_cape, 1)
    sbcin = np.round(sb_cin, 1)
    mlcape = np.round(ml_cape, 1)
    mlcin = np.round(ml_cin, 1)
    mucape = np.round(mu_cape, 1)
    pw = mpcalc.precipitable_water(sound_pres,sound_dp)
    pw = pw.to(units.inch)
    pw = round(pw, 2)

    lcl_hgt = np.round(mpcalc.pressure_to_height_std(lcl_pressure), decimals=1).to(units.meter)/1000
    lfc_hgt = np.round(mpcalc.pressure_to_height_std(lfc_pressure), decimals=1).to(units.meter)/1000
    el_hgt = np.round(mpcalc.pressure_to_height_std(el_pressure), decimals=1).to(units.meter)/1000

    parcel_prof = mpcalc.parcel_profile(pres_abv_ground,abv_sfc_temp[0].data*units.degC,abv_sfc_dewp[0])
    abv_sfc_parcel = spt.mask_below_terrain(spres,parcel_prof,sound_pres)[0]
    cape = mpcalc.cape_cin(pres_abv_ground,abv_sfc_temp.data*units.degC,abv_sfc_dewp,parcel_prof)
    capeout = int(cape[0].m)
    cinout = int(cape[1].m)

    #prof = profile.create_profile(profile='default', pres=sound_pres[::-1], hgt=z[::-1], tmpc=sound_temps[::-1], dwpc=sound_dp[::-1], wspd=spd[::-1], wdir=w_dir[::-1], missing=-9999, strictQC=False)

    u_shear01, v_shear01 = mpcalc.bulk_shear(sound_pres, u7.data, v7.data, depth = 1000 * units.meter)
    shear01 = np.round((np.sqrt(u_shear01**2 + v_shear01**2)), 1)
    shear01 = shear01.to(units.knots)
    shear01 = np.round(shear01)
    u_shear005, v_shear005 = mpcalc.bulk_shear(sound_pres, u7.data, v7.data, depth = 500 * units.meter)
    shear005 = np.round((np.sqrt(u_shear005**2 + v_shear005**2)),1)
    shear005 = shear005.to(units.knots)
    shear005 = np.round(shear005)
    u_shear015, v_shear015 = mpcalc.bulk_shear(sound_pres, u7.data, v7.data, depth = 1500 * units.meter)
    shear015 = np.round((np.sqrt(u_shear015**2 + v_shear015**2)),1)
    shear015 = shear015.to(units.knots)
    shear015 = np.round(shear015)
    u_shear02, v_shear02 = mpcalc.bulk_shear(sound_pres, u7.data, v7.data, depth = 2000 * units.meter)
    shear02 = np.round((np.sqrt(u_shear02**2 + v_shear02**2)), 1)
    shear02 = shear02.to(units.knots)
    shear02 = np.round(shear02)
    u_shear03, v_shear03 = mpcalc.bulk_shear(sound_pres, u7.data, v7.data, depth = 3000 * units.meter)
    shear03 = np.round((np.sqrt(u_shear03**2 + v_shear03**2)), 1)
    shear03 = shear03.to(units.knots)
    shear03 = np.round(shear03)
    u_shear06, v_shear06 = mpcalc.bulk_shear(sound_pres,  u7.data, v7.data, depth = 6000 * units.meter)
    shear06 = np.round((np.sqrt(u_shear06**2 + v_shear06**2)), 1)
    shear06 = shear06.to(units.knots)
    shear06 = np.round(shear06)
    rmover, lmover, mean = mpcalc.bunkers_storm_motion(sound_pres, u7.data, v7.data, z7.data)

    right_mover,left_mover,wind_mean = mpcalc.bunkers_storm_motion(sound_pres, u7.data, v7.data, z7.data)
    print(wind_mean)
    wind_mean = np.round(wind_mean)
    wind_mean = wind_mean.to(units.kt)
    #wind_mean = abs(wind_mean)
    wind_mean = np.round(wind_mean)

    srh_01_pos, srh_01_neg, srh_01_tot = mpcalc.storm_relative_helicity( z7.data, u7.data, v7.data,
                                                                        depth = 1000 * units.meter, bottom = z7.data[0], storm_u = lmover[0], storm_v = lmover[1])
    srh_01 = np.round(srh_01_neg, 1)
    srh_03_pos, srh_03_neg, srh_03_tot = mpcalc.storm_relative_helicity( z * units('m'), u * units('m/s'), v * units('m/s'),
                                                                        depth = 3000 * units.meter, bottom = z[0]*units.meter, storm_u = lmover[0], storm_v = lmover[1])
    srh_03 = np.round(srh_03_neg, 1)

    # Need to round these numbers to make the string look pretty
    tot_SRH = np.round(total_SRH)
    tot1_SRH = np.round(total1_SRH)
    tot05_SRH = np.round(total05_SRH)
    tot13_SRH = np.round(total13_SRH)
    tot2_SRH = np.round(total2_SRH)

    skew.plot(pres_abv_ground, abv_sfc_parcel, 'k', linewidth=2)
    skew.plot_barbs(pres_abv_ground, abv_sfc_u, abv_sfc_v)
    skew.shade_cin(pres_abv_ground, abv_sfc_temp.data*units.degC, parcel_prof, alpha=0.2)
    skew.shade_cape(pres_abv_ground, abv_sfc_temp.data*units.degC, parcel_prof)
    skew.plot(lcl_pressure[0], lcl_temperature[0], marker="_", color='orange', markersize=30, markeredgewidth=3, label='LCL')
    skew.plot(lfc_pressure, lfc_temperature, marker="_", color='brown', markersize=30, markeredgewidth=2, label='LFC')
    skew.plot(el_pressure, el_temperature, marker="_", color='darkblue', markersize=30, markeredgewidth=2, label='EL')
    skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
    skew.ax.set_ylim((1000,ptop))
    plt.xlabel("Temperature (C)")
    plt.ylabel("Pressure (hPa)")

    # Draw hodograph
    ax_hod = fig.add_subplot(gs[0, 1])
    h = Hodograph(ax_hod, component_range=spd.max())
    h.add_grid(increment=20, alpha=0.2)
    h.plot_colormapped(u_wind_s, v_wind_s, spd)
    origin = np.array([[0, 0, 0],[0, 0, 0]])
    #plt.quiver(*origin, wind_mean, color='grey', scale=21)
    #plt.quiver(*origin,  bunk_left_dir, bunk_left_spd, color='grey', scale=21)
    #plt.quiver(*origin, u_storm, v_storm, color='grey', scale=21)

     ######## Save the plot
    plt.savefig(output_dir+'/GFS/STL/GFS_hrly_pytpe_stl_sound_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.close()
    plt.clf()
