#from mpl_toolkits.basemap import Basemap, cm
import os.path
import sys
from matplotlib import rcParams
from matplotlib.animation import ArtistAnimation
import matplotlib
import matplotlib.pyplot as plt
#import pyart
from siphon.radarserver import RadarServer
from datetime import datetime, timedelta
from siphon.cdmr import Dataset
#import pyart
import numpy as np
import numpy.ma as ma
import netCDF4
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy.ma as ma
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from metpy.plots import StationPlot
from metpy.plots.wx_symbols import sky_cover
from metpy.calc import (bunkers_storm_motion, bulk_shear, dewpoint, dewpoint_from_relative_humidity, dry_lapse, moist_lapse, vapor_pressure, saturation_vapor_pressure,
                        wind_speed, wind_direction, pressure_to_height_std, mixing_ratio, cape_cin, wind_components,
                        height_to_pressure_std, equivalent_potential_temperature, parcel_profile, precipitable_water,
                        storm_relative_helicity, mean_pressure_weighted, most_unstable_cape_cin, most_unstable_parcel,
                        supercell_composite, significant_tornado, get_layer, relative_humidity_from_dewpoint, mixing_ratio_from_relative_humidity, specific_humidity_from_dewpoint,
                        mixing_ratio_from_specific_humidity)
#from metpy.calc.tools import log_interp, get_layer, get_layer_heights
from metpy.calc import wind_direction
from metpy.units import units
from metpy.calc import lcl
from metpy.interpolate import interpolate_1d as metinterp
import scipy.ndimage as ndimage
import matplotlib.gridspec as gridspec
from metpy.calc import wind_components, lcl, lfc, el, dry_lapse, parcel_profile, dewpoint, wet_bulb_temperature, mean_pressure_weighted
from metpy.calc import most_unstable_parcel, parcel_profile_with_lcl, bulk_shear, storm_relative_helicity, lifted_index
from metpy.calc import wind_speed, wind_direction, thermo, vapor_pressure, bunkers_storm_motion, pressure_to_height_std
from metpy.plots import SkewT, Hodograph
import metpy.calc as metcalc
from datetime import datetime
import datetime as dt
import sharppy
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo

import pickle
import csv

from siphon.catalog import TDSCatalog
import xarray as xr
from xarray.backends import NetCDF4DataStore
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def total_shear(u, v, heights, depth, bottom=0 * units.m,
                            storm_u=0 * units('m/s'), storm_v=0 * units('m/s')):
    r"""Calculate total shear. Will output both length of the hodograph over the layer and normalized shear.

    Parameters
    ----------
    u : array-like
        u component winds
    v : array-like
        v component winds
    heights : array-like
        atmospheric heights, will be converted to AGL
    depth : number
        depth of the layer
    bottom : number
        height of layer bottom AGL (default is surface)
    storm_u : number
        u component of storm motion (default is 0 m/s)
    storm_v : number
        v component of storm motion (default is 0 m/s)
    Returns
    -------
    `pint.Quantity`, `pint.Quantity`
        hodo_length, tot_shear
    """
    _, u, v, z = get_layer_heights(heights, depth, u, v, heights, with_agl=True, bottom=bottom)

    storm_relative_u = u - storm_u
    storm_relative_v = v - storm_v

    #Get storm-relative wind vectors for the middle of each layer,
    #as well as shear vectors for each layer
    sr_ushear = (storm_relative_u[1:] - storm_relative_u[:-1])
    sr_vshear = (storm_relative_v[1:] - storm_relative_v[:-1])

    layer_depths = (z[1:] - z[:-1])

    #Get magnitudes of the shear vectors, sum to get hodo length
    shear_mags = np.sqrt(sr_ushear**2 + sr_vshear**2)
    shear_mags_norm = np.sqrt(sr_ushear**2 + sr_vshear**2)/layer_depths
    hodo_length = np.sum(shear_mags)

    #Normalize by depth
    tot_shear = np.sum(shear_mags_norm)

    return (hodo_length*units('m/s'), tot_shear)

cenlat = 38.6270
cenlon = -90.1994
#########################################
#grabbing data from NOMADS
startTime=datetime.now()
year,month,day,hour = startTime.year,startTime.month,startTime.day,startTime.hour
time_start = datetime(year, month, day, hour, 0) # Our specified time
hour = time_start.hour
if hour < 10:
    hour = '0'+str(hour)
day = time_start.day
if day < 10:
    day = '0'+str(day)
month = time_start.month
if month < 10:
    month = '0'+str(month)


#Importing relevant libraries
import matplotlib.lines as lines
import matplotlib.patches as mpatches
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
import scipy.ndimage as ndimage

# Any import of metpy will activate the accessors
import metpy.calc as mpcalc
#from metpy.testing import get_test_data
from metpy.units import units

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
'''
# Create new directory
output_dir = str(runtime)
mkdir_p(output_dir)
mkdir_p(output_dir+'/synoptic')
'''

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

mdate = str(year)+str(month)+str(day)+str(hour)
# Create new directory to store output
output_dir = str(year)+str(month)+str(day)+str(hour)  #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/Soundings/GFS/STL') #create subdirectory to store GFS output like this

#Get data using siphon
for i in range(0,45):
    #Get data using siphon
    best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
    best_ds = best_gfs.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=45, south=35, east=-83, west=-100).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Geopotential_height_isobaric').variables('Convective_available_potential_energy_surface')

    data = ncss.get_data(query)

    #Parse data using MetPy
    ds = xr.open_dataset(NetCDF4DataStore(data))
    data = ds.metpy.parse_cf()
    filtered_ds = data.filter_by_attrs(standard_name='forecast_reference_time').squeeze()
    coord_names = list(filtered_ds.coords)
    runtime = filtered_ds[coord_names[0]].dt.strftime('%Y%m%d_%H%M').values


    query = ncss.query()
    query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Vertical_velocity_pressure_isobaric').variables('Convective_available_potential_energy_surface').variables('u-component_of_wind_isobaric').variables('v-component_of_wind_isobaric').variables('Storm_relative_helicity_height_above_ground_layer').variables('Pressure_surface').variables('Dewpoint_temperature_height_above_ground').variables('Temperature_height_above_ground').variables('Vertical_u-component_shear_height_above_ground_layer').variables('Vertical_v-component_shear_height_above_ground_layer').variables('Geopotential_height_isobaric').variables('Geopotential_height_surface').variables('u-component_of_wind_height_above_ground').variables('v-component_of_wind_height_above_ground').variables('Relative_humidity_isobaric').variables('Temperature_isobaric').variables('MSLP_MAPS_System_Reduction_msl')
    query.add_lonlat().lonlat_box(cenlon-2.1, cenlon +2.1, cenlat-2.1, cenlat+2.1)
    data1 = ncss.get_data(query)
    '''
    data1 = data.rename({
        'Convective_available_potential_energy_surface': 'cape'
    })
    '''
    dtime = data1.variables['Geopotential_height_isobaric'].dimensions[0]
    time = data['Geopotential_height_isobaric'].metpy.time
    dlev = data1.variables['Geopotential_height_isobaric'].dimensions[1]
    dlat = data1.variables['Geopotential_height_isobaric'].dimensions[2]
    dlon = data1.variables['Geopotential_height_isobaric'].dimensions[3]
    x, y = data['Geopotential_height_isobaric'].metpy.coordinates('x', 'y')
    CAPE = np.asarray(data1.variables['Convective_available_potential_energy_surface'][:]) * units('J/kg')
    SRH = np.asarray(data1.variables['Storm_relative_helicity_height_above_ground_layer'][:]) * units('m/s')
    SFCP = (np.asarray(data1.variables['Pressure_surface'][:])/100.) * units('hPa')
    Td = (np.asarray(data1.variables['Dewpoint_temperature_height_above_ground'][:]) * units('kelvin')).to('degC')
    T = np.asarray(data1.variables['Temperature_height_above_ground'][:]) * units('kelvin')
    #ushr = data1.variables['Vertical_u-component_shear_height_above_ground_layer'][:] * units('m/s')
    #vshr = data1.variables['Vertical_v-component_shear_height_above_ground_layer'][:] * units('m/s')
    hgt = np.asarray(data1.variables['Geopotential_height_isobaric'][:]) * units('meter')
    sfc_hgt = np.asarray(data1.variables['Geopotential_height_surface'][:]) * units('meter')
    uwnd = np.asarray(data1.variables['u-component_of_wind_isobaric'][:]) * units('m/s')
    vwnd = np.asarray(data1.variables['v-component_of_wind_isobaric'][:]) * units('m/s')
    Temp_up = np.asarray(data1.variables['Temperature_isobaric'][:]) * units('kelvin')
    VVEL_up = np.asarray(data1.variables['Vertical_velocity_pressure_isobaric'][:]) * units('Pa/s')
    RH_up = np.asarray(data1.variables['Relative_humidity_isobaric'][:])
    usfc = np.asarray(data1.variables['u-component_of_wind_height_above_ground'][:]) * units('m/s')
    vsfc = np.asarray(data1.variables['v-component_of_wind_height_above_ground'][:]) * units('m/s')
    #MSLP = (np.asarray(data1.variables['Pressure_reduced_to_MSL_msl'][:])/100.) * units('hPa')
    # Get the dimension data
    lats_r = data1.variables[dlat][:]
    lons_r= data1.variables[dlon][:]
    lev = (np.asarray(data1.variables[dlev][:])/100.) * units('hPa')

    flon = float(cenlon)
    flat = float(cenlat)
    # Set up our array of latitude and longitude values and transform to
    # the desired projection.
    crs = ccrs.PlateCarree()
    crlons, crlats = np.meshgrid(lons_r[:]*1000, lats_r[:]*1000)
    trlatlons = crs.transform_points(ccrs.LambertConformal(central_longitude=270, central_latitude=38, standard_parallels=(25.,25.)),crlons,crlats)
    trlons = trlatlons[:,:,0]
    trlats = trlatlons[:,:,1]
    dlon = np.abs(trlons - cenlon)
    dlat = np.abs(trlats - cenlat)
    ilon = np.where(dlon == np.min(dlon))
    ilat = np.where(dlat == np.min(dlat))

    RH_prof = RH_up[0,:,ilat[0][0], ilon[1][0]]
    Tdc_up = dewpoint_from_relative_humidity(Temp_up[0,:,ilat[0][0], ilon[1][0]],RH_up[0,:,ilat[0][0], ilon[1][0]]/100)
    Omega = VVEL_up[0,:,ilat[0][0], ilon[1][0]]

    p_sounding = np.sort(np.append(lev, SFCP[0,ilat[0][0], ilon[1][0]]))
    ind = np.where(p_sounding >= SFCP[0,ilat[0][0], ilon[1][0]])[0][0]
    hgt_sounding = np.insert(hgt[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, sfc_hgt[0,ilat[0][0], ilon[1][0]].magnitude) * hgt.units
    T_sounding = (np.insert(Temp_up[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, T[0,0,ilat[0][0], ilon[1][0]].magnitude) * T.units).to(Tdc_up.units)
    Td_sounding = np.insert(Tdc_up.magnitude, ind, Td[0,0,ilat[0][0], ilon[1][0]].magnitude) * Tdc_up.units
    u_sounding = np.insert(uwnd[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, usfc[0,0,ilat[0][0], ilon[1][0]].magnitude) * usfc.units
    v_sounding = np.insert(vwnd[0,:,ilat[0][0], ilon[1][0]].magnitude, ind, vsfc[0,0,ilat[0][0], ilon[1][0]].magnitude) * usfc.units

    p_skewt = p_sounding[p_sounding <= SFCP[0,ilat[0][0], ilon[1][0]]]
    hgt_skewt = hgt_sounding[p_sounding <= SFCP[0,ilat[0][0], ilon[1][0]]]
    T_skewt = T_sounding[p_sounding <= SFCP[0,ilat[0][0], ilon[1][0]]]
    Td_skewt = Td_sounding[p_sounding <= SFCP[0,ilat[0][0], ilon[1][0]]]
    u_skewt = u_sounding[p_sounding <= SFCP[0,ilat[0][0], ilon[1][0]]].to('kt')
    v_skewt = v_sounding[p_sounding <= SFCP[0,ilat[0][0], ilon[1][0]]].to('kt')
    storm_motion = bunkers_storm_motion(p_skewt[::-1], u_skewt[::-1], v_skewt[::-1], hgt_skewt[::-1])

    #Get RH for whole sounding
    RH_skewt = relative_humidity_from_dewpoint(T_skewt, Td_skewt)
    #mixr_skewt = mixing_ratio_from_relative_humidity(p_skewt, T_skewt, RH_skewt)
    spec_humid = specific_humidity_from_dewpoint(p_skewt, Td_skewt)
    mixr_skewt = mixing_ratio_from_specific_humidity(spec_humid)
    #wind_dir = wind_direction(u_skewt[:].magnitude, v_skewt[:].magnitude)
    wind_dir = []
    for i in range(len(u_skewt)):
        wdir = wind_direction(u_skewt[i], v_skewt[i])
        wind_dir.append(wdir.magnitude)
    wind_dir = np.asarray(wind_dir)
    wind_spd = wind_speed(u_skewt, v_skewt)
    prof = profile.create_profile(profile='default', pres=p_skewt[::-1], hght=hgt_skewt[::-1]-hgt_skewt[-1], tmpc=T_skewt[::-1],
                                        dwpc=Td_skewt[::-1], wspd=wind_spd[::-1], wdir=wind_dir[::-1], missing=-9999, strictQC=True)

    #Storm motion with 2014 Bunkers values
    srwind = params.bunkers_storm_motion(prof)
    print("Bunker's Storm Motion (right-mover) [deg,kts]:", utils.comp2vec(srwind[0], srwind[1]))
    print("Bunker's Storm Motion (left-mover) [deg,kts]:", utils.comp2vec(srwind[2], srwind[3]))
    mumr = thermo.mixratio(prof.pres, prof.dwpc)
    sfcpcl = params.parcelx( prof, flag=1 ) # Surface Parcel
    #fcstpcl = params.parcelx( prof, flag=2 ) # Forecast Parcel
    mupcl = params.parcelx( prof, flag=3 ) # Most-Unstable Parcel
    mlpcl = params.parcelx( prof, flag=4 ) # 100 mb Mean Layer Parcel

    #This cell calculates and prints out a bunch of parameters
    #You'll need the freezing level, which is printed out below
    #as
    ###############################################################
    #Get MU and ML CAPE
    mlcape = mlpcl.bplus
    mlcape = np.round(mlcape, 1)
    mucape = mupcl.bplus
    mucape = np.round(mucape, 1)
    sbcape = sfcpcl.bplus
    sbcape = np.round(sbcape, 1)
    cape3 = mlpcl.b3km
    cape3 = np.round(cape3, 1)
    cape6 = mlpcl.b6km
    cape6 = np.round(cape6, 1)
    #cape1 = mlpcl.b1km
    #cape1 = np.round(cape1, 1)
    sfccape6 = sfcpcl.b6km
    sfccape6 = np.round(sfccape6, 1)
    sfccape3 = sfcpcl.b3km
    sfccape3 = np.round(sfccape3, 1)
    #sfccape1 = sfcpcl.b1km
    #sfccape1 = np.round(sfccape1, 1)
    mucape6 = mupcl.b6km
    mucape6 = np.round(mucape6, 1)
    mucape3 = mupcl.b3km
    mucape3 = np.round(mucape3, 1)
    #mucape1 = mupcl.b1km
    #mucape1 = np.round(mucape1, 1)

    lapse_rate = params.lapse_rate(prof, 1000., 3000., pres=False )
    lapse_rate = np.round(lapse_rate,1)
    lr_36 = params.lapse_rate(prof, 3000., 6000., pres=False)
    lr_36 = np.round(lr_36,1)
    all_lr = params.max_lapse_rate(prof, lower=2000, upper=6000, interval=250, depth=2000)
    all_lr = np.round(all_lr,1)
    print(lapse_rate)
    print(lr_36)
    #Get sfc-based lfc height
    lfc = mlpcl.lfchght
    lfc = np.round(lfc, 1)
    el = mlpcl.elhght
    el = np.round(el, 1)
    #Get cin
    sbcin = sfcpcl.bminus
    sbcin = np.round(sbcin, 1)
    mucin = mupcl.bminus
    mucin = np.round(mucin, 1)
    mlcin = mlpcl.bminus
    mlcin = np.round(mlcin, 1)
    #Get parameters to match previous data
    sfc = prof.pres[prof.sfc]
    pp5km = interp.pres(prof, interp.to_msl(prof, 500.))
    p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
    p1p5km = interp.pres(prof, interp.to_msl(prof, 1500.))
    p2km = interp.pres(prof, interp.to_msl(prof, 2000.))
    p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
    p4km = interp.pres(prof, interp.to_msl(prof, 4000.))
    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    p9km = interp.pres(prof, interp.to_msl(prof, 9000.))
    p11km = interp.pres(prof, interp.to_msl(prof, 11000.))


    srhp5km = winds.helicity(prof, 0, 500., stu = srwind[0], stv = srwind[1])[0]
    srhp5km = np.round(srhp5km, 1)
    srh1km = winds.helicity(prof, 0, 1000., stu = srwind[0], stv = srwind[1])[0]
    srh1km = np.round(srh1km, 1)
    srh2km = winds.helicity(prof, 0, 2000., stu = srwind[0], stv = srwind[1])[0]
    srh2km = np.round(srh2km, 1)
    srh3km = winds.helicity(prof, 0, 3000., stu = srwind[0], stv = srwind[1])[0]
    srh3km = np.round(srh3km, 1)
    srh13km = winds.helicity(prof, 1000., 3000., stu = srwind[0], stv = srwind[1])[0]
    srh13km = np.round(srh13km, 1)

    #Get shear magnitudes
    sfc_p5km_shear = winds.wind_shear(prof, pbot=sfc, ptop=pp5km)
    sfc_1km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p1km)
    sfc_2km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p2km)
    sfc_3km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p3km)
    sfc_6km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
    one_3km_shear = winds.wind_shear(prof, pbot=p1km, ptop=p3km)
    sfcp5shear = utils.mag( sfc_p5km_shear[0], sfc_p5km_shear[1])
    sfcp5shear = np.round(sfcp5shear,1)
    sfc1shear = utils.mag( sfc_1km_shear[0], sfc_1km_shear[1] )
    sfc1shear = np.round(sfc1shear,1)
    sfc2shear = utils.mag( sfc_2km_shear[0], sfc_2km_shear[1] )
    sfc2shear = np.round(sfc2shear,1)
    sfc3shear = utils.mag( sfc_3km_shear[0], sfc_3km_shear[1] )
    sfc3shear = np.round(sfc3shear, 1)
    sfc6shear = utils.mag( sfc_6km_shear[0], sfc_6km_shear[1] )
    sfc6shear = np.round(sfc6shear,1)
    one3shear = utils.mag( one_3km_shear[0], one_3km_shear[1] )
    one3shear = np.round(one3shear,1)

    #Calculating Storm-Relative Winds
    sr_sfc_15 = winds.sr_wind(prof, pbot=sfc, ptop=p1p5km, stu = srwind[0], stv = srwind[1])[0]
    sr_sfc_15 = np.round(sr_sfc_15, 1)
    sr_4_6 = winds.sr_wind(prof, pbot=p4km, ptop=p6km, stu = srwind[0], stv = srwind[1])[0]
    sr_4_6 = np.round(sr_4_6, 1)
    sr_9_11 = winds.sr_wind(prof, pbot=p9km, ptop=p11km, stu = srwind[0], stv = srwind[1])[0]
    sr_9_11 = np.round(sr_9_11, 1)
    sr_2_4 = winds.sr_wind(prof, pbot=p2km, ptop=p4km, stu = srwind[0], stv = srwind[1])[0]
    sr_2_4 = np.round(sr_2_4, 1)

    #Maybe Calculate some mean wind
    mean1_3 = winds.mean_wind_npw(prof, pbot=1000., ptop=3000., stu = srwind[0], stv = srwind[1])[0]
    print(mean1_3)
    #Use MetPy to get the RH interpolations and pressure-weighted layers
    mpwrh_13 = mean_pressure_weighted(p_skewt[::-1], RH_skewt[::-1], height=hgt_skewt[::-1]-hgt_skewt[-1], bottom=1000*units('meter'), depth=2000 * units.meter)[0].magnitude
    mpwrh_36 = mean_pressure_weighted(p_skewt[::-1], RH_skewt[::-1], height=hgt_skewt[::-1]-hgt_skewt[-1], bottom=3000*units('meter'), depth=3000 * units.meter)[0].magnitude
    mpwrh_69 = mean_pressure_weighted(p_skewt[::-1], RH_skewt[::-1], height=hgt_skewt[::-1]-hgt_skewt[-1], bottom=6000*units('meter'), depth=3000 * units.meter)[0].magnitude
    #mpwrh_69 = np.nan
    #Use MetPy to interpolate to the 3, 6, and 9km levels

    rh3 = metinterp(3000*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], RH_skewt[::-1])[0].magnitude
    rh6 = metinterp(6000*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], RH_skewt[::-1])[0].magnitude
    rh9 = metinterp(9000*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], RH_skewt[::-1])[0].magnitude

    #Get 0C height from SharpPy
    hgt0c = sfcpcl.hght0c
    hgt0c = np.round(hgt0c,1)
    #Get convective temp
    ctemp = params.convective_temp(prof)
    #Add to a new profile and get the lcl from that (which becomes the ccl)
    T_skewt_c = np.copy(T_skewt.magnitude)
    T_skewt_c[-1] = ctemp
    profc = profile.create_profile(profile='default', pres=p_skewt[::-1], hght=hgt_skewt[::-1]-hgt_skewt[-1], tmpc=T_skewt_c[::-1],
                                        dwpc=Td_skewt[::-1], wspd=wind_spd[::-1], wdir=wind_dir[::-1], missing=-9999, strictQC=True)
    sfcpclc = params.parcelx( profc, flag=1 ) # Surface Parcel for ccl
    ccl = mlpcl.lclhght
    ccl = np.round(ccl, 1)
    #Interpolate temperature to the ccl height with MetPy
    cclt = metinterp(ccl*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], T_skewt[::-1])[0].magnitude

    #Get normal lcl height and t
    lcl = mlpcl.lclhght
    lclt = metinterp(lcl*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], T_skewt[::-1])[0].magnitude
    lcl = np.round(lcl,1)

    mpl = mlpcl.mplhght
    mpl = np.round(mpl,1)

    #Get sigtor
    sigtor = params.stp_fixed(sfcpcl.bplus, sfcpcl.lclhght, srh1km, utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])
    sigtor = np.round(sigtor,1)
    #Get effective inflow stuff
    eff_inflow = params.effective_inflow_layer(prof)
    ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
    etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))
    effective_srh = winds.helicity(prof, ebot_hght, etop_hght, stu = srwind[0], stv = srwind[1])[0]
    ebwd = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
    ebwspd = utils.mag( ebwd[0], ebwd[1] )

    #Get supercell composite
    scp = params.scp(mupcl.bplus, effective_srh, ebwspd)
    scp = np.round(scp,1)

    #Get EHI 0-3km
    EHI = params.ehi(prof, sfcpcl, sfc, 3000, stu = srwind[0], stv = srwind[1])
    EHI = np.round(EHI)
    #Get EHI 0-1km
    EHI1 = params.ehi(prof, sfcpcl, sfc, 1000, stu = srwind[0], stv = srwind[1])
    EHI1 = np.round(EHI1)
    #Get PBL Height(depth)
    pbl = params.pbl_top(prof)
    pbl = np.round(pbl)
    def dcape(prof):
        '''
            Downdraft CAPE (DCAPE)

            Adapted from John Hart's (SPC) DCAPE code in NSHARP donated by Rich Thompson (SPC)

            Calculates the downdraft CAPE value using the downdraft parcel source found in the lowest
            400 mb of the sounding.  This downdraft parcel is found by identifying the minimum 100 mb layer
            averaged Theta-E.

            Afterwards, this parcel is lowered to the surface moist adiabatically (w/o virtual temperature
            correction) and the energy accumulated is called the DCAPE.

            Future adaptations of this function may utilize the Parcel/DefineParcel object.

            Parameters
            ----------
            prof : profile object
                Profile object

            Returns
            -------
            dcape : number
                downdraft CAPE (J/kg)
            ttrace : array
                downdraft parcel trace temperature (C)
            ptrace : array
                downdraft parcel trace pressure (mb)
            '''

        sfc_pres = prof.pres[prof.sfc]
        prof_thetae = prof.thetae
        prof_wetbulb = prof.wetbulb
        mask1 = prof_thetae.mask
        mask2 = prof.pres.mask
        mask = np.maximum( mask1, mask2 )
        prof_thetae = prof_thetae[~mask]
        prof_wetbulb = prof_wetbulb[~mask]
        pres = prof.pres[~mask]
        hght = prof.hght[~mask]
        dwpc = prof.dwpc[~mask]
        tmpc = prof.tmpc[~mask]
        idx = np.where(pres >= sfc_pres - 400.)[0]

        # Find the minimum average theta-e in a 100 mb layer
        mine = 1000.0
        minp = -999.0
        for i in idx:
            thta_e_mean = mean_thetae(prof, pbot=pres[i], ptop=pres[i]-100.)
            if utils.QC(thta_e_mean) and thta_e_mean < mine:
                minp = pres[i] - 50.
                mine = thta_e_mean

        upper = minp
        uptr = np.where(pres >= upper)[0]
        uptr = uptr[-1]

        # Define parcel starting point
        tp1 = thermo.wetbulb(upper, interp.temp(prof, upper), interp.dwpt(prof, upper))
        pe1 = upper
        te1 = interp.temp(prof, pe1)
        h1 = interp.hght(prof, pe1)
        tote = 0
        lyre = 0

        # To keep track of the parcel trace from the downdraft
        ttrace = [tp1]
        ptrace = [upper]

        # Lower the parcel to the surface moist adiabatically and compute
        # total energy (DCAPE)
        iter_ranges = range(uptr, -1, -1)
        ttraces = ma.zeros(len(iter_ranges))
        ptraces = ma.zeros(len(iter_ranges))
        ttraces[:] = ptraces[:] = ma.masked
        for i in iter_ranges:
            pe2 = pres[i]
            te2 = tmpc[i]
            h2 = hght[i]
            tp2 = thermo.wetlift(pe1, tp1, pe2)

            if utils.QC(te1) and utils.QC(te2):
                tdef1 = (tp1 - te1) / (thermo.ctok(te1))
                tdef2 = (tp2 - te2) / (thermo.ctok(te2))
                lyrlast = lyre
                lyre = 9.8 * (tdef1 + tdef2) / 2.0 * (h2 - h1)
                tote += lyre

            ttraces[i] = tp2
            ptraces[i] = pe2

            pe1 = pe2
            te1 = te2
            h1 = h2
            tp1 = tp2
        drtemp = tp2 # Downrush temp in Celsius

        return tote, ma.concatenate((ttrace, ttraces[::-1])), ma.concatenate((ptrace, ptraces[::-1]))
    #dcape = params.dcape(prof)
    #dcape = np.round(dcape,1)

    print("MLCAPE", mlcape)
    print('lcl z', lcl)
    print('lfc z', lfc)
    print('MUCAPE', mucape)
    #print('cin', cin)
    print("ESRH", effective_srh)
    print('ESHR', ebwspd)
    print("scp", scp)
    print('sigtor',sigtor)
    print('rh3',rh3)
    print('rh6',rh6)
    print('rh9',rh9)
    print("p3",p3km)
    print("p6",p6km)
    print("p9",p9km)
    print('RH 1-3', mpwrh_13)
    print('RH 3-6', mpwrh_36)
    print('RH 6-9', mpwrh_69)
    print('0C Z',hgt0c)
    print('CCL T',cclt)
    print('LCL T',lclt)
    print('DD Depth')
    print('DD Center')
    print('Cap Top')
    print('Cap Max T')
    print('Cap Max T Height')
    print('0-1 SRH', srh1km)
    print('0-3 SRH', srh3km)
    print('0-1 Shear', sfc1shear)
    print('0-3 Shear', sfc3shear)
    print('0-6 Shear', sfc6shear)
    print("EHI", EHI)

    #Print out just the freezing level for easy access
    print('0C Z',hgt0c)

    #Extract a parcel profile from SharpPy
    profile1 = mlpcl.ttrace*units('degC')
    pprofile1 = mlpcl.ptrace*units('hPa')
    vtemp_pr = prof.vtmp

    zH5 = data['Geopotential_height_isobaric'].metpy.sel(vertical=500 * units.hPa)
    zH5_crs = zH5.metpy.cartopy_crs

    #plt.style.use('bmh')

    #plt.style.use('bmh')
    # Change default to be better for skew-T
    fig = plt.figure(figsize=(22, 12))
    #fig2, ax4 = plt.figure()
    #skew = SkewT(fig)
    gs = gridspec.GridSpec(8, 6)
    gspec = fig.add_gridspec(nrows=8, ncols=6)
    skew = SkewT(fig, rotation=45, subplot=gs[:7, :3])
    ax1 = fig.add_subplot(gs[:7, 3:]) #Hodograph
    #ax2 = fig.add_subplot(gspec[7:,0]) # wetbulb and theta-e
    ax2 = fig.add_axes([0.141,0.7,0.1,0.18])
    ax4 = fig.add_axes([0.141,0.01,0.17,0.17], projection=ccrs.PlateCarree())
    #ax5 = fig.add_axes([0.312,0.01,0.219,0.17], projection=ccrs.PlateCarree())
    #ax6 = fig.add_axes([0.450,0.01,0.09,0.17], projection=ccrs.PlateCarree())
    #fig.title('GFS Sounding for St. Louis, MO')
    skew.plot(p_skewt, T_skewt, 'r', linewidth=2)
    skew.plot(p_skewt, vtemp_pr[::-1], 'r', linestyle='--',)
    skew.plot(p_skewt, Td_skewt, 'g', linewidth=2)
    skew.plot_barbs(p_skewt, u_skewt, v_skewt)

    ax2.plot(thermo.ktoc(prof.thetae), prof.hght, 'r-', label='Theta-E')
    ax2.plot(prof.wetbulb, prof.hght, 'c-', label='Wetbulb')
    ax2.set_xlabel("Temperature [C]")
    ax2.set_ylabel("Height [m above MSL]")
    ax2.legend()
    #ax2.grid()

    ax3 = inset_axes(skew.ax, '30%', '30%', loc=1)
    #ax3.plot(RH_skewt, hgt_skewt, color='darkgreen', label='RH', linewidth=2)
    ax3.plot(mumr, hgt_skewt, color='darkgreen', label='MixRatio')
    ax3.set_xlabel("Moisture")
    ax3.set_ylabel("Height [km above MSL]")
    ax3.legend()
    #ax3.grid()


    ax4.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='#c7c783', zorder=0)
    height_contour = data['Geopotential_height_isobaric'].metpy.sel(vertical=500*units.hPa).squeeze()
    height_contour = ndimage.gaussian_filter(height_contour, sigma=3, order=0)
    h_contour = ax4.contour(x, y, height_contour, colors='#9403fc', levels=range(5400, 6000, 60), linewidths=1.5)
    #h_contour.clabel(fontsize=12, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
    #ax4.set_extent((235, 295, 25, 50), crs = cartopy.crs.PlateCarree())    # Set a title and show the plot
    height_contour7 = data['Geopotential_height_isobaric'].metpy.sel(vertical=700*units.hPa).squeeze()
    height_contour7 = ndimage.gaussian_filter(height_contour7, sigma=3, order=0)
    h_contour7 = ax4.contour(x, y, height_contour7, colors='#e303fc', levels=range(2900, 3200, 60), linewidth=1.5)
    height_contour8 = data['Geopotential_height_isobaric'].metpy.sel(vertical=850*units.hPa).squeeze()
    height_contour8 = ndimage.gaussian_filter(height_contour8, sigma=3, order=0)
    h_contour8 = ax4.contour(x, y, height_contour8, colors='#fc037b', levels=range(1300, 1500, 60), linewidth=1.5)
    CAPE1 = data['Convective_available_potential_energy_surface'].squeeze()
    capep = ax4.contourf(x, y, CAPE1, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000], alpha = 0.7, cmap='RdPu')

    '''
    ax4 = inset_axes(skew.ax, '30%', '30%', loc=3)
    ax4.plot(Omega, lev, color='yellow', label='Omega', linewidth=2)
    #ax4.plot(mixr_skewt, hgt_skewt, color='lightgreen', label='Omega')
    ax4.set_xlabel("Omega")
    ax4.set_ylabel("Height [km above MSL]")
    ax4.legend()
    ax4.grid()
    '''
    #profile1 = metcalc.parcel_profile(p_skewt[::-1], T_skewt[::-1], Td_skewt[::-1]).to('degC')
    #Plot the most unstable parcel path
    #skew.plot(p_skewt[0], mupcl, 'orange', linewidth=4, linestyle = '--')
    #Plot the parcel path
    #skew.plot(pprofile1, profile1, 'gold', linewidth=2)


    #Let's try to fill between the profile and parcel path.
    #greater = T_skewt >= prof
    #skew.ax.fill_betweenx(levc, T_skewt, prof, where=greater, facecolor='blue', alpha=0.4)
    #skew.ax.fill_betweenx(levc, T_skewt, prof, where=~greater, facecolor='red', alpha=0.4)

    skew.ax.set_ylim(1020, 120)
    # Good bounds for aspect ratio
    skew.ax.set_xlim(-30, 40)
    skew.ax.set_xticklabels([-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30], size = 14)
    skew.ax.set_xlabel("Temperature [C]", alpha=0.1)
    skew.ax.set_yticklabels([100,200,300,400,500,600,700,800,900,1000],size = 14)
    skew.ax.set_ylabel("Pressure [mb]")

    def draw_heights(skew, prof):
        trans = transforms.blended_transform_factory(ax.transAxes,ax.transData)

    # Plot the height values on the skew-t, if there's an issue, inform the user.
        for hght in [1000,2000,3000,4000,5000,6000,9000,12000,15000]:
            p = tab.interp.pres(prof, tab.interp.to_msl(prof, hght))
            try:
                skew.ax.text(0.01, p_skewt, str(hght/1000) +' km -', verticalalignment='center', fontsize=9, color='r')
            except:
                print("problem plotting height label:", h)

        skew.ax.text(0.01, prof.pres[prof.sfc], 'Sfc', verticalalignment='center', fontsize=9, transform=trans, color='r')
        return skew
    def draw_effective_inflow_layer(skew, prof):
        # Plot the effective inflow layer on the Skew-T, like with the GUI (TODO: include the effective SRH on the top like in the GUI).
        trans = transforms.blended_transform_factory(ax.transAxes,ax.transData)
        skew.plot([0.2,0.3], [prof.ebottom, prof.ebottom], color='c', lw=2, transform=trans)
        skew.plot([0.25,0.25], [prof.etop, prof.ebottom], color='c', lw=2, transform=trans)
        skew.plot([0.2,0.3], [prof.etop, prof.etop], color='c', lw=2, transform=trans)
    plt.figtext( 0.15, 0.35, 'Levels (m):', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.33, 'LCL:', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.33, f'{lcl}', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.31, 'LFC:', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.31, f'{lfc}', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.29, 'EL:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.29, f'{el}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.27, 'MPL:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.27, f'{mpl}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.25, 'CCL:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.25, f'{ccl}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.23, '0C:', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.23, f'{hgt0c}', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.15, 0.21, 'PBL:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.17, 0.21, f'{pbl}', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.313, 0.17, 'CAPE:', fontsize=12,  backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.17, 'Moisture:', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.340, 0.17, 'SBP', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.358, 0.17, 'MLP', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.381, 0.17, 'MUP', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.313, 0.15, 'CAPE:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.15, 'RH<lcl:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.13, 'RH${LCL-LFC}$:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.11, 'RH${1-3km}$:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.09, 'RH${HGZ}$:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.13, 'NCAPE:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.11, '6CAPE:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.09, '3CAPE:', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.313, 0.05, '1CAPE:', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.313, 0.08, 'CINH:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.07, 'Lapse Rates:', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.383, 0.07, 'Layer CAPE:', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.452, 0.07, 'Downdraft:', fontsize=12, backgroundcolor='white', color='black', weight='bold')
    plt.figtext( 0.340, 0.15, f'{sbcape}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.358, 0.15, f'{mlcape}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.381, 0.15, f'{mucape}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.340, 0.11, f'{sfccape6}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.358, 0.11, f'{cape6}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.383, 0.11, f'{mucape6}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.340, 0.09, f'{sfccape3}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.358, 0.09, f'{cape3}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.383, 0.09, f'{mucape3}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.05, 'LR_3-6km:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.05, 'TEI:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.452, 0.03, 'DCAPE:', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.458, 0.03, f'{dcape:P}', fontsize=12, backgroundcolor='white', color='cyan')
    plt.figtext( 0.346, 0.05, f'{lapse_rate}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.03, 'LR_1-3km:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.346, 0.03, f'{lr_36}', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.313, 0.01, 'LR_max:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.381, 0.05, '$CAPE_{HGZ}$:', fontsize=12, backgroundcolor='white', color='black')
    plt.figtext( 0.381, 0.03, '$CAPE_{MPL}$:', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.340, 0.01, f'{all_lr}', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.337, 0.07, f'{sfccape1}', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.358, 0.07, f'{cape1}', fontsize=12, backgroundcolor='white', color='black')
    #plt.figtext( 0.381, 0.07, f'{mucape1}', fontsize=12, backgroundcolor='white', color='black')

    '''
    plt.figtext( 0.41, 0.10, f'{sbcape}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.37, 0.08, 'SBCIN:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.41, 0.08, f'{sbcin}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.37, 0.06, 'MLCAPE:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.41, 0.06, f'{mlcape}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.37, 0.04, 'MLCIN:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.41, 0.04, f'{mlcin}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.37, 0.02, 'MUCAPE:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.41, 0.02, f'{mucape}', fontsize=14,  backgroundcolor='white', color='black')

    plt.figtext( 0.45, 0.10, '0-1km EHI:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.50, 0.10, f'{EHI1}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.45, 0.08, '0-3km EHI:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.50, 0.08, f'{EHI}', fontsize=14,  backgroundcolor='white', color='black')

    plt.figtext( 0.45, 0.06, 'SCP:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.49, 0.06, f'{scp}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.45, 0.04, 'STP:', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.49, 0.04, f'{sigtor}', fontsize=14,  backgroundcolor='white', color='black')
    plt.figtext( 0.45, 0.02, 'LHP:', fontsize=14, backgroundcolor='white', color='black')
    #plt.figtext( 0.49, 0.02, f'{lhp}', fontsize=14, backgroundcolor='white', color='black')
    '''
    plt.figtext( 0.79, 0.85, 'Bulk Shear', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.79, 0.83, '0-0.5km:' ,fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.83, 0.83, f'{sfcp5shear}', fontsize=14, backgroundcolor='white', color='purple')
    plt.figtext( 0.79, 0.81, '0-1 km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.83, 0.81, f'{sfc1shear}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.79, 0.79, '1-3 km:',fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.83, 0.79, f'{one3shear}', fontsize=14, color='purple')
    plt.figtext( 0.79, 0.77, '0-2 km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.83, 0.77, f'{sfc2shear}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.79, 0.75, '0-3 km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.83, 0.75, f'{sfc3shear}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.79, 0.73, '0-6 km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.83, 0.73, f'{sfc6shear}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.58, 0.48, 'SR Wind:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.58, 0.46, '0-1.5km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.63, 0.46, f'{sr_sfc_15}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.58, 0.44, '2-4km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.63, 0.44, f'{sr_2_4}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.58, 0.42, '4-6km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.63, 0.42, f'{sr_4_6}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.58, 0.40, '9-11km:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext( 0.63, 0.40, f'{sr_9_11}', fontsize=14, backgroundcolor='white', color='black')

    plt.figtext ( 0.58, 0.85, 'Critical SRH', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.58, 0.83, 'SRH 0-0.5:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.63, 0.83, f'{srhp5km}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.58, 0.81, 'SRH 0-1:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.63, 0.81, f'{srh1km}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.58, 0.79, 'SRH 1-3:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.63, 0.79, f'{srh13km}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.58, 0.77, 'SRH 0-2:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.63, 0.77, f'{srh2km}', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.58, 0.75, 'SRH 0-3:', fontsize=14, backgroundcolor='white', color='black')
    plt.figtext ( 0.63, 0.75, f'{srh3km}', fontsize=14, backgroundcolor='white', color='black')

    # Add the relevant special lines
    skew.plot_dry_adiabats(linewidth=1)
    skew.plot_moist_adiabats(linewidth=1)
    skew.plot_mixing_lines(linewidth=1)

    h = Hodograph(ax1,component_range=80.)
    h.add_grid(increment=20)
    c = h.plot_colormapped(u_skewt[::-1], v_skewt[::-1], hgt_skewt[::-1]-hgt_skewt[-1])
    #cl = h.plot_colormapped(u_skewt[::-1], v_skewt[::-1], hgt_skewt[::-1], bounds = [0, 1000, 3000, 5000, 10000] * units('meter'),
    #             colors = ['magenta', 'red', 'yellow', 'green'], linewidth = 4)
    #cl = h.plot_colormapped(u_skewt[::-1], v_skewt[::-1], np.asarray(hgt_skewt[::-1]), bounds = np.asarray([0, 1000, 3000, 5000, 10000]) * units('meter'), colors = ['magenta', 'red', 'yellow', 'green'], linewidth = 4)
    #Interpolate the wind levels
    hgt_points = metinterp([500, 1000, 3000, 5000, 10000]*units('meter'),hgt_skewt[::-1]-hgt_skewt[-1], u_skewt[::-1], v_skewt[::-1])
    h.ax.scatter(storm_motion[0][0].to('knots'), storm_motion[0][1].to('knots'), s = 38, color = 'k')
    h.ax.scatter(hgt_points[0][0], hgt_points[1][0], s = 48, color = 'black', zorder=10)
    h.ax.scatter(hgt_points[0][1], hgt_points[1][1], s = 48, color = 'black', zorder=10)
    h.ax.scatter(hgt_points[0][2], hgt_points[1][2], s = 48, color = 'black', zorder=10)
    h.ax.scatter(hgt_points[0][3], hgt_points[1][3], s = 48, color = 'black', zorder=10)
    h.ax.scatter(hgt_points[0][4], hgt_points[1][4], s = 48, color = 'black', zorder=10)

    #cb = plt.colorbar(cl, shrink = .8)
    #cb.ax.set_yticklabels(['sfc','1,000','3,000','5,000','10,000'])
    #cb.set_label('Height AGL (m)')
    h.ax.set_xticklabels([-80,-60,-40,-20,0,20,40,60,80],size = 14)
    h.ax.set_yticklabels([-80,-60,-40,-20,0,20,40,60,80],size = 14)
    #plt.savefig('RandomSounding.png')
    #plt.show()
    plt.savefig((output_dir+'/Soundings/GFS/STL/'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'_v2.png'),bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.clf()
