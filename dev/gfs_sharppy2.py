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
                        supercell_composite, significant_tornado, get_layer, relative_humidity_from_dewpoint)
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

cenlon = -90
cenlat = 33
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

mdate = str(year)+str(month)+str(day)
# Create new directory to store output
output_dir = str(year)+str(month)+str(day)  #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/Soundings/GFS/STL') #create subdirectory to store GFS output like this

#Get data using siphon
for i in range(0,45):
    #Get data using siphon
    best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml?dataset=grib/NCEP/GFS/Global_0p25deg/Best')
    best_ds = best_gfs.datasets[0]
    ncss = best_ds.subset()
    query = ncss.query()
    query.lonlat_box(north=55, south=20, east=-60, west=-120).time(datetime.utcnow()+dt.timedelta(hours=3*i))
    query.accept('netcdf4')
    query.variables('Geopotential_height_isobaric')

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
    dtime = data1.variables['Geopotential_height_isobaric'].dimensions[0]
    time = data['Geopotential_height_isobaric'].metpy.time
    dlev = data1.variables['Geopotential_height_isobaric'].dimensions[1]
    dlat = data1.variables['Geopotential_height_isobaric'].dimensions[2]
    dlon = data1.variables['Geopotential_height_isobaric'].dimensions[3]
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
    mucape = mupcl.bplus
    #Get sfc-based lfc height
    lfc = mlpcl.lfchght
    lfc = np.round(lfc,1)
    #Get cin
    cin = mlpcl.bminus
    #Get parameters to match previous data
    sfc = prof.pres[prof.sfc]

    p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
    p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    p9km = interp.pres(prof, interp.to_msl(prof, 9000.))

    srh1km = winds.helicity(prof, 0, 1000., stu = srwind[0], stv = srwind[1])[0]
    srh3km = winds.helicity(prof, 0, 3000., stu = srwind[0], stv = srwind[1])[0]

    #Get shear magnitudes
    sfc_1km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p1km)
    sfc_3km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p3km)
    sfc_6km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
    sfc1shear = utils.mag( sfc_1km_shear[0], sfc_1km_shear[1] )
    sfc3shear = utils.mag( sfc_3km_shear[0], sfc_3km_shear[1] )
    sfc6shear = utils.mag( sfc_6km_shear[0], sfc_6km_shear[1] )

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
    #Interpolate temperature to the ccl height with MetPy
    cclt = metinterp(ccl*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], T_skewt[::-1])[0].magnitude

    #Get normal lcl height and t
    lcl = mlpcl.lclhght
    lclt = metinterp(lcl*units('meter'), hgt_skewt[::-1]-hgt_skewt[-1], T_skewt[::-1])[0].magnitude
    lcl = np.round(lcl,1)

    #Get sigtor
    sigtor = params.stp_fixed(sfcpcl.bplus, sfcpcl.lclhght, srh1km, utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])

    #Get effective inflow stuff
    eff_inflow = params.effective_inflow_layer(prof)
    ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
    etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))
    effective_srh = winds.helicity(prof, ebot_hght, etop_hght, stu = srwind[0], stv = srwind[1])[0]
    ebwd = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
    ebwspd = utils.mag( ebwd[0], ebwd[1] )

    #Get suepercell composite
    scp = params.scp(mupcl.bplus, effective_srh, ebwspd)
    #Get EHI
    EHI = params.ehi(prof, sfcpcl, sfc, 3000, stu = srwind[0], stv = srwind[1])

    print("MLCAPE", mlcape)
    print('lcl z', lcl)
    print('lfc z', lfc)
    print('MUCAPE', mucape)
    print('cin', cin)
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

    #plt.style.use('bmh')

    # Change default to be better for skew-T
    fig = plt.figure(figsize=(24, 12))
    #skew = SkewT(fig)
    gs = gridspec.GridSpec(8, 4)
    skew = SkewT(fig, rotation=45, subplot=gs[:, :2])
    #skew.title('GFS Sounding for St. Louis, MO')
    skew.plot(p_skewt, T_skewt, 'r')
    skew.plot(p_skewt, vtemp_pr[::-1], 'r', linestyle='--')
    skew.plot(p_skewt, Td_skewt, 'g')
    skew.plot_barbs(p_skewt, u_skewt, v_skewt)

    #profile1 = metcalc.parcel_profile(p_skewt[::-1], T_skewt[::-1], Td_skewt[::-1]).to('degC')
    #Plot the most unstable parcel path
    #skew.plot(pres_mu, mu_profile, 'orange', linewidth=4, linestyle = '--')
    #Plot the parcel path
    skew.plot(pprofile1, profile1, 'gold', linewidth=2)


    #Let's try to fill between the profile and parcel path.
    #greater = T_skewt >= prof
    #skew.ax.fill_betweenx(levc, T_skewt, prof, where=greater, facecolor='blue', alpha=0.4)
    #skew.ax.fill_betweenx(levc, T_skewt, prof, where=~greater, facecolor='red', alpha=0.4)

    skew.ax.set_ylim(1020, 120)
    # Good bounds for aspect ratio
    skew.ax.set_xlim(-30, 40)
    skew.ax.set_xticklabels([-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30], size = 14)
    skew.ax.set_xlabel("Temperature [C]")
    skew.ax.set_yticklabels([100,200,300,400,500,600,700,800,900,1000],size = 14)
    skew.ax.set_ylabel("Pressure [mb]")

    plt.figtext( 0.15, 0.30, 'Levels (m):', fontsize=10)
    plt.figtext( 0.15, 0.28, 'LCL:', fontsize=10)
    plt.figtext( 0.17, 0.28, f'{lcl}', fontsize=10)
    plt.figtext( 0.15, 0.26, 'LFC:', fontsize=10)
    plt.figtext( 0.17, 0.26, f'{lfc}', fontsize=10)
    plt.figtext( 0.15, 0.24, '0C:', fontsize=10)
    plt.figtext( 0.17, 0.24, f'{hgt0c}', fontsize=10)

    # Add the relevant special lines
    skew.plot_dry_adiabats(linewidth=1)
    skew.plot_moist_adiabats(linewidth=1)
    skew.plot_mixing_lines(linewidth=1)

    ax1 = fig.add_subplot(gs[:6, 2:])
    h = Hodograph(ax1,component_range=80.)
    h.add_grid(increment=20)
    c = h.plot_colormapped(u_skewt[::-1], v_skewt[::-1], hgt_skewt[::-1]-hgt_skewt[-1])
    #cl = h.plot_colormapped(u_skewt[::-1], v_skewt[::-1], hgt_skewt[::-1], bounds = [0, 1000, 3000, 5000, 10000] * units('meter'),
    #             colors = ['magenta', 'red', 'yellow', 'green'], linewidth = 4)
    #cl = h.plot_colormapped(u_skewt[::-1], v_skewt[::-1], np.asarray(hgt_skewt[::-1]), bounds = np.asarray([0, 1000, 3000, 5000, 10000]) * units('meter'), colors = ['magenta', 'red', 'yellow', 'green'], linewidth = 4)
    #Interpolate the wind levels
    hgt_points = metinterp([500, 1000, 3000, 5000, 10000]*units('meter'),hgt_skewt[::-1]-hgt_skewt[-1], u_skewt[::-1], v_skewt[::-1])
    h.ax.scatter(storm_motion[0][0].to('knots'), storm_motion[0][1].to('knots'), s = 38, color = 'k')
    h.ax.scatter(hgt_points[0][0], hgt_points[1][0], s = 48, color = 'pink', zorder=10)
    h.ax.scatter(hgt_points[0][1], hgt_points[1][1], s = 48, color = 'r', zorder=10)
    h.ax.scatter(hgt_points[0][2], hgt_points[1][2], s = 48, color = 'gold', zorder=10)
    h.ax.scatter(hgt_points[0][3], hgt_points[1][3], s = 48, color = 'green', zorder=10)
    h.ax.scatter(hgt_points[0][4], hgt_points[1][4], s = 48, color = 'blue', zorder=10)

    #cb = plt.colorbar(cl, shrink = .8)
    #cb.ax.set_yticklabels(['sfc','1,000','3,000','5,000','10,000'])
    #cb.set_label('Height AGL (m)')
    h.ax.set_xticklabels([-80,-60,-40,-20,0,20,40,60,80],size = 14)
    h.ax.set_yticklabels([-80,-60,-40,-20,0,20,40,60,80],size = 14)
    #plt.savefig('RandomSounding.png')
    #plt.show()
    plt.savefig((output_dir+'/Soundings/GFS/STL/'+time[0].dt.strftime('%Y-%m-%d %H%M').item()+'.png'),bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.clf()
