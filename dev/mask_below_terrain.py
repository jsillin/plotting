import xarray as xr
from metpy.units import units
import metpy.calc as mpcalc
import numpy as np

def mask_below_terrain(spres,data,levs):
    '''
    Given a surface pressure, return data only above that pressure.

    Needs spres, a surface pressure (float)
    Needs data, a pint quantity array of temps/dew point/rh/whatever
    Needs levs, a pint quantity array of pressures
    '''
    above_ground = []
    for i in range(len(levs)):
        diff = levs[i]-spres
        if diff <0:
            above_ground.append(levs[i])
    pres_abv_ground = above_ground*units.hPa
    num_points_abv_ground = len(above_ground)
    data_abv_ground = data[-num_points_abv_ground:]
    return [data_abv_ground,pres_abv_ground]
