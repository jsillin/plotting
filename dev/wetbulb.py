import numpy as np
from metpy.units import units
import metpy.calc as mpcalc

def wetbulb_with_nan(pressure,temperature,dewpoint):
    nan_mask = np.isnan(pressure) | np.isnan(temperature) | np.isnan(dewpoint)
    idx = np.arange(pressure.size)[~nan_mask]
    wetbulb_valid_only = mpcalc.wet_bulb_temperature(pressure[idx], temperature[idx], dewpoint[idx])
    wetbulb_full = np.full(pressure.size, np.nan) * wetbulb_valid_only.units
    wetbulb_full[idx] = wetbulb_valid_only

    return wetbulb_full
