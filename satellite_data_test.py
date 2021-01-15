import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import patheffects
import matplotlib.pyplot as plt
from metpy.io import GiniFile
from metpy.plots.ctables import registry
from metpy.units import units
from netCDF4 import num2date
import scipy.ndimage as ndimage
from siphon.catalog import TDSCatalog
import xarray as xr

satellite = xr.open_dataset('https://thredds.ucar.edu/thredds/dodsC/satellite/goes/east/grb/ABI/CONUS/Channel16/current/OR_ABI-L1b-RadC-M6C16_G16_s20201961736222_e20201961739007_c20201961739110.nc')
print(satellite)
