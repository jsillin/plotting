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

era5 = TDSCatalog('https://rda.ucar.edu/thredds/catalog/files/g/ds633.0/e5.oper.an.sfc/197901/catalog.html')
print(list(era5.datasets))
best_ds = era5.datasets[0]
ncss = best_ds.subset()
query = ncss.query()
query.lonlat_box(north=43, south=35, east=-100, west=-111).time(datetime.utcnow())
query.accept('netcdf4')
query.variables('2t')
data = ncss.get_data(query)
print(list(data.variables))

#dataset = xr.open_dataset('https://rda.ucar.edu/thredds/catalog/files/g/ds633.0/e5.oper.an.sfc/197901/catalog.html?dataset=files/g/ds633.0/e5.oper.an.sfc/197901/e5.oper.an.sfc.128_167_2t.ll025sc.1979010100_1979013123.nc')
#print(dataset)
