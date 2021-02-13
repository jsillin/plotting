import data_download
import xarray as xr
import numpy as np
import datetime
import os
import wget, cfgrib

ds = data_download.get_analysis('GFS','2019-11-26T12:00', '2019-11-28T12:00', ['Z', 'T', 'U', 'V'], time_step = 3)

ds2 = data_download.get_analysis('GFS', '2015-01-24T12:00', '2015-01-29T12:00', ['Z', 'T', 'U', 'V'], time_step = 3)
times = ds2['time']
print(times)
