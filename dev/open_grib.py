import xarray as xr

grib_link = 'allan_grib.grb2'

ds = xr.open_dataset(grib_link,engine='cfgrib',backend_kwargs=dict(filter_by_keys={'typeOfLevel':'isobaricInhPa'}))
print(ds)
print(ds['gh'])
