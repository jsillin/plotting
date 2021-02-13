import xarray as xr
import numpy as np
import datetime
import os
import wget, cfgrib

GFS_wget_varnames = {'T'      : ('isobaricInhPa', 't'),
                'Z'      : ('isobaricInhPa', 'gh'),
                'RelHum' : ('isobaricInhPa', 'r'),
                'PS'     : ('surface', 'sp'),
                'PSL'    : ('meanSea', 'prmsl'),
                'U'      : ('isobaricInhPa', 'u'),
                'V'      : ('isobaricInhPa', 'v')      }

GFS_thredds_varnames = {'T'      : 'Temperature_isobaric'       ,
                'Z'      : 'Geopotential_height_isobaric'       ,
                'RelHum' : 'Relative_humidity_isobaric'         ,
                'Omega'  : 'Vertical_velocity_pressure_isobaric',
                'U'      : 'u-component_of_wind_isobaric'       ,
                'V'      : 'v-component_of_wind_isobaric'       }


E5_varnames = ['T', 'Z', 'U', 'V', 'PV', 'Mont', 'Vort', 'Pres', 'RelHum', 'PS', 'PSL', 'T2m']

run_paths = {'GFS'    : 'GFS_analyses/',
             'E5is'   : 'ERA5_isentropic/',
             'E5pl'   : 'ERA5_isobaric/',
             'E5plhr' : 'ERA5_hr_isobaric/',
             'E5sfchr' : 'ERA5_hr_sfc/',
             'NAM'    : 'NAM_analyses/',
             'NAManl' : 'NAM_historical/'}

fname_templates = {'GFS'    : 'gfsanl_4_{y}{m:02d}{d:02d}{h:02d}_{v}.nc',
                   'NAM'    : 'nam_218_{y}{m:02d}{d:02d}{h:02d}_{v}.nc',
                   'E5is'   : 'era5_is_{y}{m:02d}{d:02d}_{v}.nc',
                   'E5pl'   : 'era5_pl_{y}{m:02d}{d:02d}_{v}.nc',
                   'E5plhr' : 'era5_pl_{y}{m:02d}{d:02d}_{v}.nc',
                   'E5sfchr' : 'era5_sfc_{y}{m:02d}{d:02d}_{v}.nc',
                   'NAManl' : 'namanl_218_{y}{m:02d}{d:02d}{h:02d}_{v}.nc'}

grb_templates = {'GFS'    : 'gfsanl_4_{y}{m:02d}{d:02d}{h:02d}.grb2',
                 'NAM'    : 'nam_218_{y}{m:02d}{d:02d}{h:02d}.grb2',
                 'NAManl' : 'namanl_218_{y}{m:02d}{d:02d}{h:02d}.grb2'}

url_templates = {'GFS'    : 'https://www.ncei.noaa.gov/data/global-forecast-system/access/historical/analysis/{y}{m:02d}/{y}{m:02d}{d:02d}/gfsanl_4_{y}{m:02d}{d:02d}_{s:02d}00_00{f:d}.grb2',
                 'NAM'    : 'https://www.ncei.noaa.gov/thredds/dodsC/nam218/{y}{m:02d}/{y}{m:02d}{d:02d}/nam_218_{y}{m:02d}{d:02d}_{s:02d}00_00{f:d}.grb2',
                 'NAManl' : 'https://www.ncei.noaa.gov/thredds/dodsC/namanl/{y}{m:02d}/{y}{m:02d}{d:02d}/namanl_218_{y}{m:02d}{d:02d}_{s:02d}00_00{f:d}.grb2'}

time_res_hours  = {'GFS'    : 3,
                   'E5is'   : 24,
                   'E5pl'   : 24,
                   'E5plhr' : 24,
                   'E5sfchr' : 24,
                   'NAM'    : 1,
                   'NAManl' : 1}
def guess_cache_locations():
# {{{
   read_cache_checks = ['/scratch/eas3421/adm_data/main/',
                        '/scratch/eas3421/2020/{user}/adm_data/',
                        '/data/{user}/adm_cache/']

   write_cache_checks = ['/scratch/eas3421/2020/{user}/adm_data/',
                         '/data/{user}/adm_cache/']

   username = os.getlogin()

   print(os.environ.keys())

   if 'ADM_CACHE' in os.environ.keys():
       path = os.environ['ADM_CACHE']
       read_cache_checks.append(path)
       write_cache_checks.append(path)

   read_cache = None
   write_cache = None
   for templ in read_cache_checks:
      path = templ.format(user = username)

      if os.path.exists(path):
         read_cache = path
         print("Reading data from %s." % path, )
         break
      else:
         print("%s not found. Trying next guess." % path)

   for templ in write_cache_checks:
      path = templ.format(user = username)
      if os.path.exists(path):
         write_cache = path
         if read_cache is None:
            read_cache = path
            print("Reading and writing data from %s." % write_cache)
         else:
            print("Writing new data to %s." % path)
         break
      else:
         print("%s not found. Trying next guess." % path)

   if write_cache is None:
      if read_cache is None:
         print("Warning: no data directories are found. Set 'read_cache' and 'write_cache' manually to access data.")
      else:
         print("Warning: no write-accessible data directories are found.")

   return read_cache, write_cache
# }}}

read_cache, write_cache = guess_cache_locations()

def validate_ds(dataset):
# {{{
   dss = list(run_paths.keys())
   if not dataset in dss:
      raise ValueError('%s is not a recognized dataset. Options are: %s' % (dataset, dss))
# }}}

def generate_fname(dataset, year, month, day, hour, var, write = False):
# {{{
   validate_ds(dataset)

   if write:
      if write_cache is None:
         raise ValueError('Write-cache directory is not set. Set write_cache.')
      rpath = write_cache + run_paths[dataset]
   else:
      if read_cache is None:
         raise ValueError('Write-cache directory is not set. Set write_cache.')
      rpath = read_cache + run_paths[dataset]

   fname = fname_templates[dataset]

   path = rpath + '{y}{m:02}{d:02}/'.format(y = year, m = month, d = day)
   fn = fname.format(y = year, m = month, d = day, h = hour, v = var)
   return path, fn
# }}}

def generate_fname_wget(dataset, year, month, day, hour, var, cache = 'read'):
# {{{
   validate_ds(dataset)

   if cache == 'write':
      if write_cache is None:
         raise ValueError('Write-cache directory is not set. Set write_cache.')
      rpath = write_cache + run_paths[dataset]
      path = rpath + '{y}{m:02}{d:02}/'.format(y = year, m = month, d = day)
      fname = fname_templates[dataset]
   elif cache == 'temp':
      if write_cache is None:
         raise ValueError('Write-cache directory is not set. Set write_cache.')
      rpath = write_cache + 'temp/'
      path = rpath
      fname = grb_templates[dataset]
   else:
      if read_cache is None:
         raise ValueError('Write-cache directory is not set. Set write_cache.')
      rpath = read_cache + run_paths[dataset]
      path = rpath + '{y}{m:02}{d:02}/'.format(y = year, m = month, d = day)
      fname = fname_templates[dataset]

   fn = fname.format(y = year, m = month, d = day, h = hour, v = var)
   return path, fn

def generate_url(dataset, year, month, day, syn, fc):
# {{{
   validate_ds(dataset)
   temp = url_templates[dataset]
   return temp.format(y = year, m = month, d = day, s = syn, f = fc)
# }}}

def find_cached(dataset, year, month, day, hour, var):
# {{{
   # Check read-cache first
   path, filename = generate_fname(dataset, year, month, day, hour, var, write = False)

   if os.path.exists(path + filename):
      return path + filename

   # If not present there, check write-cache.
   path, filename = generate_fname(dataset, year, month, day, hour, var, write = True)

   if os.path.exists(path + filename):
      return path + filename

   return None
# }}}

def open_analysis(dataset, year, month, day, hour, var):
# {{{
   filename = find_cached(dataset, year, month, day, hour, var)

   return xr.open_dataset(filename)
# }}}

def get_analysis(dataset, start_date, end_date, varlist, time_step = 3, refresh = False):
# {{{
   # Validate requested dataset
   validate_ds(dataset)

   # Validate requested variables
   for v in varlist:
      if not v in GFS_thredds_varnames.keys():
         raise ValueError('Unrecognized variable %s. Options are %s.' % (v, GFS_thredds_varnames.keys()))

   # Generate list of dates to return
   start = datetime.datetime.fromisoformat(start_date)
   if end_date is None:
      nhrs = 1
   else:
      end = datetime.datetime.fromisoformat(end_date)
      dur = end - start
      nhrs = dur.days * 24 + int(dur.seconds / 3600.)

   if nhrs > 14 * 24:
      raise ValueError('No more than 14 days can be opened at once.')

   step = max(time_res_hours[dataset], time_step)

   dates = [start + datetime.timedelta(hours = h) for h in range(0, nhrs + 1, step)]
   print('Requested %g days, %d timesteps.' % (len(dates) * step / 24., len(dates)))

   # Loop through dates, fetch any missing files
   for d in dates:
      year  = d.year
      month = d.month
      day   = d.day
      hour  = d.hour

      tofetch = []
      for v in varlist:
         if refresh or find_cached(dataset, year, month, day, hour, v) is None:
            tofetch.append(v)

      if len(tofetch) > 0:
         fetch_analysis_wget(dataset, year, month, day, hour, tofetch)

   # Open individual files and merge dataset
   Vs = []
   for v in varlist:
      Vds = []
      for d in dates:
         year  = d.year
         month = d.month
         day   = d.day
         hour  = d.hour

         Vd = open_analysis(dataset, year, month, day, hour, v)[v]
         Vds.append(Vd)

      Vs.append(xr.concat(Vds, dim = 'time'))

   print('Done.')
   return xr.merge(Vs)
# }}}

def fetch_analysis(dataset, year, month, day, hour, varlist):
# {{{
   syn = 6 * int(hour / 6)
   fc = hour - syn

   url = generate_url(dataset, year, month, day, syn, fc)

   try:
      ds = xr.open_dataset(url, decode_times = True)
   except OSError as e:
      raise OSError('{ds} Analysis for {y}-{m:02d}-{d:02d} {h:02d} UTC is unavailable ({msg}).'.format(ds = dataset, y = year, m = month, d = day, h = hour, msg = e))

   for v in varlist:
      full_name = GFS_thredds_varnames[v]
      V = ds[full_name].rename(v)

      rn = {}
      for k in V.dims:
         if 'time' in k: rn[k] = 'time'
         if 'isobaric' in k: rn[k] = 'lev'

      V = V.rename(rn)

      if 'lev' in V.dims:
         V['lev'] = V.lev / 100.
         V.lev.attrs['units'] = 'hPa'
         V.lev.attrs['long_name'] = 'Pressure'

      path, filename = generate_fname(dataset, year, month, day, hour, v, write = True)

      if not os.path.exists(path):
         os.makedirs(path)

      print('Saving %s to %s.' % (v, path + filename))
      V.to_netcdf(path + filename)
# }}}

def fetch_analysis_wget(dataset, year, month, day, hour, varlist):
# {{{
   syn = 6 * int(hour / 6)
   fc = hour - syn

   if dataset in ['E5is', 'E5pl', 'E5plhr']:
      print('Unable to automatically fetch ERA5 data at present.')
      return

   url = generate_url(dataset, year, month, day, syn, fc)

   tmp_path, tmp_filename = generate_fname_wget(dataset, year, month, day, hour, '', cache = 'temp')
   if not os.path.exists(tmp_path):
      os.makedirs(tmp_path)
   tmp_fn = tmp_path + tmp_filename

   try:
      print('Downloading GRIB file %s.' % url)
      wget.download(url, out=tmp_fn)
   except OSError as e:
      raise OSError('Downloading {ds} Analysis for {y}-{m:02d}-{d:02d} {h:02d} failed ({msg}).'.format(ds = dataset, y = year, m = month, d = day, h = hour, msg = e))

   try:
       dss = cfgrib.open_datasets(tmp_fn, decode_times = True, backend_kwargs = {'indexpath':''})
      #ds['pres'] = xr.open_dataset(tmp_fn, decode_times = True, engine='cfgrib', \
            #backend_kwargs = dict(filter_by_keys={'typeOfLevel': 'isobaricInhPa'}))
      #ds['sfc'] = xr.open_dataset(tmp_fn, decode_times = True, engine='cfgrib', \
            #backend_kwargs = dict(filter_by_keys={'typeOfLevel': 'surface'}))
      #ds['msl'] = xr.open_dataset(tmp_fn, decode_times = True, engine='cfgrib', \
            #backend_kwargs = dict(filter_by_keys={'typeOfLevel': 'meanSea'}))
   except OSError as e:
      raise OSError('Grib file for {ds} Analysis for {y}-{m:02d}-{d:02d} {h:02d} UTC is not opening correctly ({msg}).'.format(ds = dataset, y = year, m = month, d = day, h = hour, msg = e))

   tax = dict(time = np.array([dss[0].valid_time.valid_time.data]))

   for v in varlist:
      dset, full_name = GFS_wget_varnames[v]
      for d in dss:
         if dset in d.coords and full_name in d.data_vars:
            V = d[full_name].rename(v)

      V = V.expand_dims(tax)

      rn = {}
      for k in V.dims:
         if 'time' in k: rn[k] = 'time'
         if 'isobaricInhPa' in k: rn[k] = 'lev'
         if 'latitude' in k: rn[k] = 'lat'
         if 'longitude' in k: rn[k] = 'lon'

      V = V.rename(rn)

      if 'lev' in V.dims:
         #V['lev'] = V.lev / 100.
         V.lev.attrs['units'] = 'hPa'
         V.lev.attrs['long_name'] = 'Pressure'

      path, filename = generate_fname_wget(dataset, year, month, day, hour, v, cache = 'write')

      if not os.path.exists(path):
         os.makedirs(path)

      print('Saving %s to %s.' % (v, path + filename))
      V.to_netcdf(path + filename)

   #if os.path.exists(tmp_fn):
      #os.remove(tmp_fn)
# }}}
