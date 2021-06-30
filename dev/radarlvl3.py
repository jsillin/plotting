# Level 3 example with multiple products
import numpy as np
import matplotlib.pyplot as plt
from numpy import ma

from metpy.cbook import get_test_data
from metpy.io.nexrad import Level3File
from metpy.plots import ctables

import wget
import matplotlib.pyplot as plt
import numpy as np
from metpy.io import Level3File
import xarray as xr
import scipy.ndimage as ndimage
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from datetime import datetime
from netCDF4 import num2date
import numpy as np
import xarray as xr
from siphon.catalog import TDSCatalog
from datetime import datetime
import datetime as dt
from xarray.backends import NetCDF4DataStore

link = ('https://tgftp.nws.noaa.gov/SL.us008001/DF.of/DC.radar/DS.141md/SI.kfws/sn.last')
wget.download(link)

import requests
print('Beginning file download with requests')

url = 'https://tgftp.nws.noaa.gov/SL.us008001/DF.of/DC.radar/DS.141md/SI.kfws/sn.last'
r = requests.get(url)

with open('url', 'wb') as f:
  f.write(r.content)

print(r.status_code)
print(r.headers['content-type'])
print(r.encoding)
# Helper code for making sense of these products. This is hidden from the slideshow
# and eventually, in some form, will make its way into MetPy proper.
def print_tab_pages(prod):
    print(('\n' + '-'*80 + '\n').join(prod.tab_pages))

def print_graph_pages(prod):
    colors = {0:'white', 3:'red', 4:'cyan'}
    for page in prod.graph_pages:
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
        ax.axesPatch.set_facecolor('black')
        for line in page:
            if 'color' in line:
                c = colors[line['color']]
                if 'text' in line:
                    ax.text(line['x'], line['y'], line['text'], color=c,
                            transform=ax.transData, verticalalignment='top',
                            horizontalalignment='left', fontdict={'family':'monospace'},
                            fontsize=8)
                else:
                    vecs = np.array(line['vectors'])
                    ax.plot(vecs[:, ::2], vecs[:, 1::2], color=c)
        ax.set_xlim(0, 639)
        ax.set_ylim(511, 0)
        ax.set_aspect('equal', 'box')
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_formatter(plt.NullFormatter())
        ax.yaxis.set_major_locator(plt.NullLocator())
        for s in ax.spines: ax.spines[s].set_color('none')

def plot_prod(prod, cmap, norm, ax=None):
    if ax is None:
        ax = plt.gca()

    data_block = prod.sym_block[0][0]
    data = np.array(data_block['data'])
    data = prod.map_data(data)
    data = np.ma.array(data, mask=np.isnan(data))
    if 'start_az' in data_block:
        az = np.array(data_block['start_az'] + [data_block['end_az'][-1]])
        rng = np.linspace(0, prod.max_range, data.shape[-1] + 1)
        x = rng * np.sin(np.deg2rad(az[:, None]))
        y = rng * np.cos(np.deg2rad(az[:, None]))
    else:
        x = np.linspace(-prod.max_range, prod.max_range, data.shape[1] + 1)
        y = np.linspace(-prod.max_range, prod.max_range, data.shape[0] + 1)
        data = data[::-1]
    pc = ax.pcolormesh(x, y, data, cmap=cmap, norm=norm)
    plt.colorbar(pc, extend='both')
    ax.set_aspect('equal', 'datalim')
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    return pc, data

def plot_points(prod, ax=None):
    if ax is None:
        ax = plt.gca()

    data_block = prod.sym_block[0]
    styles = {'MDA': dict(marker='o', markerfacecolor='None', markeredgewidth=2, size='radius'),
              'MDA (Elev.)': dict(marker='s', markerfacecolor='None', markeredgewidth=2, size='radius'),
              'TVS': dict(marker='v', markerfacecolor='red', markersize=10),
              'Storm ID': dict(text='id'),
              'HDA': dict(marker='o', markersize=10, markerfacecolor='blue', alpha=0.5)}
    artists = []
    for point in data_block:
        if 'type' in point:
            info = styles.get(point['type'], {}).copy()
            x,y = point['x'], point['y']
            text_key = info.pop('text', None)
            if text_key:
                artists.append(ax.text(x, y, point[text_key], transform=ax.transData, clip_box=ax.bbox, **info))
                artists[-1].set_clip_on(True)
            else:
                size_key = info.pop('size', None)
                if size_key:
                    info['markersize'] = np.pi * point[size_key]**2
                artists.append(ax.plot(x, y, **info))

def plot_tracks(prod, ax=None):
    if ax is None:
        ax = plt.gca()
        
    data_block = prod.sym_block[0]

    for track in data_block:
        if 'marker' in track:
            pass
        if 'track' in track:
            x,y = np.array(track['track']).T
            ax.plot(x, y, color='k')
			
# Read in a bunch of NIDS products
#tvs = Level3File(get_test_data('nids/KOUN_SDUS64_NTVFWS_201305202016'))
nmd = Level3File('sn (1).last')
#nhi = Level3File(get_test_data('nids/KOUN_SDUS64_NHIFWS_201305202016'))
#n0q = Level3File('Level3_FWS_N0Q_20210327_0222.nids')
#nst = Level3File(get_test_data('nids/KOUN_SDUS34_NSTTLX_201305202016'))

# What happens when we print one out
print(nmd)
print_tab_pages(nmd)

fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(1, 1, 1)
norm, cmap = ctables.registry.get_with_boundaries('NWSReflectivity', np.arange(0, 85, 10))
pc, data = plot_prod(n0q, cmap, norm, ax)
plot_points(tvs)
plot_points(nmd)
plot_points(nhi)
plot_tracks(nst)
ax.set_ylim(-20, 40)

plt.show()



ax.set_xlim(-50, 20)