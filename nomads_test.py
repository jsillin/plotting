#MAX CAPE from NAM
import metpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime
import datetime as dt
import xarray as xr

startTime=datetime.now()


#SET DATE AND TIME FOR MODEL RUN. I HAVE ANOTHER SCRIPT THAT CHOOSES THE LATEST KINDOF.


m_date='20200706'
m_hour='12'


#SET YOUR LEVELS FOR EXPECTED RANGE.

cape_levs=np.arange(1000,7000,50)

#SET URL FOR NETCDF4 FETCH FROM THE 12K NAM

url='https://nomads.ncep.noaa.gov:9090/dods/'\
     'nam/nam'+m_date+'/nam_'+m_hour+ 'z'


#SET LATITUDE AND LONGITUDE BOUNDARIES FOR DATA''
#THE SMALLER THE FASTER
#MAKE SURE TO NOTE UNDER THE NOMADS GRADS INFO HEADING THE MAX AND MIN VALUES AVAILABLE


#SET LAT/lON FOR THIS FETCH

slat=35.0
nlat=50.0
elon=-65.0
wlon=-80.0
lat_trim=-1

#SET NUMBER OF "PAGES"
# W IS THE FIRST, L IS THE LAST. NETCDF4 STARTS AT ZERO INDEX FOR THE FIRST RESULT.
#CHECK INFO FOR NUMBER OF AVAILABLE TIME STAMPS PER MODEL

w=0
l=28

#SET VARIABLES THIS WAY TO DEFINE THE GRID POINT LONGITUDE AND LATITUDE
#NOTE THESE NUMBERS ARE ON THE INFO PAGE AT NOMADS AND ARE DIFFRENT FOR EACH MODEL

iy1 = int((slat - 12.219908) / (61.2055625455-12.219908) * 442.0)
iy2 = int((nlat - 12.219908) / (61.2055625455-12.219908) * 442.0)

ix1 = int((wlon + 152.878623) / (152.878623-49.4726308108) * 912.0)
ix2 = int((elon + 152.878623) / (152.878623-49.4726308108) * 912.0)


#PRINT EVERYTHING OK IF YOU RUN WINDOWS AND JUST GET A BLANK BLACK WINDOW *OPTIONAL


print('Max Cape.py Init...')


#FETCH THE VARIABLES. lEV ONLY APPLIES TO WINDS OR OTHER VARIABLES LAYERED
#ADD LEV SLICE AFTER W:L IF YOU WANT WIND AT A CERTAIN HEIGHT.
# cape  =     file.variables[dat3][w:l,LEV,iy1:iy2,ix1:ix2]
#THIS WOULDNT WORK HERE, CAPE DOES NOT HAVE ADDITONAL LEVELS
#BUT FOR WINDS, YOU WOULD SUBSTITUE LEV FOR THE LEVEL OF WIND YOU WANTED, 1000 MB BEING ZERO
#NOW LETS GRAB SOME DATA

# 1000mb = lev 0
# 925mb  = lev 3
# 850mb  = lev 6
# 700mb  = lev 12
# 500mb  = lev 20
# 300mb  = lev 28
# 250mb  = lev 30
# 200mb  = lev 32

iso_lev = 20

dat3='capesfc'
dat4='hgtprs'
dat5='refcclm'
dat6='ugrdprs'
dat7='vgrdprs'

file = netCDF4.Dataset(url)
lev  = file.variables['lev'][:]
t    = file.variables['time'][:]
lat  = file.variables['lat'][iy1:iy2]
lon  = file.variables['lon'][ix1:ix2]
cape = file.variables[dat3][w:l,iy1:iy2,ix1:ix2]
gph = file.variables[dat4][w:l,iso_lev,iy1:iy2,ix1:ix2]
radar = file.variables[dat5][w:l,iy1:iy2,ix1:ix2]
u = file.variables[dat6][w:l,iso_lev,iy1:iy2,ix1:ix2]
v = file.variables[dat7][w:l,iso_lev,iy1:iy2,ix1:ix2]

u = u*1.94384449
v = v*1.94384449

gph = xr.DataArray(gph)
print(lev)

print(type(t[0]))

#NOW LETS MAKE A PLOTTING GRID
lon_2d, lat_2d = np.meshgrid(lon, lat)
#IF YOU WANT THE MAX VALUE FROM EACH TIME STAMPE DO THIS
#AND THEN JUST PLOT CAPE[M]
#m=cape[0]*0
#for xx in range((l-w)):
#    m=np.where(cape[xx]>m,cape[xx],m)

wind_slice = slice(10, -10, 10)


#ADD STATE BOUNDARIES

states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                name='admin_1_states_provinces_lakes',
                                                scale='50m',
                                                facecolor='none')
#I DONT DO SUB AXIS PLOTS. i PLOT THIS WAY. tO EACH THEIR OWN I GUESS
#BTW 12.8/6,4 RENDERS A GIF THAT TWITTER WILL EXECT IF YOU ANIMATE IT
#I USE MAGICK TO ANIMATE THESE ON A COMMANDE LINE. lOVE IT
#YOU CAN GET THE DATA FROM THE THE T VARIABLE TO PLOT THE VALD TIME OF THE PLOT
print('Plotting....')
begin = w
end = l
fhr = begin
b=6
while fhr <= end :
    fig = plt.figure(figsize=(12.8*2,6.4*2)) #for gif size 1024X512 at 80 DPI
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([wlon, elon, slat, nlat+lat_trim])
    ax.add_feature(states_provinces, edgecolor='black', linewidth=2)
    capef = ax.contourf(lon_2d, lat_2d,cape[fhr],levels=cape_levs,cmap='RdPu', transform=ccrs.PlateCarree(),alpha=.7)
    gph_c = ax.contour(lon_2d, lat_2d,gph[fhr],colors='k', levels=range(5400, 6000, 60))
    gph_c.clabel(fontsize=8, colors='k', inline=1, inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)
    radarf = ax.contourf(lon_2d, lat_2d,radar[fhr], levels=[20, 25, 30, 35, 40, 45, 50, 55, 60, 65], alpha = 0.7, cmap = 'Blues')
    ax.barbs(lon_2d[::b,::b],lat_2d[::b,::b],u[fhr][::b,::b], v[fhr][::b,::b], transform=ccrs.PlateCarree(), length=6, color = '#ff0000')
    #ax.barbs(x,y, u, v)
    #ax.contour(lon_2d,lat_2d,u)
    cb = plt.colorbar(capef, extend='neither', shrink=0.6)
    cb.set_label('CAPE')

    plt.title('NAM 12km CAPE (J/kg), Composite Reflectivity (dBZ) Forecast Hour: '+str(fhr)+' Valid '+str(t[fhr]), loc='left',size=16)
    plt.savefig('NOMADS_test_out/nam_cape_test6_f'+str(fhr)+'.png')
    plt.close()

    fhr += 1
    timeelapsed = datetime.now()-startTime
    print(timeelapsed)

print ("Done")

exit
