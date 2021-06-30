"""
Skew-T Analysis
===============

Classic skew-T/log-p plot using data from University of Wyoming.

This example uses example data from the University of Wyoming sounding
archive for 12 UTC 31 October 2016 for Minneapolis, MN (MPX) and uses
MetPy to plot the classic skew-T with Temperature, Dewpoint, and wind
barbs.

"""

from datetime import datetime

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import pandas_dataframe_to_unit_arrays, units
import numpy as np
from siphon.simplewebservice.wyoming import WyomingUpperAir
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


######################################################################
# Set time using a datetime object and station as variables
#

#year = int(input("Enter Year: "))
#month = int(input("Enter Month: "))
#day = int(input("Enter Day: "))
#hour = int(input("Enter Hour: "))
#station = input("Enter Station Code: ")
#dt = datetime(year, month, day, hour)

dt = datetime(2013,5,31,12)
station = 'OUN'

######################################################################
# Grab Remote Data
# ----------------
#
# This requires an internet connection to access the sounding data from a
# remote server at the University of Wyoming.
#

# Read remote sounding data based on time (dt) and station
df = WyomingUpperAir.request_data(dt, station)
#df.to_csv('eas2900_final_sounding2.csv')
# Create dictionary of united arrays
data = pandas_dataframe_to_unit_arrays(df)

print(df)
######################################################################
# Isolate variables and attach units
#
top = 100
# Isolate united arrays from dictionary to individual variables
p = data['pressure'][:top]
T = data['temperature'][:top]
Td = data['dewpoint'][:top]
u = data['u_wind'][:top]
v = data['v_wind'][:top]
h = data['height'][:top]

wind_speed = df['speed'][:top].values * units.knots
wind_dir = df['direction'][:top].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir)
######################################################################
# Make Skew-T Plot
# ----------------
#
# The code below makes a basic skew-T plot using the MetPy plot module
# that contains a SkewT class.
#
print(p.shape)
print(u)
print(v)

# Calculate the LCL
lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])

#Calculate PWAT
pwat = mpcalc.precipitable_water(Td,p)
pwatout = round(pwat.to(units.inch).m,2)

# Calculate the parcel profile.
parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
#print(parcel_prof)
# Change default to be better for skew-T
fig = plt.figure(figsize=(9, 11))
# Initiate the skew-T plot type from MetPy class loaded earlier
skew = SkewT(fig, rotation=45)
#Calculate CAPE
cape = mpcalc.cape_cin(p,T,Td,parcel_prof)
capeout = int(cape[0].m)
cinout = int(cape[1].m)
if capeout !=0:
    # Shade areas of CAPE and CIN
    skew.shade_cin(p, T, parcel_prof)
    skew.shade_cape(p, T, parcel_prof)





# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r')
skew.plot(p, Td, 'g')
skew.plot_barbs(p[::3], u[::3], v[::3], y_clip_radius=0.03)

# Set some appropriate axes limits for x and y
skew.ax.set_xlim(-40, 40)
skew.ax.set_ylim(1020, 100)

# Add the relevant special lines to plot throughout the figure
skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
                       alpha=0.25, color='orangered')
skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
                         alpha=0.25, color='tab:green')
skew.plot_mixing_lines(p=np.arange(1000, 99, -20) * units.hPa,
                       linestyle='dotted', color='tab:blue')

# Plot LCL temperature as black dot
skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

# Plot the parcel profile as a black line
skew.plot(p, parcel_prof, 'k', linewidth=2)

#Find 0C layer
#env_zero = np.interp(0.0*units.degC,T,h)
#print(env_zero.to(units.feet))
#env_zero_heights = np.where(T<0)
#print(env_zero_heights)
#env_zero_nearest = np.min(h[env_zero_heights])
#print(env_zero_nearest)

# Plot a zero degree isotherm
skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
#skew.ax.axvline(-20, color='b', linestyle='--',linewidth=1)

# Create a hodograph
# Create an inset axes object that is 40% width and height of the
# figure and put it in the upper right hand corner.
#ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
#h = Hodograph(ax_hod, component_range=80.)
#h.add_grid(increment=20)
#h.plot_colormapped(u, v, h)  # Plot a line colored by pressure (height)

# Add some descriptive titles
plt.title('{} Sounding'.format(station), loc='left')
plt.title('Valid Time: {}'.format(dt), loc='right')

#Read out key values
#plt.text(-55,500,"CAPE: "+str(round(cape[0],0))+' J/kg', wrap=True)

#plt.text(-50,600,"PWAT: "+str(pwat)+' Inches', wrap=True)

#skew.ax.text(x = 0.02, y = 0.02, s = "CAPE: "+str(capeout)+' J/kg' + '\nCIN: '+str(cinout)+' J/kg', fontsize = 10, weight = 'bold',horizontalalignment='left',
#               verticalalignment='bottom', bbox=dict(facecolor='white', ec = 'black'),  transform=skew.ax.transAxes, zorder = 1000)

# Show the plot
#plt.show()
plt.savefig('elreno.png')
