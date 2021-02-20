import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from datetime import datetime
import datetime as dt

current = datetime.now()

dtfs = datetime.now().strftime('%Y-%m-%d_%H%MZ')

file = 'C:/Users/Jack/Documents/kestrel_2-7-21_4pm.csv'

data = pd.read_csv(file)
print(data)

data = data.set_index('FORMATTED DATE_TIME')
print(data)

temp = data['Temperature']
dpt = data['Dew Point']
rh = data['Relative Humidity']
wbt = data['Wet Bulb Temp']

print(temp)
fig = plt.figure(figsize=(15,15))
ax1 = fig.add_subplot(111)
wbt.plot(color='cornflowerblue',legend=False,ax=ax1)
temp.plot(color='r',legend=False,ax=ax1)
dpt.plot(color='g',legend=False,ax=ax1)
rh.plot(color='teal',legend=False,ax=ax1,secondary_y=True)

cfblue = mpatches.Patch(color='cornflowerblue',label='Wet Bulb')
teal = mpatches.Patch(color='teal', label='Relative Humidity')
green = mpatches.Patch(color='g', label='Dew Point')
red = mpatches.Patch(color='r',label='Temperature')
leg = ax1.legend(handles=[red,green,cfblue,teal],loc=4,framealpha=1)
leg.set_zorder(100)

#ax1.legend(loc='lower right')
ax1.set_title('Yarmouth Maine Kestrel Observations')
ax1.set_ylabel('Degrees F')
ax1.set_xlabel('Time')
plt.savefig('raw_temp_plot'+dtfs+'.png',bbox_inches='tight',pad_inches=0.2)

wbts = wbt.rolling(6).mean()
temps = temp.rolling(6).mean()
dpts = dpt.rolling(6).mean()
rhs = rh.rolling(6).mean()

fig4 = plt.figure(figsize=(15,15))
ax4 = fig4.add_subplot(111)
ax4.set_title('Yarmouth Maine Kestrel Observations')
ax4.set_ylabel('Degrees F')
ax4.set_xlabel('Time')
wbts.plot(color='cornflowerblue',legend=False,ax=ax4)
temps.plot(color='r',legend=False,ax=ax4)
dpts.plot(color='g',legend=False,ax=ax4)
rhs.plot(color='teal',legend=False,ax=ax4,secondary_y=True)
ax4.axhline(32,color='purple',linewidth=2)
leg = ax4.legend(handles=[red,green,cfblue,teal],loc=4,framealpha=1)
leg.set_zorder(100)
plt.savefig('smoothed_temp_plot'+dtfs+'.png',bbox_inches='tight',pad_inches=0.2)


fig2 = plt.figure(figsize=(15,15))
ax2 = fig2.add_subplot(111)
rh.plot(color='darkgreen',legend=True,ax=ax2)
ax2.set_title('Yarmouth Maine Kestrel Observations')
ax2.set_ylabel('%')
ax2.set_xlabel('Time')
ax2.set_ylim(bottom=0,top=100)
plt.savefig('rh_plot'+dtfs+'.png',bbox_inches='tight',pad_inches=0.2)

wind = data['Wind Speed'].rolling(20).mean()
fig3 = plt.figure(figsize=(15,15))
ax3 = fig3.add_subplot(111)
wind.plot(color='purple',legend=True,ax=ax3)
ax3.set_title('Yarmouth Maine Kestrel Observations')
ax3.set_ylabel('mph')
ax3.set_xlabel('Time')
plt.savefig('wind_plot_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.2)
