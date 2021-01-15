import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

file = '20111418AL3120_ships.txt'
data=open(file,'r')

rows = []
for line in data:
    line=line.strip()
    rows.append(line)

print(rows[3])
title=str(rows[3])

time = rows[5].split('   ')
intensity = rows[7].split('   ')
shear = rows[11].split('   ')
ssts = rows[14].split('  ')
ml_rh = rows[20].split('   ')
lat = rows[26].split('  ')
lon = rows[27].split('  ')
ohc = rows[29].split('   ')

ri_per = rows[68][66:70].strip()
ri_shr = rows[69][66:70].strip()
ri_het = rows[70][66:70].strip()
ri_irb = rows[71][66:70].strip()
ri_mxw = rows[72][66:70].strip()
ri_ir2 = rows[73][66:70].strip()
ri_dry = rows[74][66:70].strip()
ri_pot = rows[75][66:70].strip()
ri_d20 = rows[76][66:70].strip()
ri_tpw = rows[77][66:70].strip()

ri_prob1 = rows[79][44:47].strip()
ri_prob2 = rows[80][44:47].strip()
ri_prob3 = rows[81][44:47].strip()
ri_prob4 = rows[82][44:47].strip()
ri_prob5 = rows[83][44:47].strip()
ri_prob6 = rows[84][44:47].strip()
ri_prob7 = rows[85][44:47].strip()
ri_prob8 = rows[86][44:47].strip()

ri_thres1 = rows[79][17:27].strip()
ri_thres2 = rows[80][17:27].strip()
ri_thres3 = rows[81][17:27].strip()
ri_thres4 = rows[82][17:27].strip()
ri_thres5 = rows[83][17:27].strip()
ri_thres6 = rows[84][17:27].strip()
ri_thres7 = rows[85][17:27].strip()
ri_thres8 = rows[86][17:27].strip()

ri_probs=[ri_prob1,ri_prob2,ri_prob3,ri_prob4,ri_prob5,ri_prob6,ri_prob7,ri_prob8]
ri_thres=[ri_thres1,ri_thres2,ri_thres3,ri_thres4,ri_thres5,ri_thres6,ri_thres7,ri_thres8]
print(ri_per,ri_shr,ri_het,ri_irb,ri_mxw,ri_ir2,ri_dry,ri_pot,ri_d20,ri_tpw)
ri_contributors = [ri_per,ri_shr,ri_het,ri_irb,ri_mxw,ri_ir2,ri_dry,ri_pot,ri_d20,ri_tpw]
for i in range(len(ri_contributors)):
    ri_contributors[i]=float(ri_contributors[i])
intensity_v = []
ssts_v = []
ml_rh_v = []
time_v = []
shear_v = []
lat_v = []
lon_v = []
ohc_v = []
for i in range(3,len(shear)):
    if (shear[i] != 'N/A'):
        shear_v.append(int(shear[i]))

for i in range(1,len(ml_rh)):
    if (ml_rh[i] != 'N/A'):
        ml_rh_v.append(int(ml_rh[i]))

for i in range(4, len(shear)):
    if(ssts[i] != ' N/A'):
        if ohc[i]!= 'N/A':
            ssts_v.append(float(ssts[i]))
            time_v.append(int(time[i]))
            intensity_v.append(int(intensity[i]))
            ohc_v.append(int(ohc[i]))

for i in range(4, len(lat)):
    if(lat[i] != ' N/A'):
        if(lat[i] != 'N/A'):
            lat_v.append(float(lat[i]))
            lon_v.append(float(lon[i]))

shear_va=np.array([shear_v])
ml_rh_va=np.array([ml_rh_v])
sst_va=np.array([ssts_v])
time_va=np.array([time_v])
intensity_va=np.array([intensity_v])
lat_va = np.array([lat_v])
lon_va = np.array([lon_v])
ri_va = np.array([ri_contributors])
rit_va = np.array([ri_thres])
rip_va = np.array([ri_probs])
'''
fig = plt.figure(figsize=(18,8))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)

ax1.plot(time_v,shear_v,'r',label = 'Wind Shear (kts)')
ax1.set_ylim(0,50)
ax2.plot(time_v,ml_rh_v,'g',label = 'Mid-Level Relative Humidity (%)')
ax2.set_ylim(0,100)
ax3.plot(time_v,ssts_v, 'b',label = 'Sea Surface Temperatures (C)')
ax3.set_ylim(20,35)
ax1.legend()
ax2.legend()
ax3.legend()
plt.suptitle('06z 11/1/20 SHIPS Forecast for TS Eta')
fig.text(0.1,0.4,'RI Probs',bbox=dict(boxstyle='square',ec=(1.,0.5,0.),fc=(1.,0.4,0.8)))

plt.savefig('ships_test_9_.png')
'''
fig2 = plt.figure(figsize=(18,10))
gs = fig2.add_gridspec(ncols=3,nrows=4,height_ratios=[1,1,1,3],width_ratios=[2,2,1])
gs.update(wspace=0.025,hspace=0.05)
shearax = plt.subplot(gs[0,0])
relhumax = plt.subplot(gs[1,0])
sstax = plt.subplot(gs[2,0])
riax = plt.subplot(gs[2,1])
intensitax = plt.subplot(gs[3,:])
trax = plt.subplot(gs[:2,1],projection=ccrs.PlateCarree())
ripax = plt.subplot(gs[:3,2],frameon=False)
ripax.text(.1,.85,ri_thres[0],fontsize=14)
ripax.text(.6,.85,ri_probs[0],fontsize=14)
ripax.text(.1,.75,ri_thres[1],fontsize=14)
ripax.text(.6,.75,ri_probs[1],fontsize=14)
ripax.text(.1,.65,ri_thres[2],fontsize=14)
ripax.text(.6,.65,ri_probs[2],fontsize=14)
ripax.text(.1,.55,ri_thres[3],fontsize=14)
ripax.text(.6,.55,ri_probs[3],fontsize=14)
ripax.text(.1,.45,ri_thres[4],fontsize=14)
ripax.text(.6,.45,ri_probs[4],fontsize=14)
ripax.text(.1,.35,ri_thres[5],fontsize=14)
ripax.text(.6,.35,ri_probs[5],fontsize=14)
ripax.text(.1,.25,ri_thres[6],fontsize=14)
ripax.text(.6,.25,ri_probs[6],fontsize=14)
ripax.text(.1,.15,ri_thres[7],fontsize=14)
ripax.text(.6,.15,ri_probs[7],fontsize=14)

trax.coastlines(resolution='10m')
trax.add_feature(cfeature.BORDERS.with_scale('10m'))

lomin = 360-np.min(lon_v)+3
lomax = 360-np.max(lon_v)-3
lamin = np.min(lat_v)-3
lamax = np.max(lat_v)+3
print(lomax,lomin,lamin,lamax)

for i in range(len(lat_v)):
    lo = 360.-lon_v[i]
    trax.plot([lo],[lat_v[i]],color='blue',linewidth=2,marker='o',transform=ccrs.Geodetic())

trax.set_extent((lomax,lomin,lamin,lamax))
cmap = 'RdYlGn_r'
cmap1 = 'RdYlGn'

norm_shear = mpl.colors.Normalize(vmin=0,vmax=50)
norm_rh = mpl.colors.Normalize(vmin=0,vmax=100)
norm_sst = mpl.colors.Normalize(vmin=25,vmax=30)
norm_ri = mpl.colors.Normalize(vmin=0,vmax=100)
print(ri_va)
im = shearax.imshow(shear_va, cmap=cmap,norm=norm_shear,interpolation='nearest')
im1 = relhumax.imshow(ml_rh_va, cmap=cmap1,norm=norm_rh,interpolation='nearest')
im2 = sstax.imshow(sst_va, cmap=cmap1,norm=norm_sst,interpolation='nearest')
im3 = riax.imshow(ri_va,cmap='PiYG',norm=norm_ri,interpolation='nearest')

for i in range(len(shear_v)):
    for k,j in zip(shear_v,range(len(time_v))):
        text = shearax.text(j,0,shear_va[0,j],ha='center',va='center',color='k',fontsize=11)
        text2 = relhumax.text(j,0,ml_rh_va[0,j],ha='center',va='center',color='k',fontsize=11)
        text3 = sstax.text(j,0,sst_va[0,j],ha='center',va='center',color='k',fontsize=11)

for i in range(len(ri_contributors)):
    for k,j in zip(ri_contributors,range(len(ri_contributors))):
        text4 = riax.text(j,0,ri_va[0,j],ha='center',va='center',color='k',fontsize=11)

ripax.set_title('RI Probabilities',fontsize=16)
"""for i in range(len(ri_thres)):
    for k,j in zip(ri_thres,range(len(ri_probs))):
        print(j)
        text5 = ripax.text(j,i,rit_va[0,j],ha='center',va='center',color='k',fontsize=12)
        text6 = ripax.text(j,i,rip_va[0,j],ha='center',va='center',color='k',fontsize=12)"""


shearax.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
sstax.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
relhumax.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
riax.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
ripax.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)

shearax.set_title('Environmental Indicators',fontsize=16)
riax.set_title('PER    SHR   OHC   IRB    MXW   IR2    DRY    POT    DIV    TPW')
shearax.set_ylabel('VWS')
sstax.set_ylabel('SST')
relhumax.set_ylabel('RH')
trax.set_title('Forecast Track',fontsize=16)
intensitax.set_title('Forecast Intensity',fontsize=16)
intensitax.set_ylabel('Maximum Winds (kts)')
intensitax.set_xlabel('Forecast Time')
print(time_v)
print(intensity_v)
intensitax.plot(time_va[0],intensity_va[0],linewidth=2,color='b')
intensitax.axhline(35,color='darkgray',linewidth=1,label='TS')
intensitax.axhline(64,color='darkgray',linewidth=1,label='Cat 1')
plt.suptitle('18z 11/14/20 SHIPS Forecast Guidance for TS Iota',fontsize=20)

plt.savefig('SHIPS_19.png')
