# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot altimetry of inmet stations "

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

path = '/home/mda_silv'
skip_list = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
# Select lat and lon 
ix, iy = [], []
for i in range(1, 567):
	if i in skip_list:
		continue
	iy.append(inmet[i][2])
	ix.append(inmet[i][3])

jx, jy = [], []
for j in range(1, 73):
	jy.append(smn_i[j][1])
	jx.append(smn_i[j][2])

kx, ky = [], []
for k in range(1, 110):
	ky.append(smn_ii[k][1])
	kx.append(smn_ii[k][2])
					
# Plot figure
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ax.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
ax.plot(jx, jy, 'o', color='gray', label='SMN', markersize=2)
ax.plot(kx, ky, 'o', color='gray', markersize=2)
ax.legend(loc=4, ncol=1, fontsize=8, frameon=True)

lon_bounds = [-85, -30]
lat_bounds = [-60, 15]

states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 5), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='k', ls='--', alpha=0.5)  

for label in ax.get_xticklabels() + ax.get_yticklabels():
	label.set_fontsize(font_size)
	
ax.add_feature(states_provinces, edgecolor='0.05')
ax.add_feature(cfeat.BORDERS, linewidth=0.75)
ax.coastlines(linewidth=0.75)

# CSAM
a1,b1 = (-70,-36)
a2,b2 = (-70,-16)
a3,b3 = (-45,-16)
a4,b4 = (-45,-36)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '{0}/users/FPS_SESA/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_stations_inmet_smn.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

