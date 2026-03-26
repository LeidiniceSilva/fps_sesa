# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt

from matplotlib import gridspec
from matplotlib.path import Path
from netCDF4 import Dataset as nc
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

skip_list_inmet_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
skip_list_inmet_ii = [2, 3, 4, 14, 19, 20, 21, 24, 25, 26, 27, 28, 32, 33, 34, 35, 38, 40, 41, 44, 45, 48, 52, 54, 55, 56, 59, 60, 62, 64, 68, 
70, 77, 79, 80, 82, 83, 92, 93, 96, 100, 106, 107, 111, 113, 120, 127, 130, 133, 135, 136, 140, 141, 144, 152, 154, 155, 160, 161, 163, 167, 168, 
173, 177, 180, 181, 182, 184, 186, 187, 188, 193, 197, 199, 204, 206, 207, 210, 212, 215, 216, 219, 220, 224, 225, 226, 229, 233, 237, 239, 240, 
241, 243, 248, 249, 251, 253, 254, 256, 261, 262, 264, 266, 269, 275, 276, 277, 280, 281, 282, 293, 295, 296, 298, 300, 303, 306, 308, 314, 315, 
316, 317, 319, 322, 325, 330, 331, 334, 337, 341, 344, 347, 348, 350, 353, 354, 357, 358, 360, 361, 362, 364, 370, 383, 384, 385, 389, 390, 392, 
393, 395, 396, 400, 401, 402, 404, 405, 408, 415, 416, 418, 423, 424, 427, 434, 440, 441, 443, 446, 448, 450, 451, 454, 455, 459, 465, 467, 471, 
474, 477, 481, 483, 488, 489, 492, 496, 504, 509, 513, 514, 516, 518, 519, 520, 523, 526, 528, 534, 538, 541, 544, 546, 552, 553, 557, 559]

skip_list_smn_ii = [39, 51, 55, 58, 64, 65, 66, 72, 75, 83, 86, 90, 91, 92]

ix, iy = [], []	  
for i in range(1, 567):
	lon = float(inmet[i][3])  
	lat = float(inmet[i][2]) 
	if i in skip_list_inmet_i:
		continue
	if i in skip_list_inmet_ii:
		continue
	if lon <= -48 and lat <= -16.5:
		ix.append(lon)
		iy.append(lat)

jx, jy = [], []	  
for j in range(1, 73):
	jx.append(smn_i[j][2])
	jy.append(smn_i[j][1])

kx, ky = [], []	  
for k in range(1, 110):
	if k in skip_list_smn_ii:
		continue
	kx.append(smn_ii[k][2])
	ky.append(smn_ii[k][1])
	
# Specify directories 
dirnc = '/home/mda_silv/users/FPS_SESA/database/rcm/reg_usp'
domname = 'orog_CSAM-4i_ECMWF-ERA5'

# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join(dirnc,domname+'_evaluation_r1i1p1f1-USP-RegCM471_v0.nc'), mode='r')
    
lat  = RCMf.variables['xlat'][:,:]
lon  = RCMf.variables['xlon'][:,:]
topo = RCMf.variables['topo'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
RCMf.close()
	
# Creating mask of the border
number = 29
ny,nx = topo.shape
border_mask = np.full((ny, nx), np.nan)
border_mask[:number, :] = 1
border_mask[-number:, :] = 1
border_mask[:, :number] = 1
border_mask[:, -number:] = 1

lon_bounds = [-75, -45]
lat_bounds = [-38, -15]

# Plot study area
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
font_size = 10

ct=ax.contourf(lon, lat, topo, np.arange(0, 3030, 30), cmap='terrain', extend='max')
ax.plot(ix, iy, 'o', color='white', label='INMET', markersize=5, markeredgecolor='black', markeredgewidth=0.5)
ax.plot(jx, jy, 'o', color='gray', label='SMN', markersize=5, markeredgecolor='black', markeredgewidth=0.5)	
ax.plot(kx, ky, 'o', color='gray', markersize=5, markeredgecolor='black', markeredgewidth=0.5)	
ax.text(-61.15, -16.25, u'SESA', color='black', fontsize=font_size, fontweight='bold')
ax.text(-74, -17.75, u'\u25B2 \nN', color='black', fontsize=font_size, fontweight='bold')

ax.legend(loc=4, ncol=1, fontsize=8, frameon=True)

ax.set_xlabel(u'Longitude', fontweight='bold')
ax.set_ylabel(u'Latitude', fontweight='bold')
ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 5), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 5), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='k', ls='--', alpha=0.5)  

states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
ax.add_feature(states_provinces, edgecolor='0.1')
ax.add_feature(cfeat.BORDERS, linewidth=0.75)
ax.coastlines(linewidth=0.75)
	
cbar = plt.colorbar(ct, ax=ax, orientation='vertical', shrink=0.7, pad=0.05)
cbar.set_label('Topography (m)', fontweight='bold')

# Path out to save figure
path_out = '/home/mda_silv/users/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_study_area_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
