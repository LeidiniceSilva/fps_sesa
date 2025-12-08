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

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []
kx = []
ky = []

for i in range(1, 99):

	ix.append(inmet[i][3])
	iy.append(inmet[i][2])

for j in range(1, 72):

	jx.append(smn_i[j][2])
	jy.append(smn_i[j][1])
	
for k in range(1, 86):

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
