# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot study area"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.path import Path
from netCDF4 import Dataset as nc
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from dict_inmet_stations import inmet
from dict_smn_stations import urug_smn
from mpl_toolkits.basemap import Basemap, cm

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []

for i in range(1, 101):
			
	ix.append(inmet[i][3])
	iy.append(inmet[i][2])

for j in range(1, 72):
	jx.append(urug_smn[j][2])
	jy.append(urug_smn[j][1])
	
# Specify directories 
dirnc = '/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp'
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

lat_start  = -35
lat_end    = -17
lon_start  = -75
lon_end    = -48

# Plot study area
fig = plt.figure(figsize=(7, 7))

ax = plt.subplot(1,1,1)
my_map = Basemap(ax=ax, llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end, resolution='c', area_thresh=10000., projection='cyl', lon_0=lonc, lat_0=latc, lat_ts=0)	
my_map.drawparallels(np.arange(lat_start, lat_end,  5.), labels=[1,0,0,0], fontsize=10, linewidth=1., color='black')
my_map.drawmeridians(np.arange(lon_start, lon_end, 5.), labels=[0,0,0,1], fontsize=10, linewidth=1., color='black')                  
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=1.)
x, y = my_map(lon,lat)

llevels = (1, 25, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000)
im = my_map.contourf(x, y, topo, llevels, cmap=plt.cm.terrain, extend='max')
plt.xlabel(u'Longitude', labelpad=20, fontsize=10, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=10, fontweight='bold')
plt.text(-74, -34, u'\u25B2 \nN', fontsize=10, fontweight='bold')
cbar = fig.colorbar(im, drawedges=True, fraction=0.030, pad=0.04, aspect=20)
cbar.set_label('Topography (meters)', fontsize=10, fontweight='bold')

my_map.plot(ix, iy, 'o', color='black', label='INMET', markersize=3)
my_map.plot(jx, jy, 'o', color='red', label='SMN', markersize=3)
plt.legend(loc=4, fontsize=9, frameon=False)
	
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper'
name_out = 'pyplt_maps_study_area_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

