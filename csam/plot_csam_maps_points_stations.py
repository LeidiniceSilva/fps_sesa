# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/01/2023"
__description__ = "This script plot point for each inmet automatic station over csam domain"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_csam_inmet_stations import inmet

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []
kx = []
ky = []

for i in range(1, 289):

	if i == 4:
		continue
	if i == 12:
		continue
	if i == 45:
		continue
	if i == 55:
		continue
	if i == 77:
		continue
	if i == 85:
		continue
	if i == 90:
		continue
	if i == 98:
		continue
	if i == 99:
		continue
	if i == 100:
		continue
	if i == 118:
		continue
	if i == 122:
		continue
	if i == 130:
		continue
	if i == 135:
		continue
	if i == 151:
		continue
	if i == 155:
		continue
	if i == 159:
		continue
	if i == 160:
		continue
	if i == 163:
		continue
	if i == 164:
		continue
	if i == 181:
		continue
	if i == 183:
		continue
	if i == 186:
		continue
	if i == 187:
		continue
	if i == 188:
		continue
	if i == 209:
		continue
	if i == 216:
		continue
	if i == 228:
		continue
	if i == 236:
		continue
	if i == 245:
		continue
	if i == 246:
		continue
	if i == 262:
		continue
	if i == 268:
		continue
	if i == 273:
		continue
	if i == 287:
		continue
			
	ix.append(inmet[i][3])
	iy.append(inmet[i][2])
	
# Plot my map 
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-75., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
my_map.drawmeridians(np.arange(-75.,-48.,4.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(ix, iy, 'o', color='b', label='INMET', markersize=5)
plt.legend(loc=2, fontsize=10)

path = '/home/nice/Documentos/github_projects/shp'
my_map.readshapefile('{0}/lim_pais/lim_pais'.format(path), 'world', drawbounds=True, color='black', linewidth=0.5)

plt.title('INMET automatic weather stations over CSAM domain', fontsize=10, fontweight='bold')
plt.text(-74, -34, u'\u25B2 \nN', fontsize=10, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/csam'
name_out = 'pyplt_maps_points_station_csam.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

