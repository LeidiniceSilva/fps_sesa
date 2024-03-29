# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot point for the INMET weather station"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet
from dict_urug_smn_stations import urug_smn

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
		
# Plot my map 
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(ix, iy, 'o', color='b', label='INMET', markersize=5)
my_map.plot(jx, jy, 'o', color='r', label='SMN', markersize=5)
plt.xlabel(u'Longitude', labelpad=20, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontweight='bold')
plt.legend(loc=4, fontsize=10)

path = '/home/nice/Documentos/github_projects/shp'
my_map.readshapefile('{0}/lim_pais/lim_pais'.format(path), 'world', drawbounds=True, color='black', linewidth=0.5)

plt.title('Weather stations over SESA domain', fontsize=10, fontweight='bold')
plt.text(-59.8, -34, u'\u25B2 \nN', fontsize=10, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_maps_points_station_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

