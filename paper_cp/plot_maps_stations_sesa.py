# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot stations over sesa domain"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []
kx = []
ky = []

for i in range(1, 100):

	ix.append(inmet[i][3])
	iy.append(inmet[i][2])

for j in range(1, 73):

	jx.append(smn_i[j][2])
	jy.append(smn_i[j][1])
	
for k in range(1, 68):

	kx.append(smn_ii[k][2])
	ky.append(smn_ii[k][1])
	
# Plot my map 
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-90., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=20., resolution='c')
my_map.drawmeridians(np.arange(-90.,-20.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,30.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
my_map.plot(jx, jy, 'o', color='gray', label='SMN', markersize=2)
my_map.plot(kx, ky, 'o', color='gray', markersize=2)
plt.legend(loc=1, fontsize=10)

path = '/home/nice/Documentos/github_projects/shp'
my_map.readshapefile('{0}/lim_pais/lim_pais'.format(path), 'world', drawbounds=True, color='black', linewidth=0.5)

plt.title('Automatic Weather Stations')
plt.text(-68, -40, u'SESA', fontsize=10)
plt.text(-36, -57, u'\u25B2 \nN', fontsize=10, fontweight='bold')

x1,i1 = my_map(-65,-35)
x2,i2 = my_map(-65,-17)
x3,i3 = my_map(-48,-17)
x4,i4 = my_map(-48,-35)

poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_stations.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

