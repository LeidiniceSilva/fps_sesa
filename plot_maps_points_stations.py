# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot point for each inmet automatic station over sesa domain"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []
kx = []
ky = []

for i in range(1, 289):
	ix.append(inmet[i][3])
	iy.append(inmet[i][2])

for j in range(1, 88):
	jx.append(arg_emas[j][1])
	jy.append(arg_emas[j][2])

for k in range(1, 72):
	kx.append(urug_smn[k][2])
	ky.append(urug_smn[k][1])
	
# Plot my map 
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-90., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=20., resolution='c')
my_map.drawmeridians(np.arange(-90.,-20.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,30.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(ix, iy, 'o', color='b', label='INMET', markersize=2)
my_map.plot(jx, jy, 'o', color='r', label='EMAS', markersize=2)
my_map.plot(kx, ky, 'o', color='g', label='SMN', markersize=2)
plt.legend(loc=1, fontsize=10)

path = '/home/nice/Documentos/github_projects/shp'
my_map.readshapefile('{0}/lim_pais/lim_pais'.format(path), 'world', drawbounds=True, color='black', linewidth=0.5)

plt.title('Automatic Weather Stations')
plt.text(-68, -36, u'SESA', fontsize=10)
plt.text(-36, -57, u'\u25B2 \nN', fontsize=10, fontweight='bold')

x1,i1 = my_map(-70,-38)
x2,i2 = my_map(-70,-14)
x3,i3 = my_map(-38,-14)
x4,i4 = my_map(-38,-38)

poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_maps_points_station.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

