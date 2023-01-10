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

for j in range(1, 88):
	jx.append(arg_emas[j][1])
	jy.append(arg_emas[j][2])

for k in range(1, 72):

	if k == 2:
		continue
	if k == 4:
		continue
	if k == 9:
		continue
	if k == 17:
		continue
	if k == 24:
		continue
	if k == 26:
		continue
	if k == 28:
		continue
	if k == 31:
		continue
	if k == 32:
		continue
	if k == 40:
		continue
	if k == 43:
		continue
	if k == 57:
		continue
	if k == 67:
		continue
	if k == 72:
		continue
	if k == 76:
		continue
	if k == 81:
		continue
			
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

