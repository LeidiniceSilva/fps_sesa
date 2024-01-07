# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jan 02, 2024"
__description__ = "This script plot altimetry of weather stations "

import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_smn_i_stations_tot import smn_i
from dict_inmet_stations_tot import inmet

path = '/home/nice/Documentos'

ws_list = ['A899', 'A836', 'A802', 'A811', 'A827', 'A881', 'A804', 'A838', 'A812', 'A831', 'A832', 'A801', 'A886', 'A813', 'A809', 'A826', 'A803', 'A889', 'A884', 'A879', 'A808', 'A833', 'A840', 'A897', 'A867', 'A837', 'A829', 'A894', 'A883', 'A830', 'A866', 'A853', 'A814', 'A880', 'A852', 'A815', 'A839', 'A844', 'A845', 'A856', 'A810', 'A805', 'A865', 'A870', 'A828', 'A806', 'A863', 'A854', 'A860', 'A841', 'A868', 'A858', 'A817', 'A859', 'A876', 'A875', 'A848', 'A862', 'A874', 'A823', 'A873', 'A807', 'B804', 'B806', 'A818', 'A820', 'A714', 'A715', 'A871', 'A821', 'A752', 'A835', 'A869', 'A850', 'A718', 'A705', 'A763', 'A707', 'A768', 'A743', 'A727', 'A762', 'A747', 'A734', 'A748', 'A702', 'A756', 'A729', 'A520', 'A519', 'A724', 'A507', 'A035', 'A016', 'A003', 'A026', 'A909', 'A029', 'A033']

# Select lat and lon 
ix = []		  
iy = []
iz = []
jx = []
jy = []
jz = []

for i in range(1, 567):
	if inmet[i][0] in ws_list:
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		iz.append(inmet[i][4])

for j in range(1, 73):
	jy.append(smn_i[j][1])
	jx.append(smn_i[j][2])
	jz.append(smn_i[j][3])

lat_tot = iy + jy
lon_tot = ix + jx
alt_tot = iz + jz

# Plot figure
fig = plt.figure()

ax = fig.add_subplot(1, 2, 1)
my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
my_map.drawmeridians(np.arange(-70.,-45.,5.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,-15.,5.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=2)
my_map.plot(jx, jy, 'o', color='gray', label='SMN', markersize=2)
plt.title('(a)', loc='left', fontsize=10, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=10, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=30, fontsize=10, fontweight='bold')
plt.text(-48, -39, u'\u25B2 \nN', fontsize=8, fontweight='bold')
plt.text(-64, -19, u'SESA', color='red', fontsize=8, fontweight='bold')
plt.legend(loc=3, ncol=2, fontsize=8)

# SESA
a1,b1 = (-65,-35)
a2,b2 = (-65,-17)
a3,b3 = (-48,-17)
a4,b4 = (-48,-35)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly1)

ax = fig.add_subplot(1, 2, 2)
my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
my_map.drawmeridians(np.arange(-70.,-45.,5.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,-15.,5.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

sc=my_map.scatter(lon_tot, lat_tot, 4, alt_tot, cmap='jet', marker='o', vmin=0, vmax=1000)
plt.title('(b)', loc='left', fontsize=10, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=10, fontweight='bold')
cbar=plt.colorbar(sc, cax=fig.add_axes([0.91, 0.27, 0.018, 0.45]), extend='max')
cbar.set_label('Altimetry (meters)', fontsize=10, fontweight='bold')
cbar.ax.tick_params(labelsize=10)

# SESA
a1,b1 = (-65,-35)
a2,b2 = (-65,-17)
a3,b3 = (-48,-17)
a4,b4 = (-48,-35)
poly2 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='red', linewidth=1.)
plt.gca().add_patch(poly2)

# Path out to save figure
path_out = '{0}/FPS_SESA/figs/sesa_v2'.format(path)
name_out = 'pyplt_maps_stations_sesa_alt.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

