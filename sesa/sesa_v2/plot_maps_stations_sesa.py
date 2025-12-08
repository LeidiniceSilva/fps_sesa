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
from dict_smn_i_stations import smn_i
from dict_inmet_stations import inmet

path = '/marconi/home/userexternal/mdasilva'

# Select lat and lon 
ix = []		  
iy = []
jx = []
jy = []

for i in range(1, 100):
	if i == 13:
		continue
	if i == 18:
		continue
	if i == 19:
		continue
	if i == 24:
		continue
	if i == 28:
		continue
	if i == 29:
		continue
	if i == 43:
		continue
	if i == 44:
		continue
	if i == 63:
		continue
	if i == 64:
		continue
	if i == 77:
		continue
	if i == 79:
		continue
	if i == 82:
		continue
	if i == 98:
		continue
	iy.append(inmet[i][2])
	ix.append(inmet[i][3])

for j in range(1, 73):
	jy.append(smn_i[j][1])
	jx.append(smn_i[j][2])

lat_tot = iy + jy
lon_tot = ix + jx

# Plot figure
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
my_map.drawmeridians(np.arange(-70.,-45.,5.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,-15.,5.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('{0}//github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

my_map.plot(ix, iy, 'o', color='blue', label='INMET', markersize=4)
my_map.plot(jx, jy, 'o', color='gray', label='SMN', markersize=4)
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

# Path out to save figure
path_out = '{0}/user/mdasilva/FPS_SESA/figs/sesa_v2'.format(path)
name_out = 'pyplt_maps_stations_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

