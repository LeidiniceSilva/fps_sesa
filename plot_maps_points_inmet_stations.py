# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot point for each inmet automatic station over sesa domain"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_inmet_stations_latlon import coord

# Select lat and lon 
xx = []		  
yy = []
for i in range(1, 289):
	xx.append(coord[i][1])
	yy.append(coord[i][0])

# Plot my map 
fig = plt.figure()

plt.title('INMET Automatic Stations')
my_map = Basemap(projection='cyl', llcrnrlon=-90., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=20., resolution='c')
my_map.drawmeridians(np.arange(-90.,-20.,10.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-60.,30.,10.), labels=[1,0,0,0], linewidth=0.5, color='black') 

my_map.plot(xx, yy, 'o', color='blue', label='Stations', markersize=2)
plt.legend()

path = '/home/nice/Documentos/github_projects/shp'
my_map.readshapefile('{0}/lim_pais/lim_pais'.format(path), 'world', drawbounds=True, color='black', linewidth=0.5)

plt.text(-68, -36, u'CSAM', fontsize=10)
plt.text(-36, -57, u'\u25B2 \nN', fontsize=10, fontweight='bold')

x1,i1 = my_map(-70,-38)
x2,i2 = my_map(-70,-14)
x3,i3 = my_map(-38,-14)
x4,i4 = my_map(-38,-38)

poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

# Path out to save figure
path_out = '/home/nice/Downloads/FPS_SESA/figs'
name_out = 'pyplt_maps_points_station.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

