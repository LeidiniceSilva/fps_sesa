# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot cluster analysis of the INMET weather station"

import os
import conda
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet
from dtaidistance import dtw, clustering
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering


def import_inmet(dt, type_cycle):
	
	ix = []		  
	iy = []
	clim_i = []
	clim_ii = []

	# Select lat and lon 
	for i in range(1, 155):
				
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])

		print('Reading INMET weather station:', i, inmet[i][0], inmet[i][1])
		# Reading INMET weather station	
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
	
		if type_cycle == 'diurnal_cycle':
			d_i = d_i.groupby('time.hour').mean('time')
			values_i = d_i.values*24
			clim_d_i = values_i.tolist()	
			clim_d_i = clim_d_i[3:]+clim_d_i[:3]
			clim_i.append(clim_d_i)
		else:
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			clim_i.append(values_i*24)
		
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_ii = d_ii.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_ii = d_ii.groupby('time.year').sum('time')
		values_ii = d_ii.values
		clim_ii.append(values_ii)		
			
	return iy, ix, clim_i, clim_ii


type_cycle = 'annual_cycle'	
dt = 'H_2018-01-01_2021-12-31'

print('Import latitude, longitude and database')
# Import latitude, longitude and database
iy, ix, clim_i, clim_ii = import_inmet(dt, type_cycle)			
		
lon_xx = ix
lat_yy = iy
clim_tot = clim_i
clim_tot_i = clim_ii

if type_cycle == 'diurnal_cycle':
	list_hc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 3, 1, 1, 1, 0, 3, 1, 3, 1, 1, 1, 3, 1, 3, 4, 3, 1, 3, 1, 3, 3, 3, 3, 1, 3, 0, 3, 3, 3, 3, 3, 4, 3,
	4, 0, 0, 0, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 3, 4, 3, 0, 4, 0, 4, 4, 4, 4, 4, 0, 4, 0, 0, 4, 0,
	4, 4, 4, 4, 4, 4, 0, 4, 0, 2, 4, 4, 2, 4, 4, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 0, 0, 0, 2, 0, 2, 2, 2, 0, 2,
	2, 2, 2, 2, 4, 0]
else:
	list_hc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 3, 1, 1, 1, 0, 3, 1, 3, 1, 1, 1, 3, 1, 3, 4, 3, 1, 3, 1, 3, 3, 3, 3, 1, 3, 0, 3, 3, 3, 3, 3, 4, 3,
	4, 0, 0, 0, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 3, 4, 3, 0, 4, 0, 4, 4, 4, 4, 4, 0, 4, 0, 0, 4, 0,
	4, 4, 4, 4, 4, 4, 0, 4, 0, 2, 4, 4, 2, 4, 4, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 0, 0, 0, 2, 0, 2, 2, 2, 0, 2,
	2, 2, 2, 2, 4, 0]

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

for count, idx in enumerate(list_hc):
	
	if idx == 0:
		count_i.append(count)

	if idx == 1:
		count_ii.append(count)
		
	if idx == 2:
		count_iii.append(count)
	
	if idx == 3:
		count_iv.append(count)
	
	if idx == 4:
		count_v.append(count)

pre_c_i = []
pre_c_ii = []
pre_c_iii = []
pre_c_iv = []
pre_c_v = []

pre_c_i_i = []
pre_c_ii_i = []
pre_c_iii_i = []
pre_c_iv_i = []
pre_c_v_i = []

for c_i in count_i:
	pre_c_i.append(clim_tot[c_i])
	pre_c_i_i.append(clim_tot_i[c_i])

for c_ii in count_ii:
	pre_c_ii.append(clim_tot[c_ii])
	pre_c_ii_i.append(clim_tot_i[c_ii])

for c_iii in count_iii:
	pre_c_iii.append(clim_tot[c_iii])
	pre_c_iii_i.append(clim_tot_i[c_iii])
	
for c_iv in count_iv:
	pre_c_iv.append(clim_tot[c_iv])
	pre_c_iv_i.append(clim_tot_i[c_iv])
	
for c_v in count_v:
	pre_c_v.append(clim_tot[c_v])
	pre_c_v_i.append(clim_tot_i[c_v])

cluster_i = np.nanmean(pre_c_i, axis=0)
cluster_ii = np.nanmean(pre_c_ii, axis=0)
cluster_iii = np.nanmean(pre_c_iii, axis=0)
cluster_iv = np.nanmean(pre_c_iv, axis=0)
cluster_v = np.nanmean(pre_c_v, axis=0)

cluster_i_i = np.nanmean(pre_c_i_i, axis=0)
cluster_ii_i = np.nanmean(pre_c_ii_i, axis=0)
cluster_iii_i = np.nanmean(pre_c_iii_i, axis=0)
cluster_iv_i = np.nanmean(pre_c_iv_i, axis=0)
cluster_v_i = np.nanmean(pre_c_v_i, axis=0)

print('Plot figure')
# Plot figure
fig = plt.figure()

ax = fig.add_subplot(2, 1, 1)

if type_cycle == 'diurnal_cycle':
	time = np.arange(0.5, 24 + 0.5)
	plt.xticks(time, ('00Z', '', '02Z', '', '04Z', '', '06Z', '', '08Z', '', '10Z', '', '12Z', '', '14Z', '', '16Z', '', '18Z', '', '20Z', '', '22Z', ''))
else:
	time = np.arange(0.5, 12 + 0.5)
	plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))

plt.plot(time, cluster_i, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'Cluster I')
plt.plot(time, cluster_ii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='gray', label = 'Cluster II')
plt.plot(time, cluster_iii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='green', label = 'Cluster III')
plt.plot(time, cluster_iv, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='red', label = 'Cluster IV')
plt.plot(time, cluster_v, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='yellow', label = 'Cluster V')
plt.yticks(np.arange(0, 11, 1))
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=10)
plt.legend(loc=9, ncol=5, fontsize=8, handletextpad=0.2, frameon=False)
	
ax = fig.add_subplot(2, 1, 2)
box_plot_data=[cluster_i_i, cluster_ii_i, cluster_iii_i, cluster_iv_i, cluster_v_i]
box = plt.boxplot(box_plot_data, patch_artist=True, labels=['Cluster I','Cluster II','Cluster III','Cluster IV', 'Cluster V'])
colors = ['blue', 'gray', 'green', 'red', 'yellow']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='k', linewidth=1.)
plt.yticks(np.arange(0, 2400, 200))
plt.ylabel('Annual mean precipitation (mm)', fontsize=10)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_stations_cluster_{0}_boxplot_sesa.png'.format(type_cycle)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

