# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/02/2023"
__description__ = "This script plot cluster analysis from each weather station"

import os
import conda
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn
from dtaidistance import dtw, clustering
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering


def import_inmet(dt):
	
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
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
		if i == 98:
			continue
		if i == 99:
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
		if i == 246:
			continue
		if i == 268:
			continue
		if i == 287:
			continue
			
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])

		print('Reading INMET weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet weather station	
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.hour').mean('time')
		values_i = d_i.values*24
		clim_d_i = values_i.tolist()	
		clim_d_ia = clim_d_i[3:]
		clim_d_ib = clim_d_i[:3]
		clim_d_iab = clim_d_ia+clim_d_ib
		clim_i.append(clim_d_iab)
	
	return iy, ix, clim_i


def import_urug_smn(dt):
	
	jy = []
	jx = []
	clim_j = []
	
	# Select lat and lon 
	for j in range(1, 72):
		
		print('Reading Uruguai weather station:', j, urug_smn[j][0])	

		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])		

		# Reading Uruguai weather stations
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_j = d_j.groupby('time.hour').mean('time')
		values_j = d_j.values*24
		clim_d_j = values_j.tolist()	
		clim_d_ja = clim_d_j[3:]
		clim_d_jb = clim_d_j[:3]
		clim_d_jab = clim_d_ja+clim_d_jb
		clim_j.append(clim_d_jab)
				
	return jy, jx, clim_j


def import_arg_emas(dt):
	
	ky = []
	kx = []
	clim_k = []
	
	# Select lat and lon 
	for k in range(1, 88):
		
		print('Reading Argentina weather station:', k, arg_emas[k][0])	

		ky.append(arg_emas[k][2])
		kx.append(arg_emas[k][1])		

		# Reading Argentina weather stations
		d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
		d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
		d_k = d_k.groupby('time.hour').mean('time')
		values_k = d_k.values*24
		clim_d_k = values_k.tolist()	
		clim_d_ka = clim_d_k[3:]
		clim_d_kb = clim_d_k[:3]
		clim_d_kab = clim_d_ka+clim_d_kb
		clim_k.append(clim_d_kab)
				
	return ky, kx, clim_k
	
	
dt = 'H_2018-01-01_2021-12-31'

print('Import latitude, longitude and database')
# Import latitude, longitude and 
iy, ix, clim_i = import_inmet(dt)			
jy, jx, clim_j = import_urug_smn(dt)
ky, kx, clim_k = import_arg_emas(dt)

lon_xx = ix+jx+kx
lat_yy = iy+jy+ky
clim_tot = clim_i+clim_j+clim_k

list_hc = [2, 2, 2, 0, 2, 2, 2, 0, 2, 0, 2, 0, 2, 2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 2, 1, 0, 1, 0, 
0, 1, 0, 0, 4, 2, 2, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 
1, 1, 1, 1, 4, 0, 2, 1, 1, 1, 1, 1, 0, 1, 1, 4, 1, 1, 1, 1, 1, 0, 1, 1, 0, 2, 1, 1, 4, 1, 1, 1, 1, 1, 3, 4, 1, 
1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 2, 4, 0, 1, 4, 1, 4, 1, 1, 1, 4, 1, 3, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 3, 1, 1, 
1, 3, 1, 1, 1, 1, 1, 2, 1, 3, 0, 1, 1, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 
1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 3, 1, 0, 
1, 1, 1, 1, 1, 1, 3, 1, 1, 2, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 1, 2, 1, 3, 1, 3, 1, 
1, 1, 1, 1, 1, 3, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 
2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 
3, 3, 3, 3, 3, 3, 3, 3, 2, 3]

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

for c_i in count_i:
	pre_c_i.append(clim_tot[c_i])

for c_ii in count_ii:
	pre_c_ii.append(clim_tot[c_ii])
	
for c_iii in count_iii:
	pre_c_iii.append(clim_tot[c_iii])
	
for c_iv in count_iv:
	pre_c_iv.append(clim_tot[c_iv])
	
for c_v in count_v:
	pre_c_v.append(clim_tot[c_v])
	
cluster_i_dc = np.nanmean(pre_c_i, axis=0)
cluster_ii_dc = np.nanmean(pre_c_ii, axis=0)
cluster_iii_dc = np.nanmean(pre_c_iii, axis=0)
cluster_iv_dc = np.nanmean(pre_c_iv, axis=0)
cluster_v_dc = np.nanmean(pre_c_v, axis=0)

cluster_i_b = [1642.996, 1549.528, 1043.268,  941.564]
cluster_ii_b = [1230.21637427, 1236.6128655,  1116.08888889,  887.4994152]
cluster_iii_b = [1179.91182796, 1490.2344086,   934.49462366,  887.99139785]
cluster_iv_b = [499.72258065, 531.37849462, 312.30752688, 334.30967742]
cluster_v_b = [2018.42, 2218.8,  1645.44, 1476.74]

print('Plot figure')
# Plot figure
fig = plt.figure()

ax = fig.add_subplot(2, 1, 1)
time = np.arange(0.5, 24 + 0.5)
plt.plot(time, cluster_i_dc, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'Cluster I')
plt.plot(time, cluster_ii_dc, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='gray', label = 'Cluster II')
plt.plot(time, cluster_iii_dc, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='green', label = 'Cluster III')
plt.plot(time, cluster_iv_dc, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='red', label = 'Cluster IV')
plt.plot(time, cluster_v_dc, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='yellow', label = 'Cluster V')
plt.yticks(np.arange(0, 11, 1))
plt.xticks(time, ('00Z', '', '02Z', '', '04Z', '', '06Z', '', '08Z', '', '10Z', '', '12Z', '', '14Z', '', '16Z', '', '18Z', '', '20Z', '', '22Z', ''))
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=10)
plt.legend(loc=1, fontsize=10)
	
ax = fig.add_subplot(2, 1, 2)
box_plot_data=[cluster_i_b, cluster_ii_b, cluster_iii_b, cluster_iv_b, cluster_v_b]
box = plt.boxplot(box_plot_data, patch_artist=True, labels=['Cluster I','Cluster II','Cluster III','Cluster IV', 'Cluster V'])
colors = ['blue', 'gray', 'green', 'red', 'yellow']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.yticks(np.arange(0, 2400, 200))
plt.ylabel('Annual mean precipitation (mm)', fontsize=10)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_stations_cluster_diurnal_cycle_boxplot.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()
