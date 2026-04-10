# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot dendogram"

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram

path = '/home/mda_silv/users/FPS_SESA'

skip_list_inmet_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
skip_list_inmet_ii = [2, 3, 4, 14, 19, 20, 21, 24, 25, 26, 27, 28, 32, 33, 34, 35, 38, 40, 41, 44, 45, 48, 52, 54, 55, 56, 59, 60, 62, 64, 68, 
70, 77, 79, 80, 82, 83, 92, 93, 96, 100, 106, 107, 111, 113, 120, 127, 130, 133, 135, 136, 140, 141, 144, 152, 154, 155, 160, 161, 163, 167, 168, 
173, 177, 180, 181, 182, 184, 186, 187, 188, 193, 197, 199, 204, 206, 207, 210, 212, 215, 216, 219, 220, 224, 225, 226, 229, 233, 237, 239, 240, 
241, 243, 248, 249, 251, 253, 254, 256, 261, 262, 264, 266, 269, 275, 276, 277, 280, 281, 282, 293, 295, 296, 298, 300, 303, 306, 308, 314, 315, 
316, 317, 319, 322, 325, 330, 331, 334, 337, 341, 344, 347, 348, 350, 353, 354, 357, 358, 360, 361, 362, 364, 370, 383, 384, 385, 389, 390, 392, 
393, 395, 396, 400, 401, 402, 404, 405, 408, 415, 416, 418, 423, 424, 427, 434, 440, 441, 443, 446, 448, 450, 451, 454, 455, 459, 465, 467, 471, 
474, 477, 481, 483, 488, 489, 492, 496, 504, 509, 513, 514, 516, 518, 519, 520, 523, 526, 528, 534, 538, 541, 544, 546, 552, 553, 557, 559]

skip_list_smn_ii = [39, 51, 55, 58, 64, 65, 66, 72, 75, 83, 86, 90, 91, 92]

def import_inmet():
		
	lat, lon, clim_i = [], [], []		  
	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue
		yy = float(inmet[i][2]) 
		xx = float(inmet[i][3])  
		if xx <= -48 and yy <= -16.5:
			lat.append(yy)
			lon.append(xx)

			print('Reading weather station:', i, inmet[i][0], inmet[i][1])
			d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/pre/'.format(path) + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			clim_i.append(values_i)
			
	return lat, lon, clim_i


def import_smn_i():
	
	lat, lon, clim_i = [], [], []		  
	for i in range(1, 73):
		lat.append(smn_i[i][1])
		lon.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
	
	return lat, lon, clim_i


def import_smn_ii():
	
	lat, lon, clim_i = [], [], []		  
	for i in range(1, 110):
		if i in skip_list_smn_ii:
			continue
		lat.append(smn_ii[i][1])
		lon.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])
		d_i = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i)
		
	return lat, lon, clim_i
	
	
# Import latitude, longitude and database
lat_x, lon_x, clim_i_x = import_inmet()			
lat_y, lon_y, clim_i_y = import_smn_i()			
lat_z, lon_z, clim_i_z = import_smn_ii()			

lon_xx = lon_x + lon_y + lon_z
lat_yy = lat_x + lat_y + lat_z
clim_tot = clim_i_x + clim_i_y + clim_i_z

# Create the instituction list
inst_all = (
    ['INMET'] * len(clim_i_x) +
    ['SMN I'] * len(clim_i_y) +
    ['SMN II'] * len(clim_i_z)
)

# Biulding df
df = pd.DataFrame(clim_tot)
df = df.apply(pd.to_numeric, errors='coerce')
df = df.replace([np.inf, -np.inf], np.nan).dropna()

valid_idx = df.index
lat_yy_clean = [lat_yy[i] for i in valid_idx]
lon_xx_clean = [lon_xx[i] for i in valid_idx]
inst_clean = [inst_all[i] for i in valid_idx]

# Linkage hierarchical clustering
Z = linkage(df, method='ward', metric='euclidean')
hc = AgglomerativeClustering(n_clusters=5, metric='euclidean', linkage='ward')
labels = hc.fit_predict(df)

print(hc)
print()
print(labels)

# Final df
df_clusters = pd.DataFrame({
    'Cluster': labels,
    'Latitude': lat_yy_clean,
    'Longitude': lon_xx_clean,
    'Institution': inst_clean
})

# INMET-only cluster list
df_inmet = df_clusters[df_clusters['Institution'] == 'INMET']

print('\nClusters containing INMET stations:')
print(sorted(df_inmet['Cluster'].unique()))

print('\nINMET stations per cluster:')
for cl in sorted(df_inmet['Cluster'].unique()):
    subset = df_inmet[df_inmet['Cluster'] == cl]
    print('\nCluster {}'.format(cl))
    print('Number of stations: {}'.format(len(subset)))
    print(subset[['Latitude', 'Longitude']])
    
# Plot figure  
plt.figure(figsize=(30, 10))

cluster_colors = ['blue', 'red', 'green', 'gray', 'orange']
n_leaves = df.shape[0]

node_members = {i: {i} for i in range(n_leaves)}

for i, row in enumerate(Z):
    left, right = int(row[0]), int(row[1])
    node_members[n_leaves + i] = node_members[left] | node_members[right]

def link_color_func(node_id):
    members = node_members.get(int(node_id), set())
    if not members:
        return 'k'
    vals, counts = np.unique(labels[list(members)], return_counts=True)
    return cluster_colors[int(vals[np.argmax(counts)])]

d_data = dendrogram(Z, leaf_rotation=90., leaf_font_size=8., link_color_func=link_color_func, color_threshold=None)
plt.xlabel('Weather stations (INMET + SMN)', fontsize=20)
plt.ylabel('Euclidean distance', fontsize=20)

# Create legend mapping
for i, color in enumerate(cluster_colors):
	plt.plot([], [], color=color, label='Cluster {}'.format(i + 1))
plt.legend(loc='upper right', fontsize=20)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_dendrogram_pr_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
