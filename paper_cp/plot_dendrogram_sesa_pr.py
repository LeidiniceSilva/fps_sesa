# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sept 22, 2025"
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


def import_inmet():
		
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for i in range(1, 99):
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet 
		d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/pre/'.format(path) + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i)
			
	return iy, ix, clim_i


def import_smn_i():
	
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for i in range(1, 72):
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])
		# Reading smn 
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
	
	return iy, ix, clim_i


def import_smn_ii():
	
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for i in range(1, 86):
		iy.append(smn_ii[i][1])
		ix.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])
		# Reading smn 
		d_i = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i)
		
	return iy, ix, clim_i
	

var = 'pr'
	
# Import latitude, longitude and database
lat_x, lon_x, clim_i_x = import_inmet()			
lat_y, lon_y, clim_i_y = import_smn_i()			
lat_z, lon_z, clim_i_z = import_smn_ii()			

lon_xx = lon_x + lon_y + lon_z
lat_yy = lat_x + lat_y + lat_z

# Combine all climatology data
clim_tot = clim_i_x + clim_i_y + clim_i_z
df = pd.DataFrame(clim_tot)

# Clean data
dp = df.dropna(axis=0)
dp = dp.apply(pd.to_numeric, errors='coerce')
dp = dp.replace([np.inf, -np.inf], np.nan).dropna()

valid_idx = dp.index
lat_yy_clean = [lat_yy[i] for i in valid_idx]
lon_xx_clean = [lon_xx[i] for i in valid_idx]

print("\nOriginal stations:", len(df))
print("Valid stations after cleaning:", len(dp))

# Safety check
assert len(dp) == len(lat_yy_clean) == len(lon_xx_clean), "Mismatch between data and coordinates!"

# Linkage hierarchical clustering
z = linkage(dp, method='ward', metric='euclidean')

# Agglomerative clustering
Agg_hc = AgglomerativeClustering(n_clusters=5, metric='euclidean', linkage='ward')
y_hc = Agg_hc.fit_predict(dp)

print("\nCluster labels:")
print(", ".join(map(str, y_hc)))
print()
print("Lat:", lat_yy_clean)
print()
print("Lon:", lon_xx_clean)
print()
print("Lat/Lon count:", len(lat_yy_clean), len(lon_xx_clean))
print("\nTotal clusters:", len(set(y_hc)))

# Plot figure  
plt.figure(figsize=(30, 10))

# Define 5 distinct colors (one per cluster)
cluster_colors = ['blue', 'gray', 'green', 'red', 'yellow']

n_leaves = dp.shape[0]
linkage_matrix = z.astype(float)
node_members = {i: {i} for i in range(n_leaves)}  # initial: each leaf maps to itself
for i, row in enumerate(linkage_matrix):
    left = int(row[0])
    right = int(row[1])
    new_node_id = n_leaves + i
    members = set()
    members.update(node_members[left])
    members.update(node_members[right])
    node_members[new_node_id] = members

def link_color_func(node_id):
    node_id = int(node_id)
    members = node_members.get(node_id, set())
    
    if not members:
        return 'k'
	
    labels = y_hc[list(members)]
    vals, counts = np.unique(labels, return_counts=True)
    majority_label = int(vals[np.argmax(counts)])
    
    return cluster_colors[majority_label]

d_data = dendrogram(z, leaf_rotation=90., leaf_font_size=8., link_color_func=link_color_func, no_labels=False, color_threshold=None)
plt.xlabel('Weather stations (INMET + SMN)', fontsize=20)
plt.ylabel('Euclidean distance', fontsize=20)

# Create legend mapping
for i, color in enumerate(cluster_colors):
    plt.plot([], [], color=color, label=f'Cluster {i+1}')
plt.legend(loc='upper right', fontsize=20)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_dendrogram_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
