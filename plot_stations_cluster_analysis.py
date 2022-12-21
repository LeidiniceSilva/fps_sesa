# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot cluster analysis from each weather station"

import os
import conda
import cmocean
import numpy as np
import pandas as pd
import xarray as xr
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from sklearn.datasets import load_iris
from dtaidistance import dtw, clustering
from mpl_toolkits.basemap import Basemap
from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn
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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.hour').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
		
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
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_j = d_j.groupby('time.hour').mean('time')
		values_j = d_j.values
		clim_j.append(values_j*24)
				
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
		d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
		d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
		d_k = d_k.groupby('time.hour').mean('time')
		values_k = d_k.values
		clim_k.append(values_k*24)
				
	return ky, kx, clim_k
	
	
dt = 'H_2018-01-01_2021-12-31'

# Import latitude, longitude, correlation and bias
iy, ix, clim_i = import_inmet(dt)			
jy, jx, clim_j = import_urug_smn(dt)
ky, kx, clim_k = import_arg_emas(dt)
 
lon_xx = ix+jx+kx
lat_yy = iy+jy+ky
clim_tot = clim_i+clim_j+clim_k
df = pd.DataFrame(clim_tot)

# Linkage hierarchical 
z = linkage(df, method='ward', metric='euclidean')

# Agglomerative clustering
Agg_hc = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
y_hc = Agg_hc.fit_predict(df)

plt.figure(figsize=(30,10))
dendrogram(y_hc, leaf_rotation=90., leaf_font_size=6., color_threshold=3800)
plt.title('Dendrogram', fontsize=8) 
plt.xlabel('Weather automatic stations', fontsize=8) 
plt.ylabel('Euclidean distances', fontsize=8) 

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_stations_dendrogram.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')

# Hierarchical cluster        
# ~ model1 = clustering.Hierarchical(dtw.distance_matrix_fast, {})
# ~ cluster_idx = model1.fit(series)
# ~ model2 = clustering.HierarchicalTree(model1)
# ~ cluster_idx = model2.fit(series)

# ~ fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
# ~ show_ts_label = lambda idx: "Station-" + str(idx)
# ~ model2.plot(axes=ax, show_ts_label=show_ts_label, show_tr_label=True, ts_label_margin=-10, ts_left_margin=10, ts_sample_length=1)

# Linkage clustering 
# ~ model3 = clustering.LinkageTree(dtw.distance_matrix_fast, {})
# ~ cluster_idx = model3.fit(series)

# ~ fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
# ~ show_ts_label = lambda idx: "Station-" + str(idx)
# ~ model3.plot(axes=ax, show_ts_label=show_ts_label, show_tr_label=True, ts_label_margin=-10, ts_left_margin=10, ts_sample_length=1)

# K-means clustering 
# ~ from dtaidistance.clustering import kmeans
# ~ model4 = kmeans.KMeans(k=10)
# ~ cluster_idx, performed_it = model4.fit(series)

# ~ fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
# ~ show_ts_label = lambda idx: "ts-" + str(idx)
# ~ model4.plot(axes=ax, show_ts_label=show_ts_label, show_tr_label=True, ts_label_margin=-10, ts_left_margin=10, ts_sample_length=1)
           
# ~ print('Plot figure')
# ~ # Plot figure   
# ~ fig = plt.figure()

# ~ my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
# ~ my_map.drawmeridians(np.arange(-75.,-25.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
# ~ my_map.drawparallels(np.arange(-40.,5.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
# ~ my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)

# ~ pltfig = my_map.scatter(lon_xx, lat_yy, 5, clim_tot, cmap=cm.gist_rainbow, marker='o', vmin=0, vmax=12)
# ~ plt.title('(a) Precipitation climatology of diunal cycle (mm d⁻¹)', loc='left', fontsize=8)
# ~ plt.ylabel(u'Latitude', fontsize=6, labelpad=15)
# ~ plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
# ~ cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.35, 0.019, 0.28]), extend='max')
# ~ cbar.ax.tick_params(labelsize=6)

# ~ print('Path out to save figure')
# ~ # Path out to save figure
# ~ path_out = '/home/nice/Documentos/FPS_SESA/figs'
# ~ name_out = 'pyplt_stations_cluster_analysis.png'
# ~ plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
# ~ plt.close('all')

plt.show()
exit()
