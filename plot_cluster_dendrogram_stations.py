# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/02/2023"
__description__ = "This script plot cluster analysis from each weather station"

import os
import conda
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn
from dtaidistance import dtw, clustering
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering


def import_inmet(dt, type_cycle):
		
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
			
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])

		print('Reading INMET weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet weather station	
		if type_cycle == 'diurnal_cycle':
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
			d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
			d_i = d_i.groupby('time.hour').mean('time')
			values_i = d_i.values
			clim_i.append(values_i*24)
		else:
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
			d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			clim_i.append(values_i*24)
			
	return iy, ix, clim_i


def import_urug_smn(dt, type_cycle):
	
	jy = []
	jx = []
	clim_j = []
	
	# Select lat and lon 
	for j in range(1, 72):
		
		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])		

		print('Reading Uruguai weather station:', j, urug_smn[j][0])	
		# Reading Uruguai weather stations
		if type_cycle == 'diurnal_cycle':
			d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
			d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
			d_j = d_j.groupby('time.hour').mean('time')
			values_j = d_j.values
			clim_j.append(values_j*24)
		else:
			d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
			d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
			d_j = d_j.groupby('time.month').mean('time')
			values_j = d_j.values
			clim_j.append(values_j*24)
							
	return jy, jx, clim_j


def import_arg_emas(dt, type_cycle):
	
	ky = []
	kx = []
	clim_k = []
	
	# Select lat and lon 
	for k in range(1, 88):

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
			
		ky.append(arg_emas[k][2])
		kx.append(arg_emas[k][1])		

		print('Reading Argentina weather station:', k, arg_emas[k][0])	
		# Reading Argentina weather stations
		if type_cycle == 'diurnal_cycle':
			d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
			d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
			d_k = d_k.groupby('time.hour').mean('time')
			values_k = d_k.values
			clim_k.append(values_k*24)
		else:
			d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
			d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
			d_k = d_k.groupby('time.month').mean('time')
			values_k = d_k.values
			clim_k.append(values_k*24)
				
	return ky, kx, clim_k
	

type_cycle = 'annual_cycle'	
dt = 'H_2018-01-01_2021-12-31'

print('Import latitude, longitude and database')
# Import latitude, longitude and database
iy, ix, clim_i = import_inmet(dt, type_cycle)			
jy, jx, clim_j = import_urug_smn(dt, type_cycle)
ky, kx, clim_k = import_arg_emas(dt, type_cycle)
 
lon_xx = ix+jx+kx
lat_yy = iy+jy+ky
clim_tot = clim_i+clim_j+clim_k
df = pd.DataFrame(clim_tot)

print(lon_xx)
print()
print(lat_yy)

print('Calculate cluster analisis')
# Linkage hierarchical 
z = linkage(df, method='ward', metric='euclidean')

# Agglomerative clustering
Agg_hc = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
y_hc = Agg_hc.fit_predict(df)
print(y_hc)

print('Plot figure')
# Plot figure   
plt.figure(figsize=(30,10))
dendrogram(z, leaf_rotation=90., leaf_font_size=6., color_threshold=3800)
plt.title('Dendrogram', fontsize=8) 
plt.xlabel('Weather automatic stations', fontsize=8) 
plt.ylabel('Euclidean distances', fontsize=8) 

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_stations_dendrogram_{0}.png'.format(type_cycle)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()
