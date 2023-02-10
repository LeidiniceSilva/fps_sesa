# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot cluster analysis from each weather station"

import os
import conda
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_csam_inmet_stations import inmet
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering


def import_inmet(dt, type_cycle):
		
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for i in range(1, 155):
			
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


type_cycle = 'annual_cycle'	
dt = 'H_2018-01-01_2021-12-31'

print('Import latitude, longitude and database')
# Import latitude, longitude and database
iy, ix, clim_i = import_inmet(dt, type_cycle)			

lon_xx = ix
lat_yy = iy
clim_tot = clim_i
df = pd.DataFrame(clim_tot)

print(lon_xx)
print(len(lon_xx))
print()
print(lat_yy)
print(len(lat_yy))

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
dendrogram(z, leaf_rotation=90., leaf_font_size=7., color_threshold=3800)
plt.title('Dendrogram', fontsize=12) 
plt.xlabel('Weather automatic stations', fontsize=12) 
plt.ylabel('Euclidean distances', fontsize=12) 

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/csam'
name_out = 'pyplt_stations_dendrogram_{0}_csam.png'.format(type_cycle)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()
