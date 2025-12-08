# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 16, 2023"
__description__ = "This script plot dendrogram of the weather stations"

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram

type_cycle = 'diurnal_cycle'	
path = '/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA'
# ~ skip_list = [1,2,11,30,36,43,53]

		
def import_inmet():
		
	ix = []		  
	iy = []
	clim_i = []
	# Select lat and lon 
	for i in range(1, 100):
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0])
		# Reading inmet 
		d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.hour').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
			
	return iy, ix, clim_i


def import_smn_i():
	
	ix = []		  
	iy = []
	clim_i = []
	# Select lat and lon 
	for i in range(1, 73):
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])
		# Reading smn 
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.hour').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
	
	return iy, ix, clim_i
	

# Import latitude, longitude and database
lat_x, lon_x, clim_i_x = import_inmet()			
lat_y, lon_y, clim_i_y = import_smn_i()			

lon_xx = lon_x + lon_y
lat_yy = lat_x + lat_y

clim_tot = clim_i_x + clim_i_y
df = pd.DataFrame(clim_tot)

print(df)
print()

print(lon_xx)
print(len(lon_xx))
print()
print(lat_yy)
print(len(lat_yy))
print()

# Linkage hierarchical 
z = linkage(df, method='ward', metric='euclidean')

# Agglomerative clustering
Agg_hc = AgglomerativeClustering(n_clusters=4, affinity='euclidean', linkage='ward')
y_hc = Agg_hc.fit_predict(df)
print(y_hc)
print(len(y_hc))

# Plot figure   
plt.figure(figsize=(30,10))
dendrogram(z, leaf_rotation=90., leaf_font_size=10., truncate_mode='level', color_threshold=3800)
plt.title('Dendrogram of precipitation diurnal cycle', fontsize=20) 
plt.xlabel('Weather stations (INMET+SMN)', fontsize=20) 
plt.ylabel('Euclidean distances', fontsize=20) 

# Path out to save figure
path_out = '{0}/figs/sesa_v2'.format(path)
name_out = 'pyplt_dendrogram_{0}_sesa.png'.format(type_cycle)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
