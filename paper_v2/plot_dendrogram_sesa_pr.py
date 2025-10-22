# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot dendrogram of the weather stations"

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

path = '/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA'


def import_inmet():
		
	ix = []		  
	iy = []
	clim_i = []

	# Select lat and lon 
	for i in range(1, 100):
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet 
		d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
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

clim_tot = clim_i_x + clim_i_y + clim_i_z
df = pd.DataFrame(clim_tot)
dp = df.dropna(axis=0)

print(dp)
print(len(dp))
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
Agg_hc = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
y_hc = Agg_hc.fit_predict(df)
print(y_hc)
print(len(y_hc))

# Plot figure   
plt.figure(figsize=(30,10))
dendrogram(z, leaf_rotation=90., leaf_font_size=10., color_threshold=3800)
plt.title('Dendrogram', fontsize=20) 
plt.xlabel('Weather stations (INMET+SMN)', fontsize=20) 
plt.ylabel('Euclidean distances', fontsize=20) 

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_dendrogram_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
