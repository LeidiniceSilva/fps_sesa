# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot cluster analisis"

import os
import sys
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.path import Path
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
font_size = 8
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
	
	lat, lon, clim_i, clim_ii = [], [], [], []		  
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
			d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			clim_i.append(values_i)
		
			d_ii = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/pre/'.format(path) + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.groupby('time.year').sum('time')
			values_ii = d_ii.values
			clim_ii.append(values_ii)		
			
	return lat, lon, clim_i, clim_ii


def import_smn_i():
	
	lat, lon, clim_i, clim_ii = [], [], [], []		  
	for i in range(1, 73):
		lat.append(smn_i[i][1])
		lon.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])	
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i*24)
				
		d_ii = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.groupby('time.year').sum('time')
		values_ii = d_ii.values
		clim_ii.append(values_ii)	
		
	return lat, lon, clim_i, clim_ii
	

def import_smn_ii():
	
	lat, lon, clim_i, clim_ii = [], [], [], []		  
	for i in range(1, 110):
		if i in skip_list_smn_ii:
			continue
		lat.append(smn_ii[i][1])
		lon.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])	
		d_i = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i.append(values_i)
				
		d_ii = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.groupby('time.year').sum('time')
		values_ii = d_ii.values
		clim_ii.append(values_ii)	
		
	return lat, lon, clim_i, clim_ii
	
	
# Import latitude, longitude and database
lat_x, lon_x, clim_inmet_x, clim_inmet_xx = import_inmet()			
lat_y, lon_y, clim_smn_y, clim_smn_yy = import_smn_i()
lat_z, lon_z, clim_smn_z, clim_smn_zz = import_smn_ii()

lon_xx = lon_x + lon_y + lon_z
lat_yy = lat_x + lat_y + lat_z

clim_tot = clim_inmet_x + clim_smn_y + clim_smn_z
clim_tot_i = clim_inmet_xx + clim_smn_yy + clim_smn_zz

list_hc = [1, 2, 3, 2, 0, 1, 1, 0, 2, 2, 0, 3, 0, 2, 3, 0, 1, 2, 0, 3, 0, 4, 2, 4, 3, 1, 4, 2, 4, 2, 2, 2, 1, 2, 4, 2, 2, 3, 2, 4, 4, 4, 0, 2, 4, 3, 2, 0, 0, 0, 3, 2, 2, 2, 1, 2, 4, 1, 4, 3, 4, 3, 0, 2, 0, 3, 2, 3, 2, 4, 0, 1, 4, 2, 4, 4, 0, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 3, 2, 0, 0, 0, 4, 2, 3, 2, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 4, 0, 0, 4, 0, 4, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 2, 4, 3, 1, 4, 1, 2, 1, 1, 1, 4, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 2, 2, 2, 2, 2]
latitude = [-20.444444, -29.709167, -17.339444, -28.931353, -31.347778, -20.559167, -22.358052, -29.164581, -28.126944, -26.819156, -30.545278, -16.966944, -30.807953, -29.049167, -19.53921, -29.674293, -20.447195, -29.368889, -31.403333, -19.1225, -24.78, -23.359167, -26.417222, -25.322464, -19.98586, -18.996667, -18.492778, -28.60344, -25.448611, -27.288611, -26.286562, -31.0025, -21.457778, -27.657778, -27.6025, -27.395556, -26.398611, -16.642778, -28.653333, -25.567879, -26.913611, -26.950833, -23.981944, -26.081389, -23.449444, -18.409722, -27.418333, -25.010833, -32.534722, -23.773333, -17.923611, -27.169167, -28.222381, -25.371389, -21.666111, -28.604444, -24.533333, -22.235222, -23.405278, -17.454722, -25.508889, -17.745066, -23.000556, -27.920278, -22.658333, -16.9625, -28.226805, -17.304167, -25.721944, -22.5525, -30.053611, -21.338333, -22.12, -30.368578, -22.372778, -27.678611, -21.775, -32.078889, -26.248611, -29.872222, -29.725, -27.890556, -30.750556, -29.191599, -27.854444, -28.65, -30.341389, -28.275556, -28.748611, -28.417222, -25.835556, -29.702222, -28.704722, -20.981667, -16.679722, -28.859211, -33.742222, -29.350278, -30.010278, -21.927251, -29.089444, -18.916944, -29.84, -28.5325, -28.513611, -21.319167, -26.938666, -28.416682, -31.20645, -31.061733, -29.112957, -30.9499, -30.9982, -31.231098, -31.663054, -30.392535, -29.6036, -30.719082, -30.787252, -30.299302, -29.819517, -30.751196, -30.637259, -29.98349, -30.62427, -31.337833, -30.776241, -31.09242, -30.014261, -30.580333, -30.58639, -30.613837, -31.444407, -29.787069, -31.459, -30.755734, -27.298582, -30.971026, -30.971187, -30.992189, -28.177515, -31.037601, -31.286127, -30.432802, -29.05663, -30.590296, -30.865546, -29.660671, -29.845637, -30.690307, -30.214445, -30.626772, -30.249296, -30.787214, -30.289, -29.721745, -30.205281, -30.223953, -30.476397, -31.444635, -27.153304, -31.404212, -30.390258, -31.493837, -32.135, -31.350331, -31.38978, -31.272706, -30.337372, -27.869343, -29.296744, -30.91468, -28.545073, -30.988557, -31.004877, -29.377358, -30.431886, -28.691829, -29.470444, -34.57, -34.58, -28.6, -29.88, -31.3, -31.3, -31.4, -27.45, -34.6, -34.82, -26.2, -33.0, -25.73, -30.23, -24.38, -34.55, -34.97, -22.1, -29.38, -34.13, -24.7, -32.7, -32.83, -32.88, -30.27, -23.15, -31.78, -29.68, -31.67, -27.37, -29.18, -27.45, -33.12, -32.92, -24.85, -31.57, -33.27, -33.08, -34.58, -27.77, -31.7, -22.65, -28.07, -31.95, -33.73, -28.43, -34.61, -27.65, -27.42, -27.05, -33.73, -27.1, -32.68, -29.17, -31.85, -33.93, -31.18, -29.18, -26.87, -33.02, -30.38, -34.83, -34.45, -33.35, -32.37, -33.25, -32.35, -34.86, -30.9, -34.49, -31.4, -33.71, -34.09, -34.97, -34.35, -33.54, -32.69, -26.85, -28.02, -22.03, -22.28, -22.64, -23.5, -23.44, -24.09, -24.67, -24.03, -25.24, -25.75, -25.48, -26.88, -26.67, -26.18, -26.83, -27.3]
longitude = [-52.875833, -55.525556, -53.224444, -49.49792, -54.013333, -48.545, -49.028877, -51.534202, -49.479722, -50.98552, -53.466944, -51.8175, -51.83424, -50.149722, -49.518133, -51.064042, -54.722615, -50.827222, -52.700833, -51.720833, -50.046389, -52.931944, -52.348611, -49.157733, -48.151574, -57.6375, -53.171389, -53.673597, -49.230556, -50.604167, -53.633114, -54.618056, -51.552222, -52.305833, -48.62, -53.429444, -51.353611, -49.220278, -53.111944, -51.077946, -49.268056, -48.761944, -48.885833, -48.641667, -54.181944, -49.191944, -49.646944, -50.853889, -53.375833, -50.180556, -51.7175, -51.558889, -51.512845, -52.400833, -49.734722, -48.813333, -54.019167, -49.965111, -51.932778, -52.601111, -48.808611, -49.101698, -49.843333, -53.318056, -52.134444, -50.425556, -52.403582, -48.284167, -53.748056, -55.716389, -51.174722, -48.113889, -51.408611, -56.437115, -50.974722, -49.041944, -54.528056, -52.167778, -49.574167, -52.381944, -53.720556, -54.48, -55.401389, -54.885653, -53.791111, -56.016389, -54.310833, -49.934722, -50.057778, -54.9625, -50.368889, -54.694444, -51.870833, -54.971944, -48.618056, -52.542387, -53.372222, -49.733333, -50.135833, -50.490251, -53.826667, -48.255556, -57.081944, -49.315278, -50.882778, -50.930278, -52.39809, -56.542211, -56.224017, -55.866233, -56.555522, -57.4559, -57.485838, -57.097101, -56.618974, -56.456119, -58.165178, -57.325281, -57.78309, -56.972142, -57.429814, -56.214658, -56.326675, -58.295964, -56.657495, -56.629005, -58.141889, -57.021887, -57.85945, -57.680318, -56.230152, -56.905948, -57.421509, -58.070281, -57.903, -57.055016, -54.193437, -58.247317, -57.880716, -57.913847, -55.642269, -56.635879, -57.703535, -56.787211, -57.08769, -58.468999, -56.411658, -57.743836, -57.674446, -57.820722, -57.966907, -57.982836, -57.622519, -56.786308, -57.291694, -57.082856, -57.074737, -56.66754, -57.124992, -56.84021, -53.933283, -58.002561, -57.880984, -57.120524, -57.939, -56.392531, -57.977472, -57.916594, -58.317145, -55.129216, -57.565952, -57.92773, -56.028995, -56.215761, -56.879045, -58.193534, -57.440628, -57.480465, -56.81167, -58.42, -58.48, -65.77, -61.95, -58.02, -64.2, -64.18, -58.77, -58.6, -58.53, -58.23, -58.62, -54.47, -68.75, -65.08, -60.92, -57.9, -65.6, -66.82, -63.37, -60.58, -62.15, -68.78, -68.85, -57.65, -64.32, -60.48, -57.15, -63.88, -55.97, -59.7, -59.05, -64.23, -60.78, -65.48, -68.42, -66.35, -68.42, -68.4, -64.3, -60.82, -63.81, -67.57, -65.13, -65.38, -58.92, -58.67, -55.43, -58.93, -65.42, -69.12, -61.1, -62.12, -58.02, -60.54, -60.55, -61.55, -59.7, -60.45, -60.88, -56.51, -56.01, -57.77, -56.5, -54.19, -58.07, -58.04, -56.21, -55.54, -54.31, -57.97, -54.39, -56.19, -54.95, -56.76, -56.92, -57.65, -65.1, -64.23, -60.62, -57.94, -55.83, -58.79, -57.43, -57.09, -56.46, -54.35, -57.51, -56.44, -56.38, -58.32, -57.13, -56.35, -55.33, -55.89]

print(len(clim_tot))
print(len(clim_tot_i))
print(len(longitude))
print(len(latitude))
print(len(list_hc))

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

lon_c_i = []
lon_c_ii = []
lon_c_iii = []
lon_c_iv = []
lon_c_v = []

lat_c_i = []
lat_c_ii = []
lat_c_iii = []
lat_c_iv = []
lat_c_v = []

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
	if c_i < len(clim_tot):
		lon_c_i.append(longitude[c_i])
		lat_c_i.append(latitude[c_i])
		pre_c_i.append(clim_tot[c_i])
		pre_c_i_i.append(clim_tot_i[c_i])	

for c_ii in count_ii:
	if c_ii < len(clim_tot):
		lon_c_ii.append(longitude[c_ii])
		lat_c_ii.append(latitude[c_ii])
		pre_c_ii.append(clim_tot[c_ii])
		pre_c_ii_i.append(clim_tot_i[c_ii])

for c_iii in count_iii:
	if c_iii < len(clim_tot):
		lon_c_iii.append(longitude[c_iii])
		lat_c_iii.append(latitude[c_iii])
		pre_c_iii.append(clim_tot[c_iii])
		pre_c_iii_i.append(clim_tot_i[c_iii])
	
for c_iv in count_iv:
	if c_iv < len(clim_tot):
		lon_c_iv.append(longitude[c_iv])
		lat_c_iv.append(latitude[c_iv])
		pre_c_iv.append(clim_tot[c_iv])
		pre_c_iv_i.append(clim_tot_i[c_iv])
		
for c_v in count_v:
	if c_v < len(clim_tot):
		lon_c_v.append(longitude[c_v])
		lat_c_v.append(latitude[c_v])
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

box_plot_data = [cluster_i_i, cluster_ii_i, cluster_iii_i, cluster_iv_i, cluster_v_i]

# Plot figure   
fig = plt.figure(figsize=(14, 4))
lon_bounds = [-70, -46]
lat_bounds = [-40, -16]

ax1 = fig.add_subplot(1, 3, 1, projection=ccrs.PlateCarree())
ax1.plot(lon_c_i, lat_c_i, 'o', color='blue', label='Cluster I', markeredgecolor='black', markersize=6)
ax1.plot(lon_c_ii, lat_c_ii, 'o', color='red', label='Cluster II', markeredgecolor='black', markersize=6)
ax1.plot(lon_c_iii, lat_c_iii, 'o', color='green', label='Cluster III', markeredgecolor='black', markersize=6)
ax1.plot(lon_c_iv, lat_c_iv, 'o', color='gray', label='Cluster IV', markeredgecolor='black', markersize=6) 
ax1.plot(lon_c_v, lat_c_v, 'o', color='orange', label='Cluster V', markeredgecolor='black', markersize=6)
ax1.set_xlabel(u'Longitude', fontweight='bold')
ax1.set_ylabel(u'Latitude', fontweight='bold')
ax1.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
ax1.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 5), crs=ccrs.PlateCarree())
ax1.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 5), crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(LongitudeFormatter())
ax1.yaxis.set_major_formatter(LatitudeFormatter())
ax1.grid(c='k', ls='--', alpha=0.5)  
states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')
ax1.add_feature(states_provinces, edgecolor='0.1')
ax1.add_feature(cfeat.BORDERS, linewidth=0.75)
ax1.coastlines(linewidth=0.75)
ax1.set_title('(a)', loc='left', fontsize=10, fontweight='bold')
ax1.text(-68, -19, u'\u25B2 \nN', fontsize=10, fontweight='bold')
ax1.legend(loc=3, ncol=3, fontsize=8, frameon=False)

ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(np.arange(0.5, 12 + 0.5), cluster_i, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'Cluster I')
ax2.plot(np.arange(0.5, 12 + 0.5), cluster_ii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='red', label = 'Cluster II')
ax2.plot(np.arange(0.5, 12 + 0.5), cluster_iii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='green', label = 'Cluster III')
ax2.plot(np.arange(0.5, 12 + 0.5), cluster_iv, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='gray', label = 'Cluster IV')
ax2.plot(np.arange(0.5, 12 + 0.5), cluster_v, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='orange', label = 'Cluster V')
ax2.set_title('(b)', loc='left', fontsize=10, fontweight='bold')
ax2.set_ylabel('Precipitation (mm d⁻¹)', fontsize=10, fontweight='bold')
ax2.set_xlabel(u'Months', fontweight='bold')
ax2.set_xticks(np.arange(0.5, 12 + 0.5), ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
ax2.set_yticks(np.arange(0, 11, 1))
ax2.legend(loc=9, ncol=2, fontsize=10, handletextpad=0.2, frameon=False)

ax3 = fig.add_subplot(1, 3, 3)
box = ax3.boxplot(box_plot_data, patch_artist=True, labels=['Cluster I','Cluster II','Cluster III','Cluster IV', 'Cluster V'])

# Make each box colorful
colors = ['blue', 'red', 'green', 'gray', 'orange']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
for median in box['medians']:
    median.set(color='k', linewidth=1.)
ax3.set_title('(c)', loc='left', fontsize=10, fontweight='bold')
ax3.set_ylabel('Annual mean precipitation (mm)', fontsize=10, fontweight='bold')
ax3.set_xlabel(u'Clusters', fontweight='bold')
ax3.set_yticks(np.arange(200, 1800, 100))

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_cluster_analysis_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

