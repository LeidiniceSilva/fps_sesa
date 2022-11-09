# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot point for automatic station over sesa domain"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn


def import_inmet_era5(var_i, var_ii, dt):
	
	ix = []		  
	iy = []
	bias_ii = []
	corr_ii = []

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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + '{0}_{1}_{2}.nc'.format(var_i, inmet[i][0], dt))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i = values_i*24

		# reading era5 reanalisis
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + '{0}_sesa_era5_2018-2021.nc'.format(var_ii))
		d_ii = d_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_ii = d_ii.sel(latitude=inmet[i][2], longitude=inmet[i][3], method='nearest')
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		clim_ii = values_ii*24
			
		corr_i = np.corrcoef(clim_i, clim_ii)[0][1]
		corr_ii.append(corr_i)

		bias_i = np.nanmean(clim_ii) - np.nanmean(clim_i)
		bias_ii.append(bias_i)

	return iy, ix, corr_ii, bias_ii


def import_urug_smn_era5(dt):
	
	jy = []
	jx = []
	corr_jj = []
	bias_jj = []
	
	# Select lat and lon 
	for j in range(1, 72):
		
		print('Reading Uruguai weather station:', j, urug_smn[j][0])	

		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])		

		# Reading Uruguai weather stations
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_j = d_j.groupby('time.month').mean('time')
		values_j = d_j.values
		clim_j = values_j*24
		
		# Reading ERA5 reanalisis
		d_jj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
		d_jj = d_jj.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_jj = d_jj.sel(latitude=urug_smn[j][1], longitude=urug_smn[j][2], method='nearest')
		d_jj = d_jj.groupby('time.month').mean('time')
		values_jj = d_jj.values
		clim_jj = values_jj*24
		
		corr_j = np.corrcoef(clim_j, clim_jj)[0][1]
		corr_jj.append(corr_j)

		bias_j = np.nanmean(clim_jj) - np.nanmean(clim_j)
		bias_jj.append(bias_j)
				
	return jy, jx, corr_jj, bias_jj


def import_arg_emas_era5(dt):
	
	ky = []
	kx = []
	corr_kk = []
	bias_kk = []
	
	# Select lat and lon 
	for k in range(1, 88):
		
		print('Reading Argentina weather station:', k, arg_emas[k][0])	

		ky.append(arg_emas[k][2])
		kx.append(arg_emas[k][1])		

		# Reading Argentina weather stations
		d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
		d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
		d_k = d_k.groupby('time.month').mean('time')
		values_k = d_k.values
		clim_k = values_k*24
		
		# Reading ERA5 reanalisis
		d_kk = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
		d_kk = d_kk.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_kk = d_kk.sel(latitude=arg_emas[k][2], longitude=arg_emas[k][1], method='nearest')
		d_kk = d_kk.groupby('time.month').mean('time')
		values_kk = d_kk.values
		clim_kk = values_kk*24
		
		corr_k = np.corrcoef(clim_k, clim_kk)[0][1]
		corr_kk.append(corr_k)

		bias_k = np.nanmean(clim_kk) - np.nanmean(clim_k)
		bias_kk.append(bias_k)
				
	return ky, kx, corr_kk, bias_kk
	
	
dict_var = {1: ['pre', 'Pr_1h (mm d$\mathregular{^{-1}}$)', 'tp'],
			5: ['tmp', 'Tmp_1h (°C)', 't2m'],
			11: ['uv', 'Wind_1h (m s$\mathregular{^{-1}}$)', 'uv10']}
			
idx=1
var_i = dict_var[idx][0]
var_ii = dict_var[idx][2]
dt = 'H_2018-01-01_2021-12-31'

# Import latitude, longitude, correlation and bias
iy, ix, corr_ii, bias_ii = import_inmet_era5(var_i, var_ii, dt)			
jy, jx, corr_jj, bias_jj = import_urug_smn_era5(dt)
ky, kx, corr_kk, bias_kk = import_arg_emas_era5(dt)

if idx == 1:
	lon_xx = ix+jx+kx
	lat_yy = iy+jy+ky
	corr_tot = corr_ii+corr_jj+corr_kk
	bias_tot = bias_ii+bias_jj+bias_kk
else:
	lon_xx = ix
	lat_yy = iy
	corr_tot = corr_ii 
	bias_tot = bias_ii

print('Plot figure')
# Plot figure   
fig = plt.figure()

if idx == 1:
	minb=-2
	maxb=2
	minc=-1
	maxc=1
	colorv='BrBG'
elif idx == 5:
	minb=-2
	maxb=2
	minc=-1
	maxc=1
	colorv='bwr'
else:
	minb=-2
	maxb=2
	minc=-1
	maxc=1
	colorv='PiYG'
	
ax = fig.add_subplot(1, 2, 1)
my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
pltfig = my_map.scatter(lon_xx, lat_yy, 5, bias_tot, cmap=colorv, marker='o', vmin=minb, vmax=maxb)
plt.title('(a) Bias {0}'.format(dict_var[idx][1]), loc='left', fontsize=8)
plt.ylabel(u'Latitude', fontsize=6, labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.35, 0.019, 0.28]), extend='both')
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(1, 2, 2)
my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
pltfig = my_map.scatter(lon_xx, lat_yy, 5, corr_tot, cmap=colorv, marker='o', vmin=minc, vmax=maxc)
plt.title('(b) Correlation {0}'.format(dict_var[idx][1]), loc='left', fontsize=8)
plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.99, 0.35, 0.019, 0.28]), extend='both')
cbar.ax.tick_params(labelsize=6)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_maps_stats_{0}_weather_station.png'.format(dict_var[idx][0])
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
exit()

