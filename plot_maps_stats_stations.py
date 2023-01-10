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


def import_inmet_era5(dt):
	
	ix = []		  
	iy = []
	mean_i = []
	mean_ii = []
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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		clim_i = values_i*24
		mean_i.append(np.nanmean(clim_i))

		# reading era5 reanalisis
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
		d_ii = d_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_ii = d_ii.sel(latitude=inmet[i][2], longitude=inmet[i][3], method='nearest')
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		clim_ii = values_ii*24
		mean_ii.append(np.nanmean(clim_ii))

		# calculate stats
		bias_i = np.nanmean(clim_ii) - np.nanmean(clim_i)
		bias_ii.append(bias_i)
			
		corr_i = np.corrcoef(clim_i, clim_ii)[0][1]
		corr_ii.append(corr_i)
		
	return iy, ix, mean_i, mean_ii, bias_ii, corr_ii


def import_urug_smn_era5(dt):
	
	jy = []
	jx = []
	mean_j = []
	mean_jj = []
	corr_jj = []
	bias_jj = []
	
	# Select lat and lon 
	for j in range(1, 72):
		
		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])		

		print('Reading Uruguai weather station:', j, urug_smn[j][0])	
		# Reading Uruguai weather stations
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_j = d_j.groupby('time.month').mean('time')
		values_j = d_j.values
		clim_j = values_j*24
		mean_j.append(np.nanmean(clim_j))

		# Reading ERA5 reanalisis
		d_jj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
		d_jj = d_jj.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_jj = d_jj.sel(latitude=urug_smn[j][1], longitude=urug_smn[j][2], method='nearest')
		d_jj = d_jj.groupby('time.month').mean('time')
		values_jj = d_jj.values
		clim_jj = values_jj*24
		mean_jj.append(np.nanmean(clim_jj))

		# calculate stats
		bias_j = np.nanmean(clim_jj) - np.nanmean(clim_j)
		bias_jj.append(bias_j)
						
		corr_j = np.corrcoef(clim_j, clim_jj)[0][1]
		corr_jj.append(corr_j)

	return jy, jx, mean_j, mean_jj, corr_jj, bias_jj


def import_arg_emas_era5(dt):
	
	ky = []
	kx = []
	mean_k = []
	mean_kk = []
	corr_kk = []
	bias_kk = []
	
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
		d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
		d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
		d_k = d_k.groupby('time.month').mean('time')
		values_k = d_k.values
		clim_k = values_k*24
		mean_k.append(np.nanmean(clim_k))
		
		# Reading ERA5 reanalisis
		d_kk = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
		d_kk = d_kk.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_kk = d_kk.sel(latitude=arg_emas[k][2], longitude=arg_emas[k][1], method='nearest')
		d_kk = d_kk.groupby('time.month').mean('time')
		values_kk = d_kk.values
		clim_kk = values_kk*24
		mean_kk.append(np.nanmean(clim_kk))

		# calculate stats
		bias_k = np.nanmean(clim_kk) - np.nanmean(clim_k)
		bias_kk.append(bias_k)
			
		corr_k = np.corrcoef(clim_k, clim_kk)[0][1]
		corr_kk.append(corr_k)
			
	return ky, kx, mean_k, mean_kk, corr_kk, bias_kk
	
	
dt = 'H_2018-01-01_2021-12-31'

print('Import latitude, longitude and database')
# Import latitude, longitude and database
iy, ix, mean_i, mean_ii, bias_ii, corr_ii = import_inmet_era5(dt)			
jy, jx, mean_j, mean_jj, bias_jj, corr_jj = import_urug_smn_era5(dt)
ky, kx, mean_k, mean_kk, bias_kk, corr_kk = import_arg_emas_era5(dt)

lon_xx = ix+jx+kx
lat_yy = iy+jy+ky
mean_tot_i = mean_i+mean_j+mean_k
mean_tot_ii = mean_ii+mean_jj+mean_kk
bias_tot = bias_ii+bias_jj+bias_kk
corr_tot = corr_ii+corr_jj+corr_kk

print('Plot figure')
# Plot figure   
fig = plt.figure()

ax = fig.add_subplot(2, 2, 1)
my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
pltfig = my_map.scatter(lon_xx, lat_yy, 5, mean_tot_i, cmap='Blues', marker='o', vmin=0, vmax=10)
plt.title('(a) Weather stations (mm d⁻¹)', loc='left', fontsize=8)
plt.ylabel(u'Latitude', fontsize=6, labelpad=15)

ax = fig.add_subplot(2, 2, 2)
my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
pltfig = my_map.scatter(lon_xx, lat_yy, 5, mean_tot_ii, cmap='Blues', marker='o', vmin=0, vmax=10)
plt.title('(b) ERA5 (mm d⁻¹)', loc='left', fontsize=8)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.55, 0.019, 0.28]), extend='max')
cbar.ax.tick_params(labelsize=6)
	
ax = fig.add_subplot(2, 2, 3)
my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
pltfig = my_map.scatter(lon_xx, lat_yy, 5, bias_tot, cmap='BrBG', marker='o', vmin=-2, vmax=2)
plt.title('(c) MBE (mm d⁻¹)', loc='left', fontsize=8)
plt.ylabel(u'Latitude', fontsize=6, labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.12, 0.019, 0.28]), extend='both')
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(2, 2, 4)
my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
pltfig = my_map.scatter(lon_xx, lat_yy, 5, corr_tot, cmap='PiYG', marker='o', vmin=-1, vmax=1)
plt.title('(d) PCC', loc='left', fontsize=8)
plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.99, 0.12, 0.019, 0.28]), extend='both')
cbar.ax.tick_params(labelsize=6)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_maps_stats_pre_weather_station.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()

