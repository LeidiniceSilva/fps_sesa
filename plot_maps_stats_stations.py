# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot point for each inmet automatic station over sesa domain"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_stations_inmet import inmet
from dict_stations_arg_ames import arg_ames
from dict_stations_urug_smn import urug_smn


def import_inmet_era5(var, dt):
	
	ix = []		  
	iy = []
	bias_dr = []
	corr_dr = []

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

		print('Reading inmet weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet weather station	
		dr_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + '{0}_{1}_{2}.nc'.format(var, inmet[i][0], dt))
		dr_i = dr_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		dr_i = dr_i.groupby('time.month').mean('time')
		values_dr_i = dr_i.values
		clim_dr_i = values_dr_i*24

		# reading era5 reanalisis
		dr_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + '{0}_sesa_era5_2018-2021.nc'.format(dict_var[idx][2]))
		dr_ii = dr_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
		dr_ii = dr_ii.sel(latitude=inmet[i][2], longitude=inmet[i][3], method='nearest')
		dr_ii = dr_ii.groupby('time.month').mean('time')
		values_dr_ii = dr_ii.values
		clim_dr_ii = values_dr_ii*24
			
		corr_ = np.corrcoef(clim_dr_i, clim_dr_ii)[0][1]
		corr_dr.append(corr_)

		bias_ = np.nanmean(clim_dr_ii) - np.nanmean(clim_dr_i)
		bias_dr.append(bias_)

	return iy, ix, corr_dr, bias_dr


def import_urug_smn_era5(dt):
	
	jy = []
	jx = []
	corr_ds = []
	bias_ds = []
	
	# Select lat and lon 
	for j in range(1, 72):
		
		print('Reading Uruguai weather station:', j, urug_smn[j][0])	

		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])		

		# Reading Uruguai weather stations
		ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		ds_i = ds_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		ds_i = ds_i.groupby('time.month').mean('time')
		values_ds_i = ds_i.values
		clim_ds_i = values_ds_i*24
		
		# Reading ERA5 reanalisis
		ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
		ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
		ds_ii = ds_ii.sel(latitude=urug_smn[j][1], longitude=urug_smn[j][2], method='nearest')
		ds_ii = ds_ii.groupby('time.month').mean('time')
		values_ds_ii = ds_ii.values
		clim_ds_ii = values_ds_ii*24
		
		corr_ = np.corrcoef(clim_ds_i, clim_ds_ii)[0][1]
		corr_ds.append(corr_)

		bias_ = np.nanmean(clim_ds_ii) - np.nanmean(clim_ds_i)
		bias_ds.append(bias_)
				
	return jy, jx, corr_ds, bias_ds
	
	
dict_var = {1: ['pre', 'Pr_1h (mm d$\mathregular{^{-1}}$)', 'tp'],
			5: ['tmp', 'Tmp_1h (°C)', 't2m'],
			11: ['uv', 'Wind_1h (m s$\mathregular{^{-1}}$)', 'uv10']}
idx=1
var = dict_var[idx][0]
dt = 'H_2018-01-01_2021-12-31'

# Import latitude, longitude, correlation and bias from Uruguai
iy, ix, corr_dr, bias_dr = import_inmet_era5(var, dt)			
jy, jx, corr_ds, bias_ds = import_urug_smn_era5(dt)

if idx == 1:
	lon_xx = ix+jx
	lat_yy = iy+jy
	corr_tot = corr_dr + corr_ds
	bias_tot = bias_dr + bias_ds
else:
	lon_xx = ix
	lat_yy = iy
	corr_tot = corr_dr 
	bias_tot = bias_dr

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

