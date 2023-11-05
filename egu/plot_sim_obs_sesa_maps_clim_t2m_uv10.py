# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot climatology maps from regcm5 and database"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet


def import_inmet(var, dt):
	
	ix = []		  
	iy = []
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []

	# Select lat and lon 
	for i in range(1, 101):
			
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		
		if var == 't2m':
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2019-06-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_i = d_i.groupby('time.year').mean('time')
			values_i = d_i.values
			values_i = np.nanmean(values_i)
			mean_i.append(values_i-273.15)
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'tmp_{0}_{1}.nc'.format(inmet[i][0], dt))
			d_ii = d_ii.tmp.sel(time=slice('2019-01-01','2021-12-31'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = d_ii.values
			values_ii = np.nanmean(values_ii)
			mean_ii.append(values_ii)		
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.t2m.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = d_iii.values
			values_iii = np.nanmean(values_iii)
			mean_iii.append(values_iii-273.15)
		else:
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2019-06-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_i = d_i.groupby('time.year').mean('time')
			values_i = d_i.values
			values_i = np.nanmean(values_i)
			mean_i.append(values_i)
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'uv_{0}_{1}.nc'.format(inmet[i][0], dt))
			d_ii = d_ii.uv.sel(time=slice('2019-01-01','2021-12-31'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = d_ii.values
			values_ii = np.nanmean(values_ii)
			mean_ii.append(values_ii)	
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.u10.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = d_iii.values
			values_iii = np.nanmean(values_iii)
			mean_iii.append(values_iii)
					
	return iy, ix, mean_i, mean_ii, mean_iii
		

def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], size=6, linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)

	return my_map


var = 'uv10'
dt = 'H_2018-01-01_2021-12-31'

print('Import dataset')
# Import dataset
iy, ix, mean_i, mean_ii, mean_iii = import_inmet(var, dt)			

lon_xx = ix
lat_yy = iy

regcm_tot = mean_i
inmet_tot = mean_ii
era5_tot = mean_iii

print('Plot figure')
# Plot figure   
fig = plt.figure(figsize=(4, 6))

if var == 't2m':
	color='Reds'
	v_min = 10
	v_max = 30
	legend = 'Temperature (°C)'
else:
	color='Greens'
	v_min = 0
	v_max = 8
	legend = 'Wind 10m (m s⁻¹)'
		
ax = fig.add_subplot(3, 1, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_tot, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='max')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)  

ax = fig.add_subplot(3, 1, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, inmet_tot, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) INMET', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='max')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)  

ax = fig.add_subplot(3, 1, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, era5_tot, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) ERA5', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='max')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)   
	
print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_maps_clim_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

