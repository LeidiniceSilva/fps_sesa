# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot bias maps from regcm and database"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet
from dict_urug_smn_stations import urug_smn


def import_inmet(dt):
	
	ix = []		  
	iy = []
	bias_i = []
	bias_ii = []
	bias_iii = []

	# Select lat and lon 
	for i in range(1, 101):
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		# reading regcm 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2019-01-01','2021-12-31'))
		d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_i = d_i.groupby('time.year').mean('time')
		values_i = d_i.values
		list_i = np.nanmean(values_i*86400)
		# Reading inmet 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_ii = d_ii.pre.sel(time=slice('2019-01-01','2021-12-31'))
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = d_ii.values
		list_ii = np.nanmean(values_ii*24)
		# reading cmorph 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_4km_mon_20180101-20211231.nc')
		d_iii = d_iii.cmorph.sel(time=slice('2019-01-01','2021-12-31'))
		d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = d_iii.values
		list_iii = np.nanmean(values_iii)
		# reading era5 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_iv = d_iv.tp.sel(time=slice('2019-01-01','2021-12-31'))
		d_iv = d_iv.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = d_iv.values
		list_iv = np.nanmean(values_iv)

		# calculate bias
		mean_i = list_i - list_ii
		bias_i.append(mean_i)

		mean_ii = list_i - list_iii
		bias_ii.append(mean_ii)
		
		mean_iii = list_i - list_iv
		bias_iii.append(mean_iii)
				
	return iy, ix, bias_i, bias_ii, bias_iii
		

def import_urug_smn(dt):
	
	jx = []		  
	jy = []
	bias_j = []
	bias_jj = []
	bias_jjj = []

	# Select lat and lon 
	for j in range(1, 72):
		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])	
		
		print('Reading weather station:', j, urug_smn[j][0])	
		# reading regcm 
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_j = d_j.pr.sel(time=slice('2019-01-01','2021-12-31'))
		d_j = d_j.sel(lat=urug_smn[j][1], lon=urug_smn[j][2], method='nearest')
		d_j = d_j.groupby('time.year').mean('time')
		values_j = d_j.values
		list_j = np.nanmean(values_j*86400)
		# Reading smn
		d_jj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_jj = d_jj.pre.sel(time=slice('2019-01-01','2021-12-31'))
		d_jj = d_jj.groupby('time.year').mean('time')
		values_jj = d_jj.values
		list_jj = np.nanmean(values_jj*24)
		# reading cmorph 
		d_jjj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_4km_mon_20180101-20211231.nc')
		d_jjj = d_jjj.cmorph.sel(time=slice('2019-01-01','2021-12-31'))
		d_jjj = d_jjj.sel(lat=urug_smn[j][1], lon=urug_smn[j][2], method='nearest')
		d_jjj = d_jjj.groupby('time.year').mean('time')
		values_jjj = d_jjj.values
		list_jjj = np.nanmean(values_jjj)
		# reading era5 
		d_jv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_jv = d_jv.tp.sel(time=slice('2019-01-01','2021-12-31'))
		d_jv = d_jv.sel(lat=urug_smn[j][1], lon=urug_smn[j][2], method='nearest')
		d_jv = d_jv.groupby('time.year').mean('time')
		values_jv = d_jv.values
		list_jv = np.nanmean(values_jv)
		
		# calculate bias
		mean_j = list_j - list_jj
		bias_j.append(mean_j)

		mean_jj = list_j - list_jjj
		bias_jj.append(mean_jj)
		
		mean_jjj = list_j - list_jv
		bias_jjj.append(mean_jjj)
		
	return jy, jx, bias_j, bias_jj, bias_jjj
	
	
def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], size=6, linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)

	return my_map


var = 'pr'
dt = 'H_2018-01-01_2021-12-31'

print('Import dataset')
# Import dataset
iy, ix, bias_i, bias_ii, bias_iii = import_inmet(dt)			
jy, jx, bias_j, bias_jj, bias_jjj = import_urug_smn(dt)			

lon_xx = ix+jx
lat_yy = iy+jy

bias_tot_i = bias_i+bias_j
bias_tot_ii = bias_ii+bias_jj
bias_tot_iii = bias_iii+bias_jjj

print('Plot figure')
# Plot figure   
fig = plt.figure(figsize=(4, 5))

color='BrBG'
v_min = -3
v_max = 3
legend = 'Bias of precipitation (mm d⁻¹)'

ax = fig.add_subplot(2, 1, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, bias_tot_i, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM - INMET+SMN', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)  

ax = fig.add_subplot(2, 1, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, bias_tot_iii, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) RegCM - ERA5', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)  
	
print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

