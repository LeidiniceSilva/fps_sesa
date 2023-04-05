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
		d_i = d_i.groupby('time.season').mean('time')
		values_i = d_i.values
		list_i = values_i*86400

		# Reading inmet 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_ii = d_ii.pre.sel(time=slice('2019-01-01','2021-12-31'))
		d_ii = d_ii.groupby('time.season').mean('time')
		values_ii = d_ii.values
		list_ii = values_ii*24
				
		# reading cmorph 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_4km_mon_20180101-20211231.nc')
		d_iii = d_iii.cmorph.sel(time=slice('2019-01-01','2021-12-31'))
		d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iii = d_iii.groupby('time.season').mean('time')
		values_iii = d_iii.values
		list_iii = values_iii

		# reading era5 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_iv = d_iv.tp.sel(time=slice('2019-01-01','2021-12-31'))
		d_iv = d_iv.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iv = d_iv.groupby('time.season').mean('time')
		values_iv = d_iv.values
		list_iv = values_iv

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
		
		print('Reading Uruguai weather station:', j, urug_smn[j][0])	
		# reading regcm 
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_j = d_j.pr.sel(time=slice('2019-01-01','2021-12-31'))
		d_j = d_j.sel(lat=urug_smn[j][1], lon=urug_smn[j][2], method='nearest')
		d_j = d_j.groupby('time.season').mean('time')
		values_j = d_j.values
		list_j = values_j*86400

		# Reading smn
		d_jj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_jj = d_jj.pre.sel(time=slice('2019-01-01','2021-12-31'))
		d_jj = d_jj.groupby('time.season').mean('time')
		values_jj = d_jj.values
		list_jj = values_jj*24

		# reading cmorph 
		d_jjj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_4km_mon_20180101-20211231.nc')
		d_jjj = d_jjj.cmorph.sel(time=slice('2019-01-01','2021-12-31'))
		d_jjj = d_jjj.sel(lat=urug_smn[j][1], lon=urug_smn[j][2], method='nearest')
		d_jjj = d_jjj.groupby('time.season').mean('time')
		list_jjj = d_jjj.values

		# reading era5 
		d_jv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_jv = d_jv.tp.sel(time=slice('2019-01-01','2021-12-31'))
		d_jv = d_jv.groupby('time.season').mean('time')
		d_jv = d_jv.sel(lat=urug_smn[j][1], lon=urug_smn[j][2], method='nearest')
		list_jv = d_jv.values
		
		# calculate bias
		mean_j = list_j - list_jj
		bias_j.append(mean_j)

		mean_jj = list_j - list_jjj
		bias_jj.append(mean_jj)
		
		mean_jjj = list_j - list_jv
		bias_jjj.append(mean_jjj)
		
	return jy, jx, bias_j, bias_jj, bias_jjj
	
	
def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-75., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-75.,-48.,6.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
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

print(bias_tot_i)
print(bias_tot_ii)
print(bias_tot_iii)

print(len(lon_xx))
print(len(lat_yy))
print(len(bias_tot_iii))

regcm_inmet_djf = [item[0] for item in bias_tot_i]
regcm_inmet_mam = [item[1] for item in bias_tot_i]
regcm_inmet_jja = [item[2] for item in bias_tot_i]
regcm_inmet_son = [item[3] for item in bias_tot_i]

regcm_cmorph_djf = [item[0] for item in bias_tot_ii]
regcm_cmorph_mam = [item[1] for item in bias_tot_ii]
regcm_cmorph_jja = [item[2] for item in bias_tot_ii]
regcm_cmorph_son = [item[3] for item in bias_tot_ii]

regcm_era5_djf = [item[0] for item in bias_tot_iii]
regcm_era5_mam = [item[1] for item in bias_tot_iii]
regcm_era5_jja = [item[2] for item in bias_tot_iii]
regcm_era5_son = [item[3] for item in bias_tot_iii]

print('Plot figure')
# Plot figure   
fig = plt.figure(figsize=(6, 6))

color='BrBG'
v_min = -5
v_max = 5
legend = 'Bias of precipitation (mm d⁻¹)'

ax = fig.add_subplot(4, 3, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_djf, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM47 - INMET+SMN DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_cmorph_djf, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) RegCM47 - CMORPH DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_djf, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) RegCM47 - ERA5 DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_mam, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) RegCM47 - INMET+SMN MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_cmorph_mam, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) RegCM47 - CMORPH MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_mam, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) RegCM47 - ERA5 MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_jja, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) RegCM47 - INMET+SMN JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_cmorph_jja, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(h) RegCM47 - CMORPH JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 9)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_jja, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(i) RegCM47 - ERA5 JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 10)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_son, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(j) RegCM47 - INMET+SMN SON', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 11)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_cmorph_son, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(k) RegCM47 - CMORPh SON', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 3, 12)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_son, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(l) RegCM47 - ERA5 SON', loc='left', fontsize=6, fontweight='bold')

cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
cbar = fig.colorbar(pltfig, cax=cb_ax, orientation='vertical', shrink=0.5, pad=0.5, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=8)  
	
print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

