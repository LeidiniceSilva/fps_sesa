# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/01/2023"
__description__ = "This script plot climatology maps from regcm5 and database"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_csam_inmet_stations import inmet


def import_dataset():
	
	ix = []		  
	iy = []
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []

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
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])

		# reading regcm 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-12-31'))
		d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_i = d_i.groupby('time.season').mean('time')
		values_i = d_i.values
		mean_i.append(values_i*86400)

		# Reading inmet 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2018-06-01','2021-12-31'))
		d_ii = d_ii.groupby('time.season').mean('time')
		values_ii = d_ii.values
		mean_ii.append(values_ii*24)
				
		# reading cmorph 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_4km_mon_20180101-20211231.nc')
		d_iii = d_iii.cmorph.sel(time=slice('2018-06-01','2021-12-31'))
		d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iii = d_iii.groupby('time.season').mean('time')
		values_iii = d_iii.values
		mean_iii.append(values_iii)

		# reading era5 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'mtpr_era5_csam_4km_mon_20180101-20211231.nc')
		d_iv = d_iv.mtpr.sel(time=slice('2018-06-01','2021-12-31'))
		d_iv = d_iv.groupby('time.season').mean('time')
		d_iv = d_iv.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		values_iv = d_iv.values
		mean_iv.append(values_iv*86400)
		
	return iy, ix, mean_i, mean_ii, mean_iii, mean_iv
		

def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-75., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-75.,-48.,6.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)

	return my_map
	
print('Import dataset')
# Import dataset
iy, ix, mean_i, mean_ii, mean_iii, mean_iv = import_dataset()			

lon_xx = ix
lat_yy = iy

regcm_djf = [item[0] for item in mean_i]
regcm_mam = [item[1] for item in mean_i]
regcm_jja = [item[2] for item in mean_i]
regcm_son = [item[3] for item in mean_i]

inmet_djf = [item[0] for item in mean_ii]
inmet_mam = [item[1] for item in mean_ii]
inmet_jja = [item[2] for item in mean_ii]
inmet_son = [item[3] for item in mean_ii]

cmorph_djf = [item[0] for item in mean_iii]
cmorph_mam = [item[1] for item in mean_iii]
cmorph_jja = [item[2] for item in mean_iii]
cmorph_son = [item[3] for item in mean_iii]

era5_djf = [item[0] for item in mean_iv]
era5_mam = [item[1] for item in mean_iv]
era5_jja = [item[2] for item in mean_iv]
era5_son = [item[3] for item in mean_iv]

print('Plot figure')
# Plot figure   
fig = plt.figure()

ax = fig.add_subplot(4, 4, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_djf, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(a) RegCM47 DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, inmet_djf, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(b) INMET DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, cmorph_djf, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(c) CMORPH DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, era5_djf, cmap='Blues', marker='o', vmin=0, vmax=120)
plt.title('(d) ERA5 DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_mam, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(e) RegCM47 MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, inmet_mam, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(f) INMET MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, cmorph_mam, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(g) CMORPH MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, era5_mam, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(h) ERA5 MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 9)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_jja, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(i) RegCM47 JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 10)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, inmet_jja, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(j) INMET JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 11)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, cmorph_jja, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(k) CMORPH JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 12)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, era5_jja, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(l) ERA5 JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 13)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_son, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(m) RegCM47 SON', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 14)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, inmet_son, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(n) INMET SON', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 15)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, cmorph_son, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(o) CMORPH SON', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 4, 16)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, era5_son, cmap='Blues', marker='o', vmin=0, vmax=12)
plt.title('(p) ERA SON', loc='left', fontsize=6, fontweight='bold')
	
cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
cbar = fig.colorbar(pltfig, cax=cb_ax, orientation='vertical', shrink=0.5, pad=0.5, extend='max')
cbar.set_label('Climatology of precipitation (mm d⁻¹)', fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=8)  
	
print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/csam'
name_out = 'pyplt_maps_clim_pr_csam.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()

