# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot correlation maps from regcm and database"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet


def import_dataset(var):
	
	ix = []		  
	iy = []
	
	djf_corr_i = []
	mam_corr_i = []
	jja_corr_i = []
	son_corr_i = []
	
	djf_corr_ii = []
	mam_corr_ii = []
	jja_corr_ii = []
	son_corr_ii = []

	djf_corr_iii = []
	mam_corr_iii = []
	jja_corr_iii = []
	son_corr_iii = []

	# Select lat and lon 
	for i in range(1, 155):
			
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		
		if var == 't2m':
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2019-01-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			djf_d_i = d_i.sel(time=d_i.time.dt.season=='DJF')
			mam_d_i = d_i.sel(time=d_i.time.dt.season=='MAM')
			jja_d_i = d_i.sel(time=d_i.time.dt.season=='JJA')
			son_d_i = d_i.sel(time=d_i.time.dt.season=='SON')
			djf_d_i = djf_d_i.groupby(djf_d_i.time.dt.year).mean("time")
			mam_d_i = mam_d_i.groupby(mam_d_i.time.dt.year).mean("time")
			jja_d_i = jja_d_i.groupby(jja_d_i.time.dt.year).mean("time")
			son_d_i = son_d_i.groupby(son_d_i.time.dt.year).mean("time")
			djf_list_i = djf_d_i.values-273.15
			mam_list_i = mam_d_i.values-273.15
			jja_list_i = jja_d_i.values-273.15
			son_list_i = son_d_i.values-273.15
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_ii = d_ii.tmp.sel(time=slice('2019-01-01','2021-12-31'))
			djf_d_ii = d_ii.sel(time=d_ii.time.dt.season=='DJF')
			mam_d_ii = d_ii.sel(time=d_ii.time.dt.season=='MAM')
			jja_d_ii = d_ii.sel(time=d_ii.time.dt.season=='JJA')
			son_d_ii = d_ii.sel(time=d_ii.time.dt.season=='SON')
			djf_d_ii = djf_d_ii.groupby(djf_d_ii.time.dt.year).mean("time")
			mam_d_ii = mam_d_ii.groupby(mam_d_ii.time.dt.year).mean("time")
			jja_d_ii = jja_d_ii.groupby(jja_d_ii.time.dt.year).mean("time")
			son_d_ii = son_d_ii.groupby(son_d_ii.time.dt.year).mean("time")
			djf_list_ii = djf_d_ii.values
			mam_list_ii = mam_d_ii.values
			jja_list_ii = jja_d_ii.values
			son_list_ii = son_d_ii.values	
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.t2m.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			djf_d_iii = d_iii.sel(time=d_iii.time.dt.season=='DJF')
			mam_d_iii = d_iii.sel(time=d_iii.time.dt.season=='MAM')
			jja_d_iii = d_iii.sel(time=d_iii.time.dt.season=='JJA')
			son_d_iii = d_iii.sel(time=d_iii.time.dt.season=='SON')
			djf_d_iii = djf_d_iii.groupby(djf_d_iii.time.dt.year).mean("time")
			mam_d_iii = mam_d_iii.groupby(mam_d_iii.time.dt.year).mean("time")
			jja_d_iii = jja_d_iii.groupby(jja_d_iii.time.dt.year).mean("time")
			son_d_iii = son_d_iii.groupby(son_d_iii.time.dt.year).mean("time")
			djf_list_iii = djf_d_iii.values-273.15
			mam_list_iii = mam_d_iii.values-273.15
			jja_list_iii = jja_d_iii.values-273.15
			son_list_iii = son_d_iii.values-273.15
		else:
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2019-01-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			djf_d_i = d_i.sel(time=d_i.time.dt.season=='DJF')
			mam_d_i = d_i.sel(time=d_i.time.dt.season=='MAM')
			jja_d_i = d_i.sel(time=d_i.time.dt.season=='JJA')
			son_d_i = d_i.sel(time=d_i.time.dt.season=='SON')
			djf_d_i = djf_d_i.groupby(djf_d_i.time.dt.year).mean("time")
			mam_d_i = mam_d_i.groupby(mam_d_i.time.dt.year).mean("time")
			jja_d_i = jja_d_i.groupby(jja_d_i.time.dt.year).mean("time")
			son_d_i = son_d_i.groupby(son_d_i.time.dt.year).mean("time")
			djf_list_i = djf_d_i.values-273.15
			mam_list_i = mam_d_i.values-273.15
			jja_list_i = jja_d_i.values-273.15
			son_list_i = son_d_i.values-273.15
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_ii = d_ii.uv.sel(time=slice('2019-01-01','2021-12-31'))
			djf_d_ii = d_ii.sel(time=d_ii.time.dt.season=='DJF')
			mam_d_ii = d_ii.sel(time=d_ii.time.dt.season=='MAM')
			jja_d_ii = d_ii.sel(time=d_ii.time.dt.season=='JJA')
			son_d_ii = d_ii.sel(time=d_ii.time.dt.season=='SON')
			djf_d_ii = djf_d_ii.groupby(djf_d_ii.time.dt.year).mean("time")
			mam_d_ii = mam_d_ii.groupby(mam_d_ii.time.dt.year).mean("time")
			jja_d_ii = jja_d_ii.groupby(jja_d_ii.time.dt.year).mean("time")
			son_d_ii = son_d_ii.groupby(son_d_ii.time.dt.year).mean("time")
			djf_list_ii = djf_d_ii.values
			mam_list_ii = mam_d_ii.values
			jja_list_ii = jja_d_ii.values
			son_list_ii = son_d_ii.values
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.u10.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			djf_d_iii = d_iii.sel(time=d_iii.time.dt.season=='DJF')
			mam_d_iii = d_iii.sel(time=d_iii.time.dt.season=='MAM')
			jja_d_iii = d_iii.sel(time=d_iii.time.dt.season=='JJA')
			son_d_iii = d_iii.sel(time=d_iii.time.dt.season=='SON')
			djf_d_iii = djf_d_iii.groupby(djf_d_iii.time.dt.year).mean("time")
			mam_d_iii = mam_d_iii.groupby(mam_d_iii.time.dt.year).mean("time")
			jja_d_iii = jja_d_iii.groupby(jja_d_iii.time.dt.year).mean("time")
			son_d_iii = son_d_iii.groupby(son_d_iii.time.dt.year).mean("time")
			djf_list_iii = djf_d_iii.values-273.15
			mam_list_iii = mam_d_iii.values-273.15
			jja_list_iii = jja_d_iii.values-273.15
			son_list_iii = son_d_iii.values-273.15

		# calculate correlation
		djf_corr_i.append(np.corrcoef(djf_list_ii, djf_list_i)[0][1])
		mam_corr_i.append(np.corrcoef(mam_list_ii, mam_list_i)[0][1])
		jja_corr_i.append(np.corrcoef(jja_list_ii, jja_list_i)[0][1])
		son_corr_i.append(np.corrcoef(son_list_ii, son_list_i)[0][1])

		djf_corr_ii.append(np.corrcoef(djf_list_iii, djf_list_i)[0][1])
		mam_corr_ii.append(np.corrcoef(mam_list_iii, mam_list_i)[0][1])
		jja_corr_ii.append(np.corrcoef(jja_list_iii, jja_list_i)[0][1])
		son_corr_ii.append(np.corrcoef(son_list_iii, son_list_i)[0][1])
	
	return iy, ix, djf_corr_i, mam_corr_i, jja_corr_i, son_corr_i, djf_corr_ii, mam_corr_ii, jja_corr_ii, son_corr_ii
		

def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-75., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-75.,-48.,6.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)

	return my_map


var = 't2m'

print('Import dataset')
# Import dataset
iy, ix, djf_corr_i, mam_corr_i, jja_corr_i, son_corr_i, djf_corr_ii, mam_corr_ii, jja_corr_ii, son_corr_ii = import_dataset(var)			

lon_xx = ix
lat_yy = iy

regcm_inmet_djf = djf_corr_i
regcm_inmet_mam = mam_corr_i
regcm_inmet_jja = jja_corr_i
regcm_inmet_son = son_corr_i

regcm_era5_djf = djf_corr_ii
regcm_era5_mam = mam_corr_ii
regcm_era5_jja = jja_corr_ii
regcm_era5_son = son_corr_ii

print('Plot figure')
# Plot figure   
fig = plt.figure(figsize=(5, 7))

if var == 't2m':
	color='PRGn'
	v_min = -1
	v_max = 1
	legend = 'Correlation of temperature'
else:
	color='PRGn'
	v_min = -1
	v_max = 1
	legend = 'Correlation wind 10m'

ax = fig.add_subplot(4, 2, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_djf, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM47 (INMET) DJF', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')

ax = fig.add_subplot(4, 2, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_djf, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) RegCM47 (ERA5) DJF', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 2, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_djf, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) RegCM47 (INMET) MAM', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')

ax = fig.add_subplot(4, 2, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_mam, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) RegCM47 (ERA5) MAM', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 2, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_mam, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) RegCM47 (INMET) JJA', loc='left', fontsize=6, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')

ax = fig.add_subplot(4, 2, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_mam, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) RegCM47 (ERA5) JJA', loc='left', fontsize=6, fontweight='bold')

ax = fig.add_subplot(4, 2, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_inmet_jja, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) RegCM47 (INMET) SON', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=6, labelpad=20, fontweight='bold')

ax = fig.add_subplot(4, 2, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, regcm_era5_jja, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(h) RegCM47 (ERA5) SON', loc='left', fontsize=6, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=6, labelpad=10, fontweight='bold')

cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
cbar = fig.colorbar(pltfig, cax=cb_ax, orientation='vertical', shrink=0.5, pad=0.5, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=8)  
	
print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_maps_corr_{0}_sesa.png'.format(var)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()

