# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot maps of correlation"

import os
import conda
import cmocean
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap


def import_inmet():

	iy, ix = [], []
	corr_i, corr_ii, corr_iii, corr_iv, corr_v, corr_vi, corr_vii, corr_viii, corr_ix, corr_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 100):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])		
		if var == 't2m':		
			# reading regcm usp 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			values_i = values_i-273.15
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
			d_ii = d_ii.groupby('time.month').mean('time')
			values_ii = d_ii.values
			values_ii = values_ii-273.15
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
			d_iii = d_iii.groupby('time.month').mean('time')
			values_iii = d_iii.values
			values_iii = values_iii-273.15
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
			d_iv = d_iv.groupby('time.month').mean('time')
			values_iv = d_iv.values
			values_iv = values_iv-273.15
				
			# reading wrf ucan 
			d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
			d_v = d_v.groupby('time.month').mean('time')
			values_v = d_v.values
			values_v = values_v-273.15
				
			# Reading inmet 
			d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc/tmp/' + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.tmp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.month').mean('time')
			values_vi = d_vi.values
			
			# reading era5 
			d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.t2m.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
			d_vii = d_vii.groupby('time.month').mean('time')
			values_vii = d_vii.values
			values_vii = values_vii-273.15
		else:
			# reading regcm usp 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
			d_ii = d_ii.groupby('time.month').mean('time')
			values_ii = d_ii.values
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
			d_iii = d_iii.groupby('time.month').mean('time')
			values_iii = d_iii.values
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'sfcWind_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
			d_iv = d_iv.groupby('time.month').mean('time')
			values_iv = d_iv.values
				
			# reading wrf ucan 
			d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
			d_v = d_v.groupby('time.month').mean('time')
			values_v = d_v.values
				
			# Reading inmet 
			d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc/uv/' + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.uv.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.month').mean('time')
			values_vi = d_vi.values
			
			# reading era5 
			d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.u10.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
			d_vii = d_vii.groupby('time.month').mean('time')
			values_vii = d_vii.values
		
		# calculate correlation regcm usp
		corr_i.append(np.corrcoef(values_i, values_vi)[0][1])
		corr_ii.append(np.corrcoef(values_i, values_vii)[0][1])

		# calculate correlation regcm ictp 1
		corr_iii.append(np.corrcoef(values_ii, values_vi)[0][1])
		corr_iv.append(np.corrcoef(values_ii, values_vii)[0][1])

		# calculate correlation regcm ictp 2
		corr_v.append(np.corrcoef(values_iii, values_vi)[0][1])	
		corr_vi.append(np.corrcoef(values_iii, values_vii)[0][1])

		# calculate correlation wrf ncar
		corr_vii.append(np.corrcoef(values_iv, values_vi)[0][1])
		corr_viii.append(np.corrcoef(values_iv, values_vii)[0][1])

		# calculate correlation wrf ucan
		corr_ix.append(np.corrcoef(values_v, values_vi)[0][1])
		corr_x.append(np.corrcoef(values_v, values_vii)[0][1])

	return iy, ix, corr_i, corr_ii, corr_iii, corr_iv, corr_v, corr_vi, corr_vii, corr_viii, corr_ix, corr_x


def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
	my_map.drawmeridians(np.arange(-70,-45,5), labels=[0,0,0,1], size=font_size, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-40,-15,5), labels=[1,0,0,0], size=font_size, linewidth=0.5, color='black')
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map
	

var = 't2m'

# Import dataset
lat_x, lon_x, corr_i_x, corr_ii_x, corr_iii_x, corr_iv_x, corr_v_x, corr_vi_x, corr_vii_x, corr_viii_x, corr_ix_x, corr_x_x = import_inmet()			

lat_yy = lat_x 
lon_xx = lon_x 

reg_usp_inmet_smn = corr_i_x
reg_usp_reanalise = corr_ii_x

reg_ictp_i_inmet_smn = corr_iii_x
reg_ictp_i_reanalise = corr_iv_x 

reg_ictp_ii_inmet_smn = corr_v_x 
reg_ictp_ii_reanalise = corr_vi_x 

wrf_ncar_inmet_smn = corr_vii_x 
wrf_ncar_reanalise = corr_viii_x 

wrf_ucan_inmet_smn = corr_ix_x 
wrf_ucan_reanalise = corr_x_x 

# Plot figure   
fig = plt.figure(figsize=(11, 4.5))

if var == 't2m':
	color='PRGn'
	v_min = -0.9
	v_max = 0.9
	legend = 'Correlation of temperature'
else:
	color='PRGn'
	v_min = -0.9
	v_max = 0.9
	legend = 'Correlation of wind 10m'
font_size = 7
	
ax = fig.add_subplot(2, 5, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) Reg4 vs. INMET', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.25, 0.010, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(2, 5, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) Reg5-Holt vs. INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) Reg5-UW vs. INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) WRF-NCAR vs. INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF-UCAN vs. INMET', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) Reg4 vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) Reg5-Holt vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(h) Reg5-UW vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 9)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(i) WRF-NCAR vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 10)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(j) WRF-UCAN vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_corr_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
