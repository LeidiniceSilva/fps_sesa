# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot maps of climatology bias"

import os
import cmocean
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_iii_stations import smn_iii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/afs/ictp.it/home/m/mda_silv/Documents'


def import_inmet():

	iy, ix = [], []
	bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 100):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		if var == 't2m':		
			# reading regcm usp 
			d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			values_i = np.nanmean(d_i.values)
			values_i = values_i-273.15
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = np.nanmean(d_ii.values)
			values_ii = values_ii-273.15
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = np.nanmean(d_iii.values)
			values_iii = values_iii-273.15
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = np.nanmean(d_iv.values)
			values_iv = values_iv-273.15
				
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			values_v = np.nanmean(d_v.values)
			values_v = values_v-273.15
				
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc_sesa/tmp/'.format(path) + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.tmp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			values_vi = np.nanmean(d_vi.values)
			values_vi = values_vi
			
			# reading era5 
			d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.t2m.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			values_vii = np.nanmean(d_vii.values)
			values_vii = values_vii-273.15
		else:
			# reading regcm usp 
			d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			values_i = np.nanmean(d_i.values)
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = np.nanmean(d_ii.values)
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = np.nanmean(d_iii.values)
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'sfcWind_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = np.nanmean(d_iv.values)
				
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			values_v = np.nanmean(d_v.values)
				
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc_sesa/uv/'.format(path) + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.uv.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			values_vi = np.nanmean(d_vi.values)
			
			# reading era5 
			d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.u10.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			values_vii = np.nanmean(d_vii.values)

		# calculate bias regcm usp
		bias_i.append(values_i - values_vi)	
		bias_ii.append(values_i - values_vii)

		# calculate bias regcm ictp 1
		bias_iii.append(values_ii - values_vi)	
		bias_iv.append(values_ii - values_vii)

		# calculate bias regcm ictp 2
		bias_v.append(values_iii - values_vi)	
		bias_vi.append(values_iii - values_vii)

		# calculate bias wrf ncar
		bias_vii.append(values_iv - values_vi)	
		bias_viii.append(values_iv - values_vii)

		# calculate bias wrf ucan
		bias_ix.append(values_v - values_vi)	
		bias_x.append(values_v - values_vii)

	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x


def import_smn_iii():
	
	iy, ix = [], []
	bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 102):
		yy=smn_iii[i][1]
		xx=smn_iii[i][2]
		iy.append(smn_iii[i][1])
		ix.append(smn_iii[i][2])
		
		print('Reading weather station:', i, smn_iii[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.year').mean('time')
		values_i = np.nanmean(d_i.values)
		values_i = values_i-273.15
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = np.nanmean(d_ii.values)
		values_ii = values_ii-273.15

		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = np.nanmean(d_iii.values)
		values_iii = values_iii-273.15
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.tas.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = np.nanmean(d_iv.values)
		values_iv = values_iv-273.15
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.tas.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.year').mean('time')
		values_v = np.nanmean(d_v.values)
		values_v = values_v-273.15
			
		# Reading smn 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_ii/smn_nc/tmp/'.format(path) + 'tmp_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_iii[i][0]))
		d_vi = d_vi.tmax.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		values_vi = np.nanmean(d_vi.values)
		values_vi = values_vi

		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.t2m.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.year').mean('time')
		values_vii = np.nanmean(d_vii.values)
		values_vii = values_vii-273.15

		# calculate bias regcm usp
		bias_i.append(values_i - values_vi)	
		bias_ii.append(values_i - values_vii)

		# calculate bias regcm ictp 1
		bias_iii.append(values_ii - values_vi)	
		bias_iv.append(values_ii - values_vii)

		# calculate bias regcm ictp 2
		bias_v.append(values_iii - values_vi)	
		bias_vi.append(values_iii - values_vii)

		# calculate bias wrf ncar
		bias_vii.append(values_iv - values_vi)	
		bias_viii.append(values_iv - values_vii)

		# calculate bias wrf ucan
		bias_ix.append(values_v - values_vi)	
		bias_x.append(values_v - values_vii)
			
	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x


def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
	my_map.drawmeridians(np.arange(-70,-45,5), labels=[0,0,0,1], size=font_size, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-40,-15,5), labels=[1,0,0,0], size=font_size, linewidth=0.5, color='black')
	my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map
	

var = 't2m'

# Import dataset
lat_x, lon_x, bias_i_x, bias_ii_x, bias_iii_x, bias_iv_x, bias_v_x, bias_vi_x, bias_vii_x, bias_viii_x, bias_ix_x, bias_x_x = import_inmet()			
lat_y, lon_y, bias_i_y, bias_ii_y, bias_iii_y, bias_iv_y, bias_v_y, bias_vi_y, bias_vii_y, bias_viii_y, bias_ix_y, bias_x_y = import_smn_iii()			

lat_yy = lat_x + lat_y 
lon_xx = lon_x + lon_y

reg_usp_inmet_smn = bias_i_x + bias_i_y 
reg_usp_reanalise = bias_ii_x + bias_ii_y 

reg_ictp_i_inmet_smn = bias_iii_x + bias_iii_y 
reg_ictp_i_reanalise = bias_iv_x + bias_iv_y 

reg_ictp_ii_inmet_smn = bias_v_x + bias_v_y 
reg_ictp_ii_reanalise = bias_vi_x + bias_vi_y 

wrf_ncar_inmet_smn = bias_vii_x + bias_vii_y 
wrf_ncar_reanalise = bias_viii_x + bias_viii_y 

wrf_ucan_inmet_smn = bias_ix_x + bias_ix_y 
wrf_ucan_reanalise = bias_x_x + bias_x_y 

# Plot figure   
fig = plt.figure(figsize=(11, 4.5))

if var == 't2m':
	color='bwr'
	v_min = -3
	v_max = 3
	legend = 'Bias of temperature (°C)'
	title='INMET+SMN'
else:
	color='PiYG'
	v_min = -3
	v_max = 3
	legend = 'Bias of wind 10m (m s⁻¹)'
	title='INMET'
font_size = 7
	
ax = fig.add_subplot(2, 5, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) Reg4 - {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.25, 0.010, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(2, 5, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) Reg5-Holt - {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) Reg5-UW - {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) WRF-NCAR - {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF-UCAN - {0}'.format(title), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) Reg4 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) Reg5-Holt - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(h) Reg5-UW - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 9)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(i) WRF-NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 10)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(j) WRF-UCAN - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '{0}/FPS_SESA/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


	


