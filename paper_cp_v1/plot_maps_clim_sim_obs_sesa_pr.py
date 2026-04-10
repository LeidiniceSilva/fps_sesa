# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot maps of climatology"

import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.path import Path
from matplotlib.colors import BoundaryNorm
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
font_size = 8
path = '/home/mda_silv/users/FPS_SESA'

skip_list_inmet_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
skip_list_inmet_ii = [2, 3, 4, 14, 19, 20, 21, 24, 25, 26, 27, 28, 32, 33, 34, 35, 38, 40, 41, 44, 45, 48, 52, 54, 55, 56, 59, 60, 62, 64, 68, 
70, 77, 79, 80, 82, 83, 92, 93, 96, 100, 106, 107, 111, 113, 120, 127, 130, 133, 135, 136, 140, 141, 144, 152, 154, 155, 160, 161, 163, 167, 168, 
173, 177, 180, 181, 182, 184, 186, 187, 188, 193, 197, 199, 204, 206, 207, 210, 212, 215, 216, 219, 220, 224, 225, 226, 229, 233, 237, 239, 240, 
241, 243, 248, 249, 251, 253, 254, 256, 261, 262, 264, 266, 269, 275, 276, 277, 280, 281, 282, 293, 295, 296, 298, 300, 303, 306, 308, 314, 315, 
316, 317, 319, 322, 325, 330, 331, 334, 337, 341, 344, 347, 348, 350, 353, 354, 357, 358, 360, 361, 362, 364, 370, 383, 384, 385, 389, 390, 392, 
393, 395, 396, 400, 401, 402, 404, 405, 408, 415, 416, 418, 423, 424, 427, 434, 440, 441, 443, 446, 448, 450, 451, 454, 455, 459, 465, 467, 471, 
474, 477, 481, 483, 488, 489, 492, 496, 504, 509, 513, 514, 516, 518, 519, 520, 523, 526, 528, 534, 538, 541, 544, 546, 552, 553, 557, 559]

skip_list_smn_ii = [39, 51, 55, 58, 64, 65, 66, 72, 75, 83, 86, 90, 91, 92]


def import_inmet():

	iy, ix = [], []
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], [], []
	
	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue

		yy = float(inmet[i][2]) 
		xx = float(inmet[i][3])  
		if xx <= -48 and yy <= -16.5:
			iy.append(yy)
			ix.append(xx)
					
			print('Reading weather station:', i, inmet[i][0], inmet[i][1])		
			# reading regcm usp 
			d_0 = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_0 = d_0.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_0 = d_0.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_0 = d_0.groupby('time.year').mean('time')
			d_0 = np.nanmean(d_0.values)
			mean_.append(d_0*86400)
				
			# reading regcm ictp pbl 1 3 km 
			d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
			d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			d_i = np.nanmean(d_i.values)
			mean_i.append(d_i)

			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.year').mean('time')
			d_ii = np.nanmean(d_ii.values)
			mean_ii.append(d_ii*86400)
		
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			d_iii = np.nanmean(d_iii.values)
			mean_iii.append(d_iii*86400)
						
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			d_iv = np.nanmean(d_iv.values)
			mean_iv.append(d_iv*86400)
			
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			d_v = np.nanmean(d_v.values)
			mean_v.append(d_v*86400)
		
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/pre/'.format(path) + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			d_vi = np.nanmean(d_vi.values)
			mean_vi.append(d_vi)
		
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			d_vii = np.nanmean(d_vii.values)		
			mean_vii.append(d_vii)
				
	return iy, ix, mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_i():

	iy, ix = [], []
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], [], []

	for i in range(1, 73):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])		
		# reading regcm usp 
		d_0 = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_0 = d_0.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_0 = d_0.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_0 = d_0.groupby('time.year').mean('time')
		d_0 = np.nanmean(d_0.values)
		mean_.append(d_0*86400)
				
		# reading regcm ictp pbl 1 3 km 
		d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.year').mean('time')
		d_i = np.nanmean(d_i.values)
		mean_i.append(d_i)
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.year').mean('time')
		d_ii = np.nanmean(d_ii.values)
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.year').mean('time')
		d_iii = np.nanmean(d_iii.values)
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.year').mean('time')
		d_iv = np.nanmean(d_iv.values)
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.year').mean('time')
		d_v = np.nanmean(d_v.values)
		mean_v.append(d_v*86400)
		
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		d_vi = np.nanmean(d_vi.values)
		mean_vi.append(d_vi*24)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.year').mean('time')
		d_vii = np.nanmean(d_vii.values)
		mean_vii.append(d_vii)
				
	return iy, ix, mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_ii():
		
	iy, ix = [], []
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], [], []

	for i in range(1, 110):
		if i in skip_list_smn_ii:
			continue
		yy=smn_ii[i][1]
		xx=smn_ii[i][2]
		iy.append(smn_ii[i][1])
		ix.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])		
		# reading regcm usp 
		d_0 = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_0 = d_0.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_0 = d_0.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_0 = d_0.groupby('time.year').mean('time')
		d_0 = np.nanmean(d_0.values)
		mean_.append(d_0*86400)
				
		# reading regcm ictp pbl 1 3 km 
		d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.year').mean('time')
		d_i = np.nanmean(d_i.values)
		mean_i.append(d_i)
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.year').mean('time')
		d_ii = np.nanmean(d_ii.values)
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.year').mean('time')
		d_iii = np.nanmean(d_iii.values)
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.year').mean('time')
		d_iv = np.nanmean(d_iv.values)
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.year').mean('time')
		d_v = np.nanmean(d_v.values)
		mean_v.append(d_v*86400)
		
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		d_vi = np.nanmean(d_vi.values)
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.year').mean('time')
		d_vii = np.nanmean(d_vii.values)
		mean_vii.append(d_vii)
				
	return iy, ix, mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def configure_subplot(ax):

	lon_bounds = [-70, -46]
	lat_bounds = [-36, -16]

	states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

	ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.5)  


	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(font_size)
		
	ax.add_feature(states_provinces, edgecolor='0.05')
	ax.add_feature(cfeat.BORDERS, linewidth=0.75)
	ax.coastlines(linewidth=0.75)
	

# Import dataset
lat_x, lon_x, clim_0_x, clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			
lat_y, lon_y, clim_0_y, clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y = import_smn_i()			
lat_z, lon_z, clim_0_z, clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z, clim_vi_z, clim_vii_z = import_smn_ii()			

lat_yy = lat_x + lat_y + lat_z
lon_xx = lon_x + lon_y + lon_z

reg_usp = clim_0_x + clim_0_y + clim_0_z
reg_ictp = clim_i_x + clim_i_y + clim_i_z
reg_ictp_i = clim_ii_x + clim_ii_y + clim_ii_z
reg_ictp_ii = clim_iii_x + clim_iii_y + clim_iii_z
wrf_ncar = clim_iv_x + clim_iv_y + clim_iv_z
wrf_ucan = clim_v_x + clim_v_y + clim_v_z
inmet_smn = clim_vi_x + clim_vi_y + clim_vi_z
era5 = clim_vii_x + clim_vii_y + clim_vii_z

# Plot figure   
fig, axes = plt.subplots(3,3, figsize=(10, 8), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
fig.delaxes(ax9)

vmin = 0
vmax = 6
n_classes = 16
bins = np.linspace(vmin, vmax, n_classes + 1)
cmap = cm.get_cmap('terrain_r', n_classes)
norm = BoundaryNorm(bins, cmap.N)
	
legend = 'Precipitation (mm d⁻¹)'
states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

st1 = ax1.scatter(lon_xx, lat_yy, 20, inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax1.set_title('(a) INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

st2 = ax2.scatter(lon_xx, lat_yy, 20, era5, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax2.set_title('(b) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

st3 = ax3.scatter(lon_xx, lat_yy, 20, reg_usp, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax3.set_title('(c) Reg4', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

st4 = ax4.scatter(lon_xx, lat_yy, 20, reg_ictp, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax4.set_title('(d) Reg5-Holt3', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

st5 = ax5.scatter(lon_xx, lat_yy, 20, reg_ictp_i, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax5.set_title('(e) Reg5-Holt', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

st6 = ax6.scatter(lon_xx, lat_yy, 20, reg_ictp_ii, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax6.set_title('(f) Reg5-UW', loc='left', fontsize=font_size, fontweight='bold')
ax6.set_xlabel(u'Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

st7 = ax7.scatter(lon_xx, lat_yy, 20, wrf_ncar, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax7.set_title('(f) WRF-NCAR', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
ax7.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

st8 = ax8.scatter(lon_xx, lat_yy, 20, wrf_ucan, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax8.set_title('(g) WRF-UCAN', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cbar = plt.colorbar(st8, cax=fig.add_axes([0.91, 0.25, 0.019, 0.50]), extend='max')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_clim_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

