# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot maps of bias"

import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

from matplotlib.path import Path
from dict_inmet_stations import inmet
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from netCDF4 import Dataset as nc
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'sfcWind'
dict_var = {'tas': ['tmp', 't2m'], 'sfcWind': ['uv', 'ws10']}

font_size = 8
path = '/home/mda_silv/clima-archive2-b/FPS-SESA'

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
	mean_i, mean_ii, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x, bias_xi, bias_xii = [], [], [], [], [], [], [], [], [], [], [], [], [], []

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
			# Reading inmet 
			d_i = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/inmet/inmet_br/inmet_nc/daily/{0}/'.format(dict_var[var][0]) + '{0}_{1}_D_2018-01-01_2021-12-31.nc'.format(dict_var[var][0], inmet[i][0]))
			d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.groupby('time.year').mean('time')
			values_i = np.nanmean(d_i.values)
			
			# reading era5 
			d_ii = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/era5/' + '{0}_CSAM-4i_ERA5_mon_2018060100-2021053123.nc'.format(dict_var[var][1]))
			d_ii = d_ii[dict_var[var][1]].sel(valid_time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('valid_time.year').mean('valid_time')
			values_ii = np.nanmean(d_ii.values)
					
			# reading regcm usp 
			d_iii = xr.open_dataset('{0}/rcm/reg_usp/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_USP-RegCM471_v2_mon_2018060100-2021053123.nc'.format(var))
			d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = np.nanmean(d_iii.values)
				
			# reading regcm ictp pbl 1 3 km 
			d_iv = xr.open_dataset('{0}/rcm/reg_ictp/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_v1_mon_2018060100-2021053123.nc'.format(var))
			d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = np.nanmean(d_iv.values)

			# reading regcm ictp pbl 1 
			d_v = xr.open_dataset('{0}/rcm/reg_ictp_pbl1/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_2018060100-2021053123.nc'.format(var))
			d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			values_v = np.nanmean(d_v.values)
			
			# reading regcm ictp pbl 2
			d_vi = xr.open_dataset('{0}/rcm/reg_ictp_pbl2/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_2018060100-2021053123.nc'.format(var))
			d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vi = d_vi.groupby('time.year').mean('time')
			values_vi = np.nanmean(d_vi.values)
							
			# reading wrf ncar 
			d_vii = xr.open_dataset('{0}/rcm/wrf_ncar/'.format(path) + '{0}_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_2018060100-2021053123.nc'.format(var))
			d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			values_vii = np.nanmean(d_vii.values)
				
			# reading wrf ucan 
			d_viii = xr.open_dataset('{0}/rcm/wrf_ucan/'.format(path) + '{0}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_2018060100-2021053123.nc'.format(var))
			d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_viii = d_viii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_viii = d_viii.groupby('time.year').mean('time')
			values_viii = np.nanmean(d_viii.values)			

			# calculate bias regcm usp
			bias_i.append(values_iii - values_i)	
			bias_ii.append(values_iii - values_ii)

			# calculate bias regcm ictp
			bias_iii.append(values_iv - values_i)	
			bias_iv.append(values_iv - values_ii)

			# calculate bias regcm ictp 1
			bias_v.append(values_v - values_i)	
			bias_vi.append(values_v - values_ii)

			# calculate bias regcm ictp 2
			bias_vii.append(values_vi - values_i)	
			bias_viii.append(values_vi - values_ii)

			# calculate bias wrf ncar
			bias_ix.append(values_vii - values_i)	
			bias_x.append(values_vii - values_ii)

			# calculate bias wrf ucan
			bias_xi.append(values_viii - values_i)	
			bias_xii.append(values_viii - values_ii)
			
			mean_i.append(values_i)
			mean_ii.append(values_ii)

	return iy, ix, mean_i, mean_ii, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x, bias_xi, bias_xii


def configure_subplot(ax):

	lon_bounds = [-62, -46]
	lat_bounds = [-36, -16]

	states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

	ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.5)  

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(font_size)
		
	ax.add_feature(cfeature.OCEAN, facecolor='#a6cee3')
	ax.add_feature(cfeature.BORDERS, linewidth=0.75)
	ax.add_feature(states_provinces, edgecolor='0.05')
	ax.coastlines(linewidth=0.75)
	

# Specify directories 
dirnc = '/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp'
domname = 'CSAM-3'

# RegCM file
if len(sys.argv) > 1:
    RCMf = nc(sys.argv[1], mode='r')
else:
    RCMf = nc(os.path.join(dirnc,domname+'_DOMAIN000.nc'), mode='r')
    
lat  = RCMf.variables['xlat'][:,:]
lon  = RCMf.variables['xlon'][:,:]
topo = RCMf.variables['topo'][:,:]
lonc = RCMf.longitude_of_projection_origin
latc = RCMf.latitude_of_projection_origin
topo_masked = np.ma.masked_where(topo <= 0, topo)
RCMf.close()
	
# Import dataset
lat_yy, lon_xx, ws_inmet_smn, reanalise_era5, reg_usp_inmet_smn, reg_usp_reanalise, reg_ictp_inmet_smn, reg_ictp_reanalise, reg_ictp_i_inmet_smn, reg_ictp_i_reanalise, reg_ictp_ii_inmet_smn, reg_ictp_ii_reanalise, wrf_ncar_inmet_smn, wrf_ncar_reanalise, wrf_ucan_inmet_smn, wrf_ucan_reanalise = import_inmet()			

print(len(lat_yy))
print(len(lon_xx))
print(len(ws_inmet_smn))
print(ws_inmet_smn)

bm_reg_usp_inmet_smn = np.nanmean(reg_usp_inmet_smn)
bm_reg_usp_reanalise = np.nanmean(reg_usp_reanalise)

bm_reg_ictp_inmet_smn = np.nanmean(reg_ictp_inmet_smn)
bm_reg_ictp_reanalise = np.nanmean(reg_ictp_reanalise)

bm_reg_ictp_i_inmet_smn = np.nanmean(reg_ictp_i_inmet_smn)
bm_reg_ictp_i_reanalise = np.nanmean(reg_ictp_i_reanalise)

bm_reg_ictp_ii_inmet_smn = np.nanmean(reg_ictp_ii_inmet_smn)
bm_reg_ictp_ii_reanalise = np.nanmean(reg_ictp_ii_reanalise)

bm_wrf_ncar_inmet_smn = np.nanmean(wrf_ncar_inmet_smn)
bm_wrf_ncar_reanalise = np.nanmean(wrf_ncar_reanalise)

bm_wrf_ucan_inmet_smn = np.nanmean(wrf_ucan_inmet_smn)
bm_wrf_ucan_reanalise = np.nanmean(wrf_ucan_reanalise)

# Plot figure   
fig, axes = plt.subplots(4,4, figsize=(10, 12), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16) = axes
fig.delaxes(ax3)
fig.delaxes(ax4)

if var == 'tas':
	bins_ = np.linspace(4, 28, 16 + 1)
	cmap_ = cm.get_cmap('jet', 16)
	norm_ = BoundaryNorm(bins_, cmap_.N)
	legend_ = 'Temperature (°C)'
	bins = np.linspace(-3, 3, 16 + 1)
	cmap = cm.get_cmap('bwr', 16)
	norm = BoundaryNorm(bins, cmap.N)
	legend = 'Bias of temperature (°C)'
else:
	bins_ = np.linspace(0, 6, 16 + 1)
	cmap_ = cm.get_cmap('viridis_r', 16)
	norm_ = BoundaryNorm(bins_, cmap_.N)
	legend_ = 'Wind 10m (m s⁻¹)'
	bins = np.linspace(-3, 3, 16 + 1)
	cmap = cm.get_cmap('PRGn', 16)
	norm = BoundaryNorm(bins, cmap.N)
	legend = 'Bias of wind 10m (m s⁻¹)'

ct1 = ax1.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st1 = ax1.scatter(lon_xx, lat_yy, 20, ws_inmet_smn, cmap=cmap_, norm=norm_, marker='o', edgecolor='black', linewidth=0.5)
ax1.set_title('(a) INMET', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

ct2 = ax2.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st2 = ax2.scatter(lon_xx, lat_yy, 20, reanalise_era5, cmap=cmap_, norm=norm_, marker='o', edgecolor='black', linewidth=0.5)
ax2.set_title('(b) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

ct5 = ax5.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st5 = ax5.scatter(lon_xx, lat_yy, 20, reg_usp_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax5.text(-49, -35, '{0}'.format(round(bm_reg_usp_inmet_smn, 1)), color='black', fontsize=font_size, fontweight='bold')
ax5.set_title('(c) Reg4 - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax5.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

ct6 = ax6.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st6 = ax6.scatter(lon_xx, lat_yy, 20, reg_usp_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax6.text(-49, -35, '{0}'.format(round(bm_reg_usp_reanalise, 1)), color='black', fontsize=font_size, fontweight='bold')
ax6.set_title('(d) Reg4 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

ct7 = ax7.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st7 = ax7.scatter(lon_xx, lat_yy, 20, reg_ictp_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax7.text(-49, -35, '{0}'.format(round(bm_reg_ictp_inmet_smn, 1)), color='black', fontsize=font_size, fontweight='bold')
ax7.set_title('(e) Reg5-holt3 - INMET', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

ct8 = ax8.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st8 = ax8.scatter(lon_xx, lat_yy, 20, reg_ictp_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax8.text(-49, -35, '{0}'.format(round(bm_reg_ictp_reanalise, 1)), color='black', fontsize=font_size, fontweight='bold')
ax8.set_title('(f) Reg5-holt3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

ct9 = ax9.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st9 = ax9.scatter(lon_xx, lat_yy, 20, reg_ictp_i_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax9.text(-49, -35, '{0}'.format(round(bm_reg_ictp_i_inmet_smn, 1)), color='black', fontsize=font_size, fontweight='bold')
ax9.set_title('(g) Reg5-holt - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

ct10 = ax10.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st10 = ax10.scatter(lon_xx, lat_yy, 20, reg_ictp_i_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax10.text(-49, -35, '{0}'.format(round(bm_reg_ictp_i_reanalise, 1)), color='black', fontsize=font_size, fontweight='bold')
ax10.set_title('(h) Reg5-holt - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax10)

ct11 = ax11.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st11 = ax11.scatter(lon_xx, lat_yy, 20, reg_ictp_ii_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax11.text(-49, -35, '{0}'.format(round(bm_reg_ictp_ii_inmet_smn, 1)), color='black', fontsize=font_size, fontweight='bold')
ax11.set_title('(i) Reg5-UW - INMET', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax11)

ct12 = ax12.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st12 = ax12.scatter(lon_xx, lat_yy, 20, reg_ictp_ii_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax12.text(-49, -35, '{0}'.format(round(bm_reg_ictp_ii_reanalise, 1)), color='black', fontsize=font_size, fontweight='bold')
ax12.set_title('(j) Reg5-UW - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax12)

ct13 = ax13.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st13 = ax13.scatter(lon_xx, lat_yy, 20, wrf_ncar_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax13.text(-49, -35, '{0}'.format(round(bm_wrf_ncar_inmet_smn, 1)), color='black', fontsize=font_size, fontweight='bold')
ax13.set_title('(k) WRF-NCAR - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax13.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
ax13.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax13)

ct14 = ax14.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st14 = ax14.scatter(lon_xx, lat_yy, 20, wrf_ncar_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax14.text(-49, -35, '{0}'.format(round(bm_wrf_ncar_reanalise, 1)), color='black', fontsize=font_size, fontweight='bold')
ax14.set_title('(l) WRF-NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax14.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax14)

ct15 = ax15.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st15 = ax15.scatter(lon_xx, lat_yy, 20, wrf_ucan_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax15.text(-49, -35, '{0}'.format(round(bm_wrf_ucan_inmet_smn, 1)), color='black', fontsize=font_size, fontweight='bold')
ax15.set_title('(m) WRF-UCAN - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax15.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax15)

ct16 = ax16.contourf(lon, lat, topo_masked, np.arange(0, 2020, 20), cmap='gray_r', extend='max')
st16 = ax16.scatter(lon_xx, lat_yy, 20, wrf_ucan_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax16.text(-49, -35, '{0}'.format(round(bm_wrf_ucan_reanalise, 1)), color='black', fontsize=font_size, fontweight='bold')
ax16.set_title('(n) WRF-UCAN - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax16.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax16)

cbar = plt.colorbar(st1, cax=fig.add_axes([0.91, 0.25, 0.015, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend_), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

cbar = plt.colorbar(st16, cax=fig.add_axes([0.98, 0.25, 0.015, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

cbar_topo = plt.colorbar(ct1, cax=fig.add_axes([0.25, 0.06, 0.53, 0.015]), orientation='horizontal', extend='both')
cbar_topo.set_label('Topography (m)', fontsize=font_size, fontweight='bold')
cbar_topo.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '/home/mda_silv/users/FPS_SESA/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_bias_{0}_sesa_topo.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()






	


