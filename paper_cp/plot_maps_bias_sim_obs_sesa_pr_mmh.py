# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot maps of bias"

import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.feature as cfeature

from matplotlib.path import Path
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from netCDF4 import Dataset as nc
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

parser = argparse.ArgumentParser(description='Process variable')
parser.add_argument('--var', required=True, choices=['pr'], help='Variable name')
args = parser.parse_args()
var = args.var

dict_var = {'pr': ['pre', 'tp']}

font_size = 7
path = '/home/mda_silv/users/FPS_SESA'

skip_list_inmet_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
skip_list_inmet_ii = [2, 3, 4, 14, 19, 20, 21, 24, 25, 26, 27, 28, 32, 33, 34, 35, 38, 40, 41, 44, 45, 48, 52, 54, 55, 56, 59, 60, 62, 64, 68, 
70, 77, 79, 80, 82, 83, 92, 93, 96, 100, 106, 107, 111, 113, 120, 127, 130, 133, 135, 136, 140, 141, 144, 152, 154, 155, 160, 161, 163, 167, 168, 
173, 177, 180, 181, 182, 184, 186, 187, 188, 193, 197, 199, 204, 206, 207, 210, 212, 215, 216, 219, 220, 224, 225, 226, 229, 233, 237, 239, 240, 
241, 243, 248, 249, 251, 253, 254, 256, 261, 262, 264, 266, 269, 275, 276, 277, 280, 281, 282, 293, 295, 296, 298, 300, 303, 306, 308, 314, 315, 
316, 317, 319, 322, 325, 330, 331, 334, 337, 341, 344, 347, 348, 350, 353, 354, 357, 358, 360, 361, 362, 364, 370, 383, 384, 385, 389, 390, 392, 
393, 395, 396, 400, 401, 402, 404, 405, 408, 415, 416, 418, 423, 424, 427, 434, 440, 441, 443, 446, 448, 450, 451, 454, 455, 459, 465, 467, 471, 
474, 477, 481, 483, 488, 489, 492, 496, 504, 509, 513, 514, 516, 518, 519, 520, 523, 526, 528, 534, 538, 541, 544, 546, 552, 553, 557, 559]

    
def import_inmet():

	iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii, mean_ix = [], [], [], [], [], [], [], [], [], [], []

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
			station_code = inmet[i][0]
			station_name = inmet[i][1]
			print(i, station_code, station_name)
		
			# Reading inmet 
			d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/{1}/'.format(path, dict_var[var][0]) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(dict_var[var][0], station_code))
			d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
			mean_i.append(d_i.values)

			# Reading era5 
			d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
			d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
			mean_ii.append(d_ii.values)

			# Reading regcm usp
			d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_iii.append(d_iii.values/24)
			
			# Reading regcm ictp 
			d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_iv.append(d_iv.values)
	
			# Reading regcm ictp pbl 1
			d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_v.append(d_v.values/24)

			# Reading regcm ictp pbl 2
			d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_vi.append(d_vi.values/24)

			# Reading wrf ncar
			d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_vii.append(d_vii.values/24)
		
			# Reading wrf ucan
			d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_viii.append(d_viii.values/24)
		
			# Reading wrf cima
			d_ix = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_cima/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_ix = d_ix[var].sel(time=slice('2018-06-01','2021-05-31'))
			mean_ix.append(d_ix.values)
					
	return iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii, mean_ix


def import_smn_i():

	iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii, mean_ix = [], [], [], [], [], [], [], [], [], [], []

	for i in range(1, 73):
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		station_code = f'SMN{i:03d}'
		station_name = smn_i[i][0]		
		print(i, station_code, station_name)
		
		# Reading smn 
		d_i = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_i/smn_nc/'.format(path) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(dict_var[var][0], smn_i[i][0]))
		d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
		mean_i.append(d_i.values)

		# Reading era5 
		d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
		d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
		mean_ii.append(d_ii.values)

		# Reading regcm usp
		d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_iii.append(d_iii.values/24)
			
		# Reading regcm ictp 
		d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_iv.append(d_iv.values)
	
		# Reading regcm ictp pbl 1
		d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_v.append(d_v.values/24)

		# Reading regcm ictp pbl 2
		d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_vi.append(d_vi.values/24)

		# Reading wrf ncar
		d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_vii.append(d_vii.values/24)
		
		# Reading wrf ucan
		d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_viii.append(d_viii.values/24)
			
		# Reading wrf cima
		d_ix = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_cima/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_ix = d_ix[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_ix.append(d_ix.values)
					
	return iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii, mean_ix


def compute_stats(stations_data, wet_threshold=0.1, valid_range=(0, 500)):

    mean_list, p99_list, freq_list, int_list = [], [], [], []

    for da in stations_data:
        da = np.asarray(da)
        total_days = len(da)
        
        da = da[np.isfinite(da)]
        da = da[(da >= valid_range[0]) & (da <= valid_range[1])]

        if len(da) == 0:
            mean_list.append(np.nan); p99_list.append(np.nan)
            freq_list.append(np.nan); int_list.append(np.nan)
            continue

        # Mean and p99 (all values)
        mean_list.append(float(np.nanmean(da)))
        p99_list.append(float(np.nanpercentile(da, 99.9)))

        # Wet days (values >= wet_threshold)
        wet_days = da[da >= wet_threshold]
        
        # Frequency (fraction of wet days over total days)
        freq_list.append(float(len(wet_days) / total_days * 100))
        
        # Intensity (mean of wet days)
        int_list.append(float(np.nanmean(wet_days)) if len(wet_days) > 0 else np.nan)

    return mean_list, p99_list, freq_list, int_list
    
    
def configure_subplot(ax):

	lon_bounds = [-62, -46]
	lat_bounds = [-36, -16]

	ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
		
	ax.add_feature(cfeature.OCEAN, facecolor='#a6cee3')
	ax.add_feature(cfeature.LAND, facecolor='lightgray') 
	ax.add_feature(cfeature.BORDERS, linewidth=0.75)
	ax.coastlines(linewidth=0.75)
	

# Import dataset
lat_x, lon_x, mean_i_x, mean_ii_x, mean_iii_x, mean_iv_x, mean_v_x, mean_vi_x, mean_vii_x, mean_viii_x, mean_ix_x = import_inmet()			
lat_y, lon_y, mean_i_y, mean_ii_y, mean_iii_y, mean_iv_y, mean_v_y, mean_vi_y, mean_vii_y, mean_viii_y, mean_ix_y = import_smn_i()			

lat_yy = lat_x + lat_y 
lon_xx = lon_x + lon_y 

inmet_smn   = mean_i_x    + mean_i_y    
era5        = mean_ii_x   + mean_ii_y   
reg_usp     = mean_iii_x  + mean_iii_y  
reg_ictp    = mean_iv_x   + mean_iv_y   
reg_ictp_i  = mean_v_x    + mean_v_y    
reg_ictp_ii = mean_vi_x   + mean_vi_y   
wrf_ncar    = mean_vii_x  + mean_vii_y  
wrf_ucan    = mean_viii_x + mean_viii_y 
wrf_cima    = mean_ix_x   + mean_ix_y   
	
mean_inmet_smn,   perc_inmet_smn,   freq_inmet_smn,   int_inmet_smn   = compute_stats(inmet_smn)
mean_era5,        perc_era5,        freq_era5,        int_era5        = compute_stats(era5)
mean_reg_usp,     perc_reg_usp,     freq_reg_usp,     int_reg_usp     = compute_stats(reg_usp)
mean_reg_ictp,    perc_reg_ictp,    freq_reg_ictp,    int_reg_ictp    = compute_stats(reg_ictp)
mean_reg_ictp_i,  perc_reg_ictp_i,  freq_reg_ictp_i,  int_reg_ictp_i  = compute_stats(reg_ictp_i)
mean_reg_ictp_ii, perc_reg_ictp_ii, freq_reg_ictp_ii, int_reg_ictp_ii = compute_stats(reg_ictp_ii)
mean_wrf_ncar,    perc_wrf_ncar,    freq_wrf_ncar,    int_wrf_ncar    = compute_stats(wrf_ncar)
mean_wrf_ucan,    perc_wrf_ucan,    freq_wrf_ucan,    int_wrf_ucan    = compute_stats(wrf_ucan)                                                                                                                           
mean_wrf_cima,    perc_wrf_cima,    freq_wrf_cima,    int_wrf_cima    = compute_stats(wrf_cima)                                                                                                                           

bias_mean_era5        = np.array(mean_era5)        - np.array(mean_inmet_smn)
bias_mean_reg_usp     = np.array(mean_reg_usp)     - np.array(mean_inmet_smn)
bias_mean_reg_ictp    = np.array(mean_reg_ictp)    - np.array(mean_inmet_smn)
bias_mean_reg_ictp_i  = np.array(mean_reg_ictp_i)  - np.array(mean_inmet_smn)
bias_mean_reg_ictp_ii = np.array(mean_reg_ictp_ii) - np.array(mean_inmet_smn)
bias_mean_wrf_ncar    = np.array(mean_wrf_ncar)    - np.array(mean_inmet_smn)
bias_mean_wrf_ucan    = np.array(mean_wrf_ucan)    - np.array(mean_inmet_smn)
bias_mean_wrf_cima    = np.array(mean_wrf_cima)    - np.array(mean_inmet_smn)

bias_perc_era5        = np.array(perc_era5)        - np.array(perc_inmet_smn)
bias_perc_reg_usp     = np.array(perc_reg_usp)	   - np.array(perc_inmet_smn)
bias_perc_reg_ictp    = np.array(perc_reg_ictp)    - np.array(perc_inmet_smn)
bias_perc_reg_ictp_i  = np.array(perc_reg_ictp_i)  - np.array(perc_inmet_smn)
bias_perc_reg_ictp_ii = np.array(perc_reg_ictp_ii) - np.array(perc_inmet_smn)
bias_perc_wrf_ncar    = np.array(perc_wrf_ncar)    - np.array(perc_inmet_smn)
bias_perc_wrf_ucan    = np.array(perc_wrf_ucan)    - np.array(perc_inmet_smn)
bias_perc_wrf_cima    = np.array(perc_wrf_cima)    - np.array(perc_inmet_smn)

bias_freq_era5        = np.array(freq_era5)        - np.array(freq_inmet_smn)
bias_freq_reg_usp     = np.array(freq_reg_usp)     - np.array(freq_inmet_smn)
bias_freq_reg_ictp    = np.array(freq_reg_ictp)    - np.array(freq_inmet_smn)
bias_freq_reg_ictp_i  = np.array(freq_reg_ictp_i)  - np.array(freq_inmet_smn)
bias_freq_reg_ictp_ii = np.array(freq_reg_ictp_ii) - np.array(freq_inmet_smn)
bias_freq_wrf_ncar    = np.array(freq_wrf_ncar)    - np.array(freq_inmet_smn)
bias_freq_wrf_ucan    = np.array(freq_wrf_ucan)    - np.array(freq_inmet_smn)
bias_freq_wrf_cima    = np.array(freq_wrf_cima)    - np.array(freq_inmet_smn)

bias_int_era5        = np.array(int_era5)        - np.array(int_inmet_smn)
bias_int_reg_usp     = np.array(int_reg_usp)     - np.array(int_inmet_smn)
bias_int_reg_ictp    = np.array(int_reg_ictp)    - np.array(int_inmet_smn)
bias_int_reg_ictp_i  = np.array(int_reg_ictp_i)  - np.array(int_inmet_smn)
bias_int_reg_ictp_ii = np.array(int_reg_ictp_ii) - np.array(int_inmet_smn)
bias_int_wrf_ncar    = np.array(int_wrf_ncar)    - np.array(int_inmet_smn)
bias_int_wrf_ucan    = np.array(int_wrf_ucan)    - np.array(int_inmet_smn)
bias_int_wrf_cima    = np.array(int_wrf_cima)    - np.array(int_inmet_smn)

print(bias_mean_wrf_cima)
print()
print(bias_perc_wrf_cima)
print()
print(bias_freq_wrf_cima)
print()
print(bias_int_wrf_cima)
print()

# Plot figure   
fig, axes = plt.subplots(8,3, figsize=(3.5, 10), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12), (ax13, ax14, ax15), (ax16, ax17, ax18), (ax19, ax20, ax21), (ax22, ax23, ax24)  = axes

cmap = cm.get_cmap('BrBG', 20)
norm_i = BoundaryNorm(np.linspace(-20, 20, 20 + 1), cmap.N)
norm_ii = BoundaryNorm(np.linspace(-10, 10, 20 + 1), cmap.N)
norm_iii = BoundaryNorm(np.linspace(-10, 10, 20 + 1), cmap.N)
legend = 'mm h⁻¹'

ct1 = ax1.scatter(lon_xx, lat_yy, 20, bias_perc_era5, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax1.set_title('(a)', loc='left', fontweight='bold', fontsize=font_size)
ax1.set_ylabel(u'ERA5 - WS', fontsize=font_size)
configure_subplot(ax1)

ct2 = ax2.scatter(lon_xx, lat_yy, 20, bias_freq_era5, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax2.set_title('(b)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax2)

ct3 = ax3.scatter(lon_xx, lat_yy, 20, bias_int_era5, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax3.set_title('(c)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax3)

ct4 = ax4.scatter(lon_xx, lat_yy, 20, bias_perc_reg_usp, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax4.set_title('(d)', loc='left', fontweight='bold', fontsize=font_size)
ax4.set_ylabel('Reg4 - WS', rotation='vertical', fontsize=font_size)
configure_subplot(ax4)

ct5 = ax5.scatter(lon_xx, lat_yy, 20, bias_freq_reg_usp, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax5.set_title('(e)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax5)

ct6 = ax6.scatter(lon_xx, lat_yy, 20, bias_int_reg_usp, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax6.set_title('(f)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax6)

ct7 = ax7.scatter(lon_xx, lat_yy, 20, bias_perc_reg_ictp, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax7.set_title('(g)', loc='left', fontweight='bold', fontsize=font_size)
ax7.set_ylabel(u'Reg5 - WS', fontsize=font_size)
configure_subplot(ax7)

ct8 = ax8.scatter(lon_xx, lat_yy, 20, bias_freq_reg_ictp, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax8.set_title('(h)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax8)

ct9 = ax9.scatter(lon_xx, lat_yy, 20, bias_int_reg_ictp, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax9.set_title('(i)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax9)

ct10 = ax10.scatter(lon_xx, lat_yy, 20, bias_perc_reg_ictp_i, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax10.set_title('(j)', loc='left', fontweight='bold', fontsize=font_size)
ax10.set_ylabel('Reg5-Holt - WS', rotation='vertical', fontsize=font_size)
configure_subplot(ax10)

ct11 = ax11.scatter(lon_xx, lat_yy, 20, bias_freq_reg_ictp_i, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax11.set_title('(k)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax11)

ct12 = ax12.scatter(lon_xx, lat_yy, 20, bias_int_reg_ictp_i, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax12.set_title('(l)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax12)

ct13 = ax13.scatter(lon_xx, lat_yy, 20, bias_perc_reg_ictp_ii, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax13.set_title('(m)', loc='left', fontweight='bold', fontsize=font_size)
ax13.set_ylabel('Reg5-UW - WS', rotation='vertical', fontsize=font_size)
configure_subplot(ax13)

ct14 = ax14.scatter(lon_xx, lat_yy, 20, bias_freq_reg_ictp_ii, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax14.set_title('(n)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax14)

ct15 = ax15.scatter(lon_xx, lat_yy, 20, bias_int_reg_ictp_ii, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax15.set_title('(o)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax15)

ct16 = ax16.scatter(lon_xx, lat_yy, 20, bias_perc_wrf_ncar, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax16.set_title('(p)', loc='left', fontweight='bold', fontsize=font_size)
ax16.set_ylabel('WRF-NCAR - WS', rotation='vertical', fontsize=font_size)
configure_subplot(ax16)

ct17 = ax17.scatter(lon_xx, lat_yy, 20, bias_freq_wrf_ncar, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax17.set_title('(q)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax17)

ct18 = ax18.scatter(lon_xx, lat_yy, 20, bias_int_wrf_ncar, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax18.set_title('(r)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax18)

ct19 = ax19.scatter(lon_xx, lat_yy, 20, bias_perc_wrf_ucan, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax19.set_title('(s)', loc='left', fontweight='bold', fontsize=font_size)
ax19.set_ylabel('WRF-UCAN - WS', rotation='vertical', fontsize=font_size)
configure_subplot(ax19)

ct20 = ax20.scatter(lon_xx, lat_yy, 20, bias_freq_wrf_ucan, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax20.set_title('(t)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax20)

ct21 = ax21.scatter(lon_xx, lat_yy, 20, bias_int_wrf_ucan, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax21.set_title('(u)', loc='left', fontweight='bold', fontsize=font_size)
configure_subplot(ax21)

ct22 = ax22.scatter(lon_xx, lat_yy, 20, bias_perc_wrf_cima, cmap=cmap, norm=norm_i, marker='o', edgecolor='black', linewidth=0.5)
ax22.set_title('(v)', loc='left', fontweight='bold', fontsize=font_size)
ax22.set_xlabel('P99.9 ({0})'.format(legend), loc='center', fontsize=font_size)
ax22.set_ylabel('WRF-CIMA - WS', rotation='vertical', fontsize=font_size)
configure_subplot(ax22)

ct23 = ax23.scatter(lon_xx, lat_yy, 20, bias_freq_wrf_cima, cmap=cmap, norm=norm_ii, marker='o', edgecolor='black', linewidth=0.5)
ax23.set_title('(w)', loc='left', fontweight='bold', fontsize=font_size)
ax23.set_xlabel('FREQUENCY (%)', loc='center', fontsize=font_size)
configure_subplot(ax23)

ct24 = ax24.scatter(lon_xx, lat_yy, 20, bias_int_wrf_cima, cmap=cmap, norm=norm_iii, marker='o', edgecolor='black', linewidth=0.5)
ax24.set_title('(x)', loc='left', fontweight='bold', fontsize=font_size)
ax24.set_xlabel('INTENSITY ({0})'.format(legend), loc='center', fontsize=font_size)
configure_subplot(ax24)

cbar = plt.colorbar(ct1, cax=fig.add_axes([0.92, 0.68, 0.02, 0.20]), extend='neither')
cbar.ax.tick_params(labelsize=font_size)
cbar.set_label('P99 ({0})'.format(legend), fontsize=font_size)

cbar = plt.colorbar(ct2, cax=fig.add_axes([0.92, 0.40, 0.02, 0.20]), extend='neither')
cbar.ax.tick_params(labelsize=font_size)
cbar.set_label('FREQUENCY (%)', fontsize=font_size)

cbar = plt.colorbar(ct3, cax=fig.add_axes([0.92, 0.12, 0.02, 0.20]), extend='neither')
cbar.ax.tick_params(labelsize=font_size)
cbar.set_label('INTENSITY ({0})'.format(legend), fontsize=font_size)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_bias_sesa_{0}_mmh.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()






