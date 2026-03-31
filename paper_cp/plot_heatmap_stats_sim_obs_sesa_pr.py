# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot boxplot"

import os
import argparse
import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from matplotlib.patches import Polygon
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

parser = argparse.ArgumentParser(description='Process variable')
parser.add_argument('--var', required=True, choices=['tas', 'sfcWind'], help='Variable name')
args = parser.parse_args()
var = args.var

dict_var = {'tas': ['tmp', 't2m'], 'sfcWind': ['uv', 'ws10']}

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
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii = [], [], [], [], [], [], [], []
	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue
		if inmet[i][3] <= -48 and inmet[i][2] <= -16.5:
			station_code = inmet[i][0]
			station_name = inmet[i][1]
			print(station_code, station_name)
		
			# Reading inmet 
			d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/{1}/'.format(path, dict_var[var][0]) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(dict_var[var][0], station_code))
			d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.resample(time='1D').mean()
			d_i = d_i.values
			mean_i.append(d_i)

			# Reading era5 
			d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
			d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.resample(time='1D').mean()
			d_ii = d_ii.values 
			mean_ii.append(d_ii)
			
			# Reading regcm usp
			d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.resample(time='1D').mean()
			d_iii = d_iii.values 
			mean_iii.append(d_iii)

			# Reading regcm ictp 
			d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.resample(time='1D').mean()
			d_iv = d_iv.values 
			mean_iv.append(d_iv)
	
			# Reading regcm ictp pbl 1
			d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.resample(time='1D').mean()
			d_v = d_v.values 
			mean_v.append(d_v)

			# Reading regcm ictp pbl 2
			d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.resample(time='1D').mean()
			d_vi = d_vi.values 
			mean_vi.append(d_vi)

			# Reading wrf ncar
			d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.resample(time='1D').mean()
			d_vii = d_vii.values 
			mean_vii.append(d_vii)
		
			# Reading wrf ucan
			d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_viii = d_viii.resample(time='1D').mean()
			d_viii = d_viii.values 
			mean_viii.append(d_viii)
		
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii
	
					
def compute_stats(stations_data, threshold=None, valid_range=(-50, 60)):

	mean_list = []
	p95_list  = []
	freq_list = []
	int_list  = []

	for da in stations_data:
		# Ensure numpy array
		da = np.asarray(da)

		# Remove NaNs and unreasonable values
		da = da[np.isfinite(da)]
		da = da[(da >= valid_range[0]) & (da <= valid_range[1])]
		
		if len(da) == 0:
			# Handle empty station
			mean_list.append(float('nan'))
			p95_list.append(float('nan'))
			freq_list.append(float('nan'))
			int_list.append(float('nan'))
			continue

		# Annual mean
		mean_val = np.nanmean(da)
		mean_list.append(mean_val)
	
		# 95th percentile
		p95_val = np.nanpercentile(da, 95)
		p95_list.append(p95_val)
		
		#Threshold for freq/intensity
		thr = threshold if threshold is not None else p95_val
		exceed = da >= thr

		# Frequency (% of days above threshold)
		freq_val = np.sum(exceed) / len(da) * 100	
		freq_list.append(freq_val)

		# Intensity (mean of values above threshold)
		int_val = np.nanmean(da[exceed]) if np.any(exceed) else np.nan
		int_list.append(int_val)
	
	return mean_list, p95_list, freq_list, int_list


def compute_rbias(obs_, sim_):

	eps = 1e-10  
	rbias = (sim_ - obs_) / (obs_ + eps) * 100

	return rbias
    

# Import dataset
inmet_smn, era5, reg_usp, reg_ictp, reg_ictp_i, reg_ictp_ii, wrf_ncar, wrf_ucan = import_inmet()				

mean_inmet_smn,   perc_inmet_smn,   freq_inmet_smn,   int_inmet_smn   = compute_stats(inmet_smn)
mean_era5,        perc_era5,        freq_era5,        int_era5        = compute_stats(era5)
mean_reg_usp,     perc_reg_usp,     freq_reg_usp,     int_reg_usp     = compute_stats(reg_usp)
mean_reg_ictp,    perc_reg_ictp,    freq_reg_ictp,    int_reg_ictp    = compute_stats(reg_ictp)
mean_reg_ictp_i,  perc_reg_ictp_i,  freq_reg_ictp_i,  int_reg_ictp_i  = compute_stats(reg_ictp_i)
mean_reg_ictp_ii, perc_reg_ictp_ii, freq_reg_ictp_ii, int_reg_ictp_ii = compute_stats(reg_ictp_ii)
mean_wrf_ncar,    perc_wrf_ncar,    freq_wrf_ncar,    int_wrf_ncar    = compute_stats(wrf_ncar)
mean_wrf_ucan,    perc_wrf_ucan,    freq_wrf_ucan,    int_wrf_ucan    = compute_stats(wrf_ucan)                                                                                                                           

print(len(mean_inmet_smn))
print(len(perc_inmet_smn))
print(len(freq_inmet_smn))
print(len(int_inmet_smn))

bias_mean_era5        = compute_rbias(np.array(mean_era5),        np.array(mean_inmet_smn))
bias_mean_reg_usp     = compute_rbias(np.array(mean_reg_usp),     np.array(mean_inmet_smn))
bias_mean_reg_ictp    = compute_rbias(np.array(mean_reg_ictp),    np.array(mean_inmet_smn))
bias_mean_reg_ictp_i  = compute_rbias(np.array(mean_reg_ictp_i),  np.array(mean_inmet_smn))
bias_mean_reg_ictp_ii = compute_rbias(np.array(mean_reg_ictp_ii), np.array(mean_inmet_smn))
bias_mean_wrf_ncar    = compute_rbias(np.array(mean_wrf_ncar),    np.array(mean_inmet_smn))
bias_mean_wrf_ucan    = compute_rbias(np.array(mean_wrf_ucan),    np.array(mean_inmet_smn))

bias_perc_era5        = compute_rbias(np.array(perc_era5),	  np.array(perc_inmet_smn))
bias_perc_reg_usp     = compute_rbias(np.array(perc_reg_usp),     np.array(perc_inmet_smn))
bias_perc_reg_ictp    = compute_rbias(np.array(perc_reg_ictp),    np.array(perc_inmet_smn))
bias_perc_reg_ictp_i  = compute_rbias(np.array(perc_reg_ictp_i),  np.array(perc_inmet_smn))
bias_perc_reg_ictp_ii = compute_rbias(np.array(perc_reg_ictp_ii), np.array(perc_inmet_smn))
bias_perc_wrf_ncar    = compute_rbias(np.array(perc_wrf_ncar),    np.array(perc_inmet_smn))
bias_perc_wrf_ucan    = compute_rbias(np.array(perc_wrf_ucan),    np.array(perc_inmet_smn))

bias_freq_era5        = compute_rbias(np.array(freq_era5),	  np.array(freq_inmet_smn))
bias_freq_reg_usp     = compute_rbias(np.array(freq_reg_usp),     np.array(freq_inmet_smn))
bias_freq_reg_ictp    = compute_rbias(np.array(freq_reg_ictp),	  np.array(freq_inmet_smn))
bias_freq_reg_ictp_i  = compute_rbias(np.array(freq_reg_ictp_i),  np.array(freq_inmet_smn))
bias_freq_reg_ictp_ii = compute_rbias(np.array(freq_reg_ictp_ii), np.array(freq_inmet_smn))
bias_freq_wrf_ncar    = compute_rbias(np.array(freq_wrf_ncar),    np.array(freq_inmet_smn))
bias_freq_wrf_ucan    = compute_rbias(np.array(freq_wrf_ucan),    np.array(freq_inmet_smn))

bias_int_era5        = compute_rbias(np.array(int_era5),	np.array(int_inmet_smn))
bias_int_reg_usp     = compute_rbias(np.array(int_reg_usp),     np.array(int_inmet_smn))
bias_int_reg_ictp    = compute_rbias(np.array(int_reg_ictp),    np.array(int_inmet_smn))
bias_int_reg_ictp_i  = compute_rbias(np.array(int_reg_ictp_i),  np.array(int_inmet_smn))
bias_int_reg_ictp_ii = compute_rbias(np.array(int_reg_ictp_ii), np.array(int_inmet_smn))
bias_int_wrf_ncar    = compute_rbias(np.array(int_wrf_ncar),    np.array(int_inmet_smn))
bias_int_wrf_ucan    = compute_rbias(np.array(int_wrf_ucan),    np.array(int_inmet_smn))

print()
print(len(bias_perc_era5))
print(bias_freq_era5)
print()

list_hc = [1, 2, 3, 2, 0, 1, 1, 0, 2, 2, 0, 3, 0, 2, 3, 0, 1, 2, 0, 3, 0, 4, 2, 4, 3, 1, 4, 2, 4, 2, 2, 2, 1, 2, 4, 2, 2, 3, 2, 4, 4, 4, 0, 2, 4, 3, 2, 0, 0, 0, 3, 2, 2, 2, 1, 2, 4, 1, 4, 3, 4, 3, 0, 2, 0, 3, 2, 3, 2, 4, 0, 1, 4, 2, 4, 4, 0, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 3, 2, 0, 0, 0, 4, 2, 3, 2, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 4, 0, 0, 4, 0, 4, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 2, 4, 3, 1, 4, 1, 2, 1, 1, 1, 4, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 2, 2, 2, 2, 2]
list_hc = list_hc[:len(inmet_smn)]

print(len(inmet_smn))
print(len(list_hc))

count_i, count_ii, count_iii, count_iv, count_v = [], [], [], [], []
for count, idx in enumerate(list_hc):
	if idx == 0:
		count_i.append(count)
	if idx == 1:
		count_ii.append(count)
	if idx == 2:
		count_iii.append(count)
	if idx == 3:
		count_iv.append(count)
	if idx == 4:
		count_v.append(count)

mean_era5_ci, mean_era5_cii, mean_era5_ciii, mean_era5_civ, mean_era5_cv = [], [], [], [], []
perc_era5_ci, perc_era5_cii, perc_era5_ciii, perc_era5_civ, perc_era5_cv = [], [], [], [], []
freq_era5_ci, freq_era5_cii, freq_era5_ciii, freq_era5_civ, freq_era5_cv = [], [], [], [], []
int_era5_ci,  int_era5_cii,  int_era5_ciii,  int_era5_civ,  int_era5_cv  = [], [], [], [], []

mean_reg_usp_ci, mean_reg_usp_cii, mean_reg_usp_ciii, mean_reg_usp_civ, mean_reg_usp_cv = [], [], [], [], []
perc_reg_usp_ci, perc_reg_usp_cii, perc_reg_usp_ciii, perc_reg_usp_civ, perc_reg_usp_cv = [], [], [], [], []
freq_reg_usp_ci, freq_reg_usp_cii, freq_reg_usp_ciii, freq_reg_usp_civ, freq_reg_usp_cv = [], [], [], [], []
int_reg_usp_ci,  int_reg_usp_cii,  int_reg_usp_ciii,  int_reg_usp_civ,  int_reg_usp_cv  = [], [], [], [], []

mean_reg_ictp_ci, mean_reg_ictp_cii, mean_reg_ictp_ciii, mean_reg_ictp_civ, mean_reg_ictp_cv = [], [], [], [], []
perc_reg_ictp_ci, perc_reg_ictp_cii, perc_reg_ictp_ciii, perc_reg_ictp_civ, perc_reg_ictp_cv = [], [], [], [], []
freq_reg_ictp_ci, freq_reg_ictp_cii, freq_reg_ictp_ciii, freq_reg_ictp_civ, freq_reg_ictp_cv = [], [], [], [], []
int_reg_ictp_ci,  int_reg_ictp_cii,  int_reg_ictp_ciii,  int_reg_ictp_civ,  int_reg_ictp_cv  = [], [], [], [], []

mean_reg_ictp_i_ci, mean_reg_ictp_i_cii, mean_reg_ictp_i_ciii, mean_reg_ictp_i_civ, mean_reg_ictp_i_cv = [], [], [], [], []
perc_reg_ictp_i_ci, perc_reg_ictp_i_cii, perc_reg_ictp_i_ciii, perc_reg_ictp_i_civ, perc_reg_ictp_i_cv = [], [], [], [], []
freq_reg_ictp_i_ci, freq_reg_ictp_i_cii, freq_reg_ictp_i_ciii, freq_reg_ictp_i_civ, freq_reg_ictp_i_cv = [], [], [], [], []
int_reg_ictp_i_ci,  int_reg_ictp_i_cii,  int_reg_ictp_i_ciii,  int_reg_ictp_i_civ,  int_reg_ictp_i_cv  = [], [], [], [], []

mean_reg_ictp_ii_ci, mean_reg_ictp_ii_cii, mean_reg_ictp_ii_ciii, mean_reg_ictp_ii_civ, mean_reg_ictp_ii_cv = [], [], [], [], []
perc_reg_ictp_ii_ci, perc_reg_ictp_ii_cii, perc_reg_ictp_ii_ciii, perc_reg_ictp_ii_civ, perc_reg_ictp_ii_cv = [], [], [], [], []
freq_reg_ictp_ii_ci, freq_reg_ictp_ii_cii, freq_reg_ictp_ii_ciii, freq_reg_ictp_ii_civ, freq_reg_ictp_ii_cv = [], [], [], [], []
int_reg_ictp_ii_ci,  int_reg_ictp_ii_cii,  int_reg_ictp_ii_ciii,  int_reg_ictp_ii_civ,  int_reg_ictp_ii_cv  = [], [], [], [], []

mean_wrf_ncar_ci, mean_wrf_ncar_cii, mean_wrf_ncar_ciii, mean_wrf_ncar_civ, mean_wrf_ncar_cv = [], [], [], [], []
perc_wrf_ncar_ci, perc_wrf_ncar_cii, perc_wrf_ncar_ciii, perc_wrf_ncar_civ, perc_wrf_ncar_cv = [], [], [], [], []
freq_wrf_ncar_ci, freq_wrf_ncar_cii, freq_wrf_ncar_ciii, freq_wrf_ncar_civ, freq_wrf_ncar_cv = [], [], [], [], []
int_wrf_ncar_ci,  int_wrf_ncar_cii,  int_wrf_ncar_ciii,  int_wrf_ncar_civ,  int_wrf_ncar_cv  = [], [], [], [], []

mean_wrf_ucan_ci, mean_wrf_ucan_cii, mean_wrf_ucan_ciii, mean_wrf_ucan_civ, mean_wrf_ucan_cv = [], [], [], [], []
perc_wrf_ucan_ci, perc_wrf_ucan_cii, perc_wrf_ucan_ciii, perc_wrf_ucan_civ, perc_wrf_ucan_cv = [], [], [], [], []
freq_wrf_ucan_ci, freq_wrf_ucan_cii, freq_wrf_ucan_ciii, freq_wrf_ucan_civ, freq_wrf_ucan_cv = [], [], [], [], []
int_wrf_ucan_ci,  int_wrf_ucan_cii,  int_wrf_ucan_ciii,  int_wrf_ucan_civ,  int_wrf_ucan_cv  = [], [], [], [], []

for c_i in count_i:
	mean_era5_ci.append(bias_mean_era5[c_i])
	perc_era5_ci.append(bias_perc_era5[c_i])
	freq_era5_ci.append(bias_freq_era5[c_i])
	int_era5_ci.append(bias_int_era5[c_i])

	mean_reg_usp_ci.append(bias_mean_reg_usp[c_i])
	perc_reg_usp_ci.append(bias_perc_reg_usp[c_i])
	freq_reg_usp_ci.append(bias_freq_reg_usp[c_i])
	int_reg_usp_ci.append(bias_int_reg_usp[c_i])

	mean_reg_ictp_ci.append(bias_mean_reg_ictp[c_i])
	perc_reg_ictp_ci.append(bias_perc_reg_ictp[c_i])
	freq_reg_ictp_ci.append(bias_freq_reg_ictp[c_i])
	int_reg_ictp_ci.append(bias_int_reg_ictp[c_i])
	
	mean_reg_ictp_i_ci.append(bias_mean_reg_ictp_i[c_i])
	perc_reg_ictp_i_ci.append(bias_perc_reg_ictp_i[c_i])
	freq_reg_ictp_i_ci.append(bias_freq_reg_ictp_i[c_i])
	int_reg_ictp_i_ci.append(bias_int_reg_ictp_i[c_i])
	
	mean_reg_ictp_ii_ci.append(bias_mean_reg_ictp_ii[c_i])
	perc_reg_ictp_ii_ci.append(bias_perc_reg_ictp_ii[c_i])
	freq_reg_ictp_ii_ci.append(bias_freq_reg_ictp_ii[c_i])
	int_reg_ictp_ii_ci.append(bias_int_reg_ictp_ii[c_i])

	mean_wrf_ncar_ci.append(bias_mean_wrf_ncar[c_i])
	perc_wrf_ncar_ci.append(bias_perc_wrf_ncar[c_i])
	freq_wrf_ncar_ci.append(bias_freq_wrf_ncar[c_i])
	int_wrf_ncar_ci.append(bias_int_wrf_ncar[c_i])
	
	mean_wrf_ucan_ci.append(bias_mean_wrf_ucan[c_i])
	perc_wrf_ucan_ci.append(bias_perc_wrf_ucan[c_i])
	freq_wrf_ucan_ci.append(bias_freq_wrf_ucan[c_i])
	int_wrf_ucan_ci.append(bias_int_wrf_ucan[c_i])

for c_ii in count_ii:
	mean_era5_cii.append(bias_mean_era5[c_ii])
	perc_era5_cii.append(bias_perc_era5[c_ii])
	freq_era5_cii.append(bias_freq_era5[c_ii])
	int_era5_cii.append(bias_int_era5[c_ii])

	mean_reg_usp_cii.append(bias_mean_reg_usp[c_ii])
	perc_reg_usp_cii.append(bias_perc_reg_usp[c_ii])
	freq_reg_usp_cii.append(bias_freq_reg_usp[c_ii])
	int_reg_usp_cii.append(bias_int_reg_usp[c_ii])

	mean_reg_ictp_cii.append(bias_mean_reg_ictp[c_ii])
	perc_reg_ictp_cii.append(bias_perc_reg_ictp[c_ii])
	freq_reg_ictp_cii.append(bias_freq_reg_ictp[c_ii])
	int_reg_ictp_cii.append(bias_int_reg_ictp[c_ii])

	mean_reg_ictp_i_cii.append(bias_mean_reg_ictp_i[c_ii])
	perc_reg_ictp_i_cii.append(bias_perc_reg_ictp_i[c_ii])
	freq_reg_ictp_i_cii.append(bias_freq_reg_ictp_i[c_ii])
	int_reg_ictp_i_cii.append(bias_int_reg_ictp_i[c_ii])

	mean_reg_ictp_ii_cii.append(bias_mean_reg_ictp_ii[c_ii])
	perc_reg_ictp_ii_cii.append(bias_perc_reg_ictp_ii[c_ii])
	freq_reg_ictp_ii_cii.append(bias_freq_reg_ictp_ii[c_ii])
	int_reg_ictp_ii_cii.append(bias_int_reg_ictp_ii[c_ii])

	mean_wrf_ncar_cii.append(bias_mean_wrf_ncar[c_ii])
	perc_wrf_ncar_cii.append(bias_perc_wrf_ncar[c_ii])
	freq_wrf_ncar_cii.append(bias_freq_wrf_ncar[c_ii])
	int_wrf_ncar_cii.append(bias_int_wrf_ncar[c_ii])

	mean_wrf_ucan_cii.append(bias_mean_wrf_ucan[c_ii])
	perc_wrf_ucan_cii.append(bias_perc_wrf_ucan[c_ii])
	freq_wrf_ucan_cii.append(bias_freq_wrf_ucan[c_ii])
	int_wrf_ucan_cii.append(bias_int_wrf_ucan[c_ii])
	
for c_iii in count_iii:
	mean_era5_ciii.append(bias_mean_era5[c_iii])
	perc_era5_ciii.append(bias_perc_era5[c_iii])
	freq_era5_ciii.append(bias_freq_era5[c_iii])
	int_era5_ciii.append(bias_int_era5[c_iii])

	mean_reg_usp_ciii.append(bias_mean_reg_usp[c_iii])
	perc_reg_usp_ciii.append(bias_perc_reg_usp[c_iii])
	freq_reg_usp_ciii.append(bias_freq_reg_usp[c_iii])
	int_reg_usp_ciii.append(bias_int_reg_usp[c_iii])

	mean_reg_ictp_ciii.append(bias_mean_reg_ictp[c_iii])
	perc_reg_ictp_ciii.append(bias_perc_reg_ictp[c_iii])
	freq_reg_ictp_ciii.append(bias_freq_reg_ictp[c_iii])
	int_reg_ictp_ciii.append(bias_int_reg_ictp[c_iii])

	mean_reg_ictp_i_ciii.append(bias_mean_reg_ictp_i[c_iii])
	perc_reg_ictp_i_ciii.append(bias_perc_reg_ictp_i[c_iii])
	freq_reg_ictp_i_ciii.append(bias_freq_reg_ictp_i[c_iii])
	int_reg_ictp_i_ciii.append(bias_int_reg_ictp_i[c_iii])

	mean_reg_ictp_ii_ciii.append(bias_mean_reg_ictp_ii[c_iii])
	perc_reg_ictp_ii_ciii.append(bias_perc_reg_ictp_ii[c_iii])
	freq_reg_ictp_ii_ciii.append(bias_freq_reg_ictp_ii[c_iii])
	int_reg_ictp_ii_ciii.append(bias_int_reg_ictp_ii[c_iii])

	mean_wrf_ncar_ciii.append(bias_mean_wrf_ncar[c_iii])
	perc_wrf_ncar_ciii.append(bias_perc_wrf_ncar[c_iii])
	freq_wrf_ncar_ciii.append(bias_freq_wrf_ncar[c_iii])
	int_wrf_ncar_ciii.append(bias_int_wrf_ncar[c_iii])

	mean_wrf_ucan_ciii.append(bias_mean_wrf_ucan[c_iii])
	perc_wrf_ucan_ciii.append(bias_perc_wrf_ucan[c_iii])
	freq_wrf_ucan_ciii.append(bias_freq_wrf_ucan[c_iii])
	int_wrf_ucan_ciii.append(bias_int_wrf_ucan[c_iii])
	
for c_iv in count_iv:
	mean_era5_civ.append(bias_mean_era5[c_iv])
	perc_era5_civ.append(bias_perc_era5[c_iv])
	freq_era5_civ.append(bias_freq_era5[c_iv])
	int_era5_civ.append(bias_int_era5[c_iv])

	mean_reg_usp_civ.append(bias_mean_reg_usp[c_iv])
	perc_reg_usp_civ.append(bias_perc_reg_usp[c_iv])
	freq_reg_usp_civ.append(bias_freq_reg_usp[c_iv])
	int_reg_usp_civ.append(bias_int_reg_usp[c_iv])

	mean_reg_ictp_civ.append(bias_mean_reg_ictp[c_iv])
	perc_reg_ictp_civ.append(bias_perc_reg_ictp[c_iv])
	freq_reg_ictp_civ.append(bias_freq_reg_ictp[c_iv])
	int_reg_ictp_civ.append(bias_int_reg_ictp[c_iv])

	mean_reg_ictp_i_civ.append(bias_mean_reg_ictp_i[c_iv])
	perc_reg_ictp_i_civ.append(bias_perc_reg_ictp_i[c_iv])
	freq_reg_ictp_i_civ.append(bias_freq_reg_ictp_i[c_iv])
	int_reg_ictp_i_civ.append(bias_int_reg_ictp_i[c_iv])

	mean_reg_ictp_ii_civ.append(bias_mean_reg_ictp_ii[c_iv])
	perc_reg_ictp_ii_civ.append(bias_perc_reg_ictp_ii[c_iv])
	freq_reg_ictp_ii_civ.append(bias_freq_reg_ictp_ii[c_iv])
	int_reg_ictp_ii_civ.append(bias_int_reg_ictp_ii[c_iv])

	mean_wrf_ncar_civ.append(bias_mean_wrf_ncar[c_iv])
	perc_wrf_ncar_civ.append(bias_perc_wrf_ncar[c_iv])
	freq_wrf_ncar_civ.append(bias_freq_wrf_ncar[c_iv])
	int_wrf_ncar_civ.append(bias_int_wrf_ncar[c_iv])

	mean_wrf_ucan_civ.append(bias_mean_wrf_ucan[c_iv])
	perc_wrf_ucan_civ.append(bias_perc_wrf_ucan[c_iv])
	freq_wrf_ucan_civ.append(bias_freq_wrf_ucan[c_iv])
	int_wrf_ucan_civ.append(bias_int_wrf_ucan[c_iv])
	
for c_v in count_v:
	mean_era5_cv.append(bias_mean_era5[c_v])
	perc_era5_cv.append(bias_perc_era5[c_v])
	freq_era5_cv.append(bias_freq_era5[c_v])
	int_era5_cv.append(bias_int_era5[c_v])

	mean_reg_usp_cv.append(bias_mean_reg_usp[c_v])
	perc_reg_usp_cv.append(bias_perc_reg_usp[c_v])
	freq_reg_usp_cv.append(bias_freq_reg_usp[c_v])
	int_reg_usp_cv.append(bias_int_reg_usp[c_v])

	mean_reg_ictp_cv.append(bias_mean_reg_ictp[c_v])
	perc_reg_ictp_cv.append(bias_perc_reg_ictp[c_v])
	freq_reg_ictp_cv.append(bias_freq_reg_ictp[c_v])
	int_reg_ictp_cv.append(bias_int_reg_ictp[c_v])

	mean_reg_ictp_i_cv.append(bias_mean_reg_ictp_i[c_v])
	perc_reg_ictp_i_cv.append(bias_perc_reg_ictp_i[c_v])
	freq_reg_ictp_i_cv.append(bias_freq_reg_ictp_i[c_v])
	int_reg_ictp_i_cv.append(bias_int_reg_ictp_i[c_v])

	mean_reg_ictp_ii_cv.append(bias_mean_reg_ictp_ii[c_v])
	perc_reg_ictp_ii_cv.append(bias_perc_reg_ictp_ii[c_v])
	freq_reg_ictp_ii_cv.append(bias_freq_reg_ictp_ii[c_v])
	int_reg_ictp_ii_cv.append(bias_int_reg_ictp_ii[c_v])

	mean_wrf_ncar_cv.append(bias_mean_wrf_ncar[c_v])
	perc_wrf_ncar_cv.append(bias_perc_wrf_ncar[c_v])
	freq_wrf_ncar_cv.append(bias_freq_wrf_ncar[c_v])
	int_wrf_ncar_cv.append(bias_int_wrf_ncar[c_v])

	mean_wrf_ucan_cv.append(bias_mean_wrf_ucan[c_v])
	perc_wrf_ucan_cv.append(bias_perc_wrf_ucan[c_v])
	freq_wrf_ucan_cv.append(bias_freq_wrf_ucan[c_v])
	int_wrf_ucan_cv.append(bias_int_wrf_ucan[c_v])

# Group I
mean_era5_c_i = np.nanmean(mean_era5_ci, axis=0)
perc_era5_c_i = np.nanmean(perc_era5_ci, axis=0)
freq_era5_c_i = np.nanmean(freq_era5_ci, axis=0)
int_era5_c_i  = np.nanmean(int_era5_ci, axis=0)

mean_reg_usp_c_i = np.nanmean(mean_reg_usp_ci, axis=0)
perc_reg_usp_c_i = np.nanmean(perc_reg_usp_ci, axis=0)
freq_reg_usp_c_i = np.nanmean(freq_reg_usp_ci, axis=0)
int_reg_usp_c_i  = np.nanmean(int_reg_usp_ci, axis=0)

mean_reg_ictp_c_i = np.nanmean(mean_reg_ictp_ci, axis=0)
perc_reg_ictp_c_i = np.nanmean(perc_reg_ictp_ci, axis=0)
freq_reg_ictp_c_i = np.nanmean(freq_reg_ictp_ci, axis=0)
int_reg_ictp_c_i  = np.nanmean(int_reg_ictp_ci, axis=0)

mean_reg_ictp_i_c_i = np.nanmean(mean_reg_ictp_i_ci, axis=0)
perc_reg_ictp_i_c_i = np.nanmean(perc_reg_ictp_i_ci, axis=0)
freq_reg_ictp_i_c_i = np.nanmean(freq_reg_ictp_i_ci, axis=0)
int_reg_ictp_i_c_i  = np.nanmean(int_reg_ictp_i_ci, axis=0)

mean_reg_ictp_ii_c_i = np.nanmean(mean_reg_ictp_ii_ci, axis=0)
perc_reg_ictp_ii_c_i = np.nanmean(perc_reg_ictp_ii_ci, axis=0)
freq_reg_ictp_ii_c_i = np.nanmean(freq_reg_ictp_ii_ci, axis=0)
int_reg_ictp_ii_c_i  = np.nanmean(int_reg_ictp_ii_ci, axis=0)

mean_wrf_ncar_c_i = np.nanmean(mean_wrf_ncar_ci, axis=0)
perc_wrf_ncar_c_i = np.nanmean(perc_wrf_ncar_ci, axis=0)
freq_wrf_ncar_c_i = np.nanmean(freq_wrf_ncar_ci, axis=0)
int_wrf_ncar_c_i  = np.nanmean(int_wrf_ncar_ci, axis=0)

mean_wrf_ucan_c_i = np.nanmean(mean_wrf_ucan_ci, axis=0)
perc_wrf_ucan_c_i = np.nanmean(perc_wrf_ucan_ci, axis=0)
freq_wrf_ucan_c_i = np.nanmean(freq_wrf_ucan_ci, axis=0)
int_wrf_ucan_c_i  = np.nanmean(int_wrf_ucan_ci, axis=0)


print(mean_era5_ci)
print(len(mean_era5_ci))
print(mean_era5_c_i)

# Group II
mean_era5_c_ii = np.nanmean(mean_era5_cii, axis=0)
perc_era5_c_ii = np.nanmean(perc_era5_cii, axis=0)
freq_era5_c_ii = np.nanmean(freq_era5_cii, axis=0)
int_era5_c_ii  = np.nanmean(int_era5_cii, axis=0)

mean_reg_usp_c_ii = np.nanmean(mean_reg_usp_cii, axis=0)
perc_reg_usp_c_ii = np.nanmean(perc_reg_usp_cii, axis=0)
freq_reg_usp_c_ii = np.nanmean(freq_reg_usp_cii, axis=0)
int_reg_usp_c_ii  = np.nanmean(int_reg_usp_cii, axis=0)

mean_reg_ictp_c_ii = np.nanmean(mean_reg_ictp_cii, axis=0)
perc_reg_ictp_c_ii = np.nanmean(perc_reg_ictp_cii, axis=0)
freq_reg_ictp_c_ii = np.nanmean(freq_reg_ictp_cii, axis=0)
int_reg_ictp_c_ii  = np.nanmean(int_reg_ictp_cii, axis=0)

mean_reg_ictp_i_c_ii = np.nanmean(mean_reg_ictp_i_cii, axis=0)
perc_reg_ictp_i_c_ii = np.nanmean(perc_reg_ictp_i_cii, axis=0)
freq_reg_ictp_i_c_ii = np.nanmean(freq_reg_ictp_i_cii, axis=0)
int_reg_ictp_i_c_ii  = np.nanmean(int_reg_ictp_i_cii, axis=0)

mean_reg_ictp_ii_c_ii = np.nanmean(mean_reg_ictp_ii_cii, axis=0)
perc_reg_ictp_ii_c_ii = np.nanmean(perc_reg_ictp_ii_cii, axis=0)
freq_reg_ictp_ii_c_ii = np.nanmean(freq_reg_ictp_ii_cii, axis=0)
int_reg_ictp_ii_c_ii  = np.nanmean(int_reg_ictp_ii_cii, axis=0)

mean_wrf_ncar_c_ii = np.nanmean(mean_wrf_ncar_cii, axis=0)
perc_wrf_ncar_c_ii = np.nanmean(perc_wrf_ncar_cii, axis=0)
freq_wrf_ncar_c_ii = np.nanmean(freq_wrf_ncar_cii, axis=0)
int_wrf_ncar_c_ii  = np.nanmean(int_wrf_ncar_cii, axis=0)

mean_wrf_ucan_c_ii = np.nanmean(mean_wrf_ucan_cii, axis=0)
perc_wrf_ucan_c_ii = np.nanmean(perc_wrf_ucan_cii, axis=0)
freq_wrf_ucan_c_ii = np.nanmean(freq_wrf_ucan_cii, axis=0)
int_wrf_ucan_c_ii  = np.nanmean(int_wrf_ucan_cii, axis=0)

# Group III
mean_era5_c_iii = np.nanmean(mean_era5_ciii, axis=0)
perc_era5_c_iii = np.nanmean(perc_era5_ciii, axis=0)
freq_era5_c_iii = np.nanmean(freq_era5_ciii, axis=0)
int_era5_c_iii  = np.nanmean(int_era5_ciii, axis=0)

mean_reg_usp_c_iii = np.nanmean(mean_reg_usp_ciii, axis=0)
perc_reg_usp_c_iii = np.nanmean(perc_reg_usp_ciii, axis=0)
freq_reg_usp_c_iii = np.nanmean(freq_reg_usp_ciii, axis=0)
int_reg_usp_c_iii  = np.nanmean(int_reg_usp_ciii, axis=0)

mean_reg_ictp_c_iii = np.nanmean(mean_reg_ictp_ciii, axis=0)
perc_reg_ictp_c_iii = np.nanmean(perc_reg_ictp_ciii, axis=0)
freq_reg_ictp_c_iii = np.nanmean(freq_reg_ictp_ciii, axis=0)
int_reg_ictp_c_iii  = np.nanmean(int_reg_ictp_ciii, axis=0)

mean_reg_ictp_i_c_iii = np.nanmean(mean_reg_ictp_i_ciii, axis=0)
perc_reg_ictp_i_c_iii = np.nanmean(perc_reg_ictp_i_ciii, axis=0)
freq_reg_ictp_i_c_iii = np.nanmean(freq_reg_ictp_i_ciii, axis=0)
int_reg_ictp_i_c_iii  = np.nanmean(int_reg_ictp_i_ciii, axis=0)

mean_reg_ictp_ii_c_iii = np.nanmean(mean_reg_ictp_ii_ciii, axis=0)
perc_reg_ictp_ii_c_iii = np.nanmean(perc_reg_ictp_ii_ciii, axis=0)
freq_reg_ictp_ii_c_iii = np.nanmean(freq_reg_ictp_ii_ciii, axis=0)
int_reg_ictp_ii_c_iii  = np.nanmean(int_reg_ictp_ii_ciii, axis=0)

mean_wrf_ncar_c_iii = np.nanmean(mean_wrf_ncar_ciii, axis=0)
perc_wrf_ncar_c_iii = np.nanmean(perc_wrf_ncar_ciii, axis=0)
freq_wrf_ncar_c_iii = np.nanmean(freq_wrf_ncar_ciii, axis=0)
int_wrf_ncar_c_iii  = np.nanmean(int_wrf_ncar_ciii, axis=0)

mean_wrf_ucan_c_iii = np.nanmean(mean_wrf_ucan_ciii, axis=0)
perc_wrf_ucan_c_iii = np.nanmean(perc_wrf_ucan_ciii, axis=0)
freq_wrf_ucan_c_iii = np.nanmean(freq_wrf_ucan_ciii, axis=0)
int_wrf_ucan_c_iii  = np.nanmean(int_wrf_ucan_ciii, axis=0)

# Group IV
mean_era5_c_iv = np.nanmean(mean_era5_civ, axis=0)
perc_era5_c_iv = np.nanmean(perc_era5_civ, axis=0)
freq_era5_c_iv = np.nanmean(freq_era5_civ, axis=0)
int_era5_c_iv  = np.nanmean(int_era5_civ, axis=0)

mean_reg_usp_c_iv = np.nanmean(mean_reg_usp_civ, axis=0)
perc_reg_usp_c_iv = np.nanmean(perc_reg_usp_civ, axis=0)
freq_reg_usp_c_iv = np.nanmean(freq_reg_usp_civ, axis=0)
int_reg_usp_c_iv  = np.nanmean(int_reg_usp_civ, axis=0)

mean_reg_ictp_c_iv = np.nanmean(mean_reg_ictp_civ, axis=0)
perc_reg_ictp_c_iv = np.nanmean(perc_reg_ictp_civ, axis=0)
freq_reg_ictp_c_iv = np.nanmean(freq_reg_ictp_civ, axis=0)
int_reg_ictp_c_iv  = np.nanmean(int_reg_ictp_civ, axis=0)

mean_reg_ictp_i_c_iv = np.nanmean(mean_reg_ictp_i_civ, axis=0)
perc_reg_ictp_i_c_iv = np.nanmean(perc_reg_ictp_i_civ, axis=0)
freq_reg_ictp_i_c_iv = np.nanmean(freq_reg_ictp_i_civ, axis=0)
int_reg_ictp_i_c_iv  = np.nanmean(int_reg_ictp_i_civ, axis=0)

mean_reg_ictp_ii_c_iv = np.nanmean(mean_reg_ictp_ii_civ, axis=0)
perc_reg_ictp_ii_c_iv = np.nanmean(perc_reg_ictp_ii_civ, axis=0)
freq_reg_ictp_ii_c_iv = np.nanmean(freq_reg_ictp_ii_civ, axis=0)
int_reg_ictp_ii_c_iv  = np.nanmean(int_reg_ictp_ii_civ, axis=0)

mean_wrf_ncar_c_iv = np.nanmean(mean_wrf_ncar_civ, axis=0)
perc_wrf_ncar_c_iv = np.nanmean(perc_wrf_ncar_civ, axis=0)
freq_wrf_ncar_c_iv = np.nanmean(freq_wrf_ncar_civ, axis=0)
int_wrf_ncar_c_iv  = np.nanmean(int_wrf_ncar_civ, axis=0)

mean_wrf_ucan_c_iv = np.nanmean(mean_wrf_ucan_civ, axis=0)
perc_wrf_ucan_c_iv = np.nanmean(perc_wrf_ucan_civ, axis=0)
freq_wrf_ucan_c_iv = np.nanmean(freq_wrf_ucan_civ, axis=0)
int_wrf_ucan_c_iv  = np.nanmean(int_wrf_ucan_civ, axis=0)

# Group V
mean_era5_c_v = np.nanmean(mean_era5_cv, axis=0)
perc_era5_c_v = np.nanmean(perc_era5_cv, axis=0)
freq_era5_c_v = np.nanmean(freq_era5_cv, axis=0)
int_era5_c_v  = np.nanmean(int_era5_cv, axis=0)

mean_reg_usp_c_v = np.nanmean(mean_reg_usp_cv, axis=0)
perc_reg_usp_c_v = np.nanmean(perc_reg_usp_cv, axis=0)
freq_reg_usp_c_v = np.nanmean(freq_reg_usp_cv, axis=0)
int_reg_usp_c_v  = np.nanmean(int_reg_usp_cv, axis=0)

mean_reg_ictp_c_v = np.nanmean(mean_reg_ictp_cv, axis=0)
perc_reg_ictp_c_v = np.nanmean(perc_reg_ictp_cv, axis=0)
freq_reg_ictp_c_v = np.nanmean(freq_reg_ictp_cv, axis=0)
int_reg_ictp_c_v  = np.nanmean(int_reg_ictp_cv, axis=0)

mean_reg_ictp_i_c_v = np.nanmean(mean_reg_ictp_i_cv, axis=0)
perc_reg_ictp_i_c_v = np.nanmean(perc_reg_ictp_i_cv, axis=0)
freq_reg_ictp_i_c_v = np.nanmean(freq_reg_ictp_i_cv, axis=0)
int_reg_ictp_i_c_v  = np.nanmean(int_reg_ictp_i_cv, axis=0)

mean_reg_ictp_ii_c_v = np.nanmean(mean_reg_ictp_ii_cv, axis=0)
perc_reg_ictp_ii_c_v = np.nanmean(perc_reg_ictp_ii_cv, axis=0)
freq_reg_ictp_ii_c_v = np.nanmean(freq_reg_ictp_ii_cv, axis=0)
int_reg_ictp_ii_c_v  = np.nanmean(int_reg_ictp_ii_cv, axis=0)

mean_wrf_ncar_c_v = np.nanmean(mean_wrf_ncar_cv, axis=0)
perc_wrf_ncar_c_v = np.nanmean(perc_wrf_ncar_cv, axis=0)
freq_wrf_ncar_c_v = np.nanmean(freq_wrf_ncar_cv, axis=0)
int_wrf_ncar_c_v  = np.nanmean(int_wrf_ncar_cv, axis=0)

mean_wrf_ucan_c_v = np.nanmean(mean_wrf_ucan_cv, axis=0)
perc_wrf_ucan_c_v = np.nanmean(perc_wrf_ucan_cv, axis=0)
freq_wrf_ucan_c_v = np.nanmean(freq_wrf_ucan_cv, axis=0)
int_wrf_ucan_c_v  = np.nanmean(int_wrf_ucan_cv, axis=0)

data_c_i = np.array([
    [mean_era5_c_i,        perc_era5_c_i,        freq_era5_c_i,        int_era5_c_i],
    [mean_reg_usp_c_i,     perc_reg_usp_c_i,     freq_reg_usp_c_i,     int_reg_usp_c_i],
    [mean_reg_ictp_c_i,    perc_reg_ictp_c_i,    freq_reg_ictp_c_i,    int_reg_ictp_c_i],
    [mean_reg_ictp_i_c_i,  perc_reg_ictp_i_c_i,  freq_reg_ictp_i_c_i,  int_reg_ictp_i_c_i],
    [mean_reg_ictp_ii_c_i, perc_reg_ictp_ii_c_i, freq_reg_ictp_ii_c_i, int_reg_ictp_ii_c_i],
    [mean_wrf_ncar_c_i,    perc_wrf_ncar_c_i,    freq_wrf_ncar_c_i,    int_wrf_ncar_c_i],
    [mean_wrf_ucan_c_i,    perc_wrf_ucan_c_i,    freq_wrf_ucan_c_i,    int_wrf_ucan_c_i]])

data_c_ii = np.array([
    [mean_era5_c_ii,        perc_era5_c_ii,        freq_era5_c_ii,        int_era5_c_ii],
    [mean_reg_usp_c_ii,     perc_reg_usp_c_ii,     freq_reg_usp_c_ii,     int_reg_usp_c_ii],
    [mean_reg_ictp_c_ii,    perc_reg_ictp_c_ii,    freq_reg_ictp_c_ii,    int_reg_ictp_c_ii],
    [mean_reg_ictp_i_c_ii,  perc_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii,  int_reg_ictp_i_c_ii],
    [mean_reg_ictp_ii_c_ii, perc_reg_ictp_ii_c_ii, freq_reg_ictp_ii_c_ii, int_reg_ictp_ii_c_ii],
    [mean_wrf_ncar_c_ii,    perc_wrf_ncar_c_ii,    freq_wrf_ncar_c_ii,    int_wrf_ncar_c_ii],
    [mean_wrf_ucan_c_ii,    perc_wrf_ucan_c_ii,    freq_wrf_ucan_c_ii,    int_wrf_ucan_c_ii]
])

data_c_iii = np.array([
    [mean_era5_c_iii,        perc_era5_c_iii,        freq_era5_c_iii,        int_era5_c_iii],
    [mean_reg_usp_c_iii,     perc_reg_usp_c_iii,     freq_reg_usp_c_iii,     int_reg_usp_c_iii],
    [mean_reg_ictp_c_iii,    perc_reg_ictp_c_iii,    freq_reg_ictp_c_iii,    int_reg_ictp_c_iii],
    [mean_reg_ictp_i_c_iii,  perc_reg_ictp_i_c_iii,  freq_reg_ictp_i_c_iii,  int_reg_ictp_i_c_iii],
    [mean_reg_ictp_ii_c_iii, perc_reg_ictp_ii_c_iii, freq_reg_ictp_ii_c_iii, int_reg_ictp_ii_c_iii],
    [mean_wrf_ncar_c_iii,    perc_wrf_ncar_c_iii,    freq_wrf_ncar_c_iii,    int_wrf_ncar_c_iii],
    [mean_wrf_ucan_c_iii,    perc_wrf_ucan_c_iii,    freq_wrf_ucan_c_iii,    int_wrf_ucan_c_iii]
])

data_c_iv = np.array([
    [mean_era5_c_iv,        perc_era5_c_iv,        freq_era5_c_iv,        int_era5_c_iv],
    [mean_reg_usp_c_iv,     perc_reg_usp_c_iv,     freq_reg_usp_c_iv,     int_reg_usp_c_iv],
    [mean_reg_ictp_c_iv,    perc_reg_ictp_c_iv,    freq_reg_ictp_c_iv,    int_reg_ictp_c_iv],
    [mean_reg_ictp_i_c_iv,  perc_reg_ictp_i_c_iv,  freq_reg_ictp_i_c_iv,  int_reg_ictp_i_c_iv],
    [mean_reg_ictp_ii_c_iv, perc_reg_ictp_ii_c_iv, freq_reg_ictp_ii_c_iv, int_reg_ictp_ii_c_iv],
    [mean_wrf_ncar_c_iv,    perc_wrf_ncar_c_iv,    freq_wrf_ncar_c_iv,    int_wrf_ncar_c_iv],
    [mean_wrf_ucan_c_iv,    perc_wrf_ucan_c_iv,    freq_wrf_ucan_c_iv,    int_wrf_ucan_c_iv]
])

data_c_v = np.array([
    [mean_era5_c_v,        perc_era5_c_v,        freq_era5_c_v,        int_era5_c_v],
    [mean_reg_usp_c_v,     perc_reg_usp_c_v,     freq_reg_usp_c_v,     int_reg_usp_c_v],
    [mean_reg_ictp_c_v,    perc_reg_ictp_c_v,    freq_reg_ictp_c_v,    int_reg_ictp_c_v],
    [mean_reg_ictp_i_c_v,  perc_reg_ictp_i_c_v,  freq_reg_ictp_i_c_v,  int_reg_ictp_i_c_v],
    [mean_reg_ictp_ii_c_v, perc_reg_ictp_ii_c_v, freq_reg_ictp_ii_c_v, int_reg_ictp_ii_c_v],
    [mean_wrf_ncar_c_v,    perc_wrf_ncar_c_v,    freq_wrf_ncar_c_v,    int_wrf_ncar_c_v],
    [mean_wrf_ucan_c_v,    perc_wrf_ucan_c_v,    freq_wrf_ucan_c_v,    int_wrf_ucan_c_v]
])

print(data_c_i)

# Plot figure
fig = plt.figure(figsize=(12, 5))
font_size = 7

labels = ['WRF-UCAN', 'WRF-NCAR', 'Reg5-UW', 'Reg5-Holt', 'Reg5-Holt3', 'Reg4', 'ERA5']
metrics = ['MEAN', 'P95', 'FREQ', 'INT']

if var == 'tas':
	cmap = 'bwr'
	mvmin = -20
	mvmax = 20
else:
	cmap = 'PRGn'
	mvmin = -50
	mvmax = 50
	
ax = fig.add_subplot(1, 5, 1)
sns.heatmap(data_c_i, ax=ax, annot=True, fmt='.1f', cmap=cmap, vmin=mvmin, vmax=mvmax, xticklabels=metrics, yticklabels=labels, cbar_kws={"orientation": "horizontal", "label": "(%)"}, linewidths=.5, linecolor='gray', square=True)
ax.set_title('CLUSTER I', fontsize=font_size)
ax.set_xlabel('METRICS', fontsize=font_size)
ax.set_ylabel('MODELS', fontsize=font_size)
ax.tick_params(axis='x', labelsize=font_size)
ax.tick_params(axis='y', labelsize=font_size)

ax = fig.add_subplot(1, 5, 2)
sns.heatmap(data_c_ii, ax=ax, annot=True, fmt='.1f', cmap=cmap, vmin=mvmin, vmax=mvmax, xticklabels=metrics, yticklabels=labels, cbar_kws={"orientation": "horizontal", "label": "(%)"}, linewidths=.5, linecolor='gray', square=True)
ax.set_title('CLUSTER II', fontsize=font_size)
ax.set_xlabel('METRICS', fontsize=font_size)
ax.tick_params(axis='x', labelsize=font_size)
ax.tick_params(axis='y', labelsize=font_size)
ax.set_yticklabels([])
ax.set_ylabel('')

ax = fig.add_subplot(1, 5, 3)
sns.heatmap(data_c_iii, ax=ax, annot=True, fmt='.1f', cmap=cmap, vmin=mvmin, vmax=mvmax, xticklabels=metrics, yticklabels=labels, cbar_kws={"orientation": "horizontal", "label": "(%)"}, linewidths=.5, linecolor='gray', square=True)
ax.set_title('CLUSTER III', fontsize=font_size)
ax.set_xlabel('METRICS', fontsize=font_size)
ax.tick_params(axis='x', labelsize=font_size)
ax.tick_params(axis='y', labelsize=font_size)
ax.set_yticklabels([])
ax.set_ylabel('')

ax = fig.add_subplot(1, 5, 4)
sns.heatmap(data_c_iv, ax=ax, annot=True, fmt='.1f', cmap=cmap, vmin=mvmin, vmax=mvmax, xticklabels=metrics, yticklabels=labels, cbar_kws={"orientation": "horizontal", "label": "(%)"}, linewidths=.5, linecolor='gray', square=True)
ax.set_title('CLUSTER IV', fontsize=font_size)
ax.set_xlabel('METRICS', fontsize=font_size)
ax.tick_params(axis='x', labelsize=font_size)
ax.tick_params(axis='y', labelsize=font_size)
ax.set_yticklabels([])
ax.set_ylabel('')

ax = fig.add_subplot(1, 5, 5)
sns.heatmap(data_c_v, ax=ax, annot=True, fmt='.1f', cmap=cmap, vmin=mvmin, vmax=mvmax, xticklabels=metrics, yticklabels=labels, cbar_kws={"orientation": "horizontal", "label": "(%)"}, linewidths=.5, linecolor='gray', square=True)
ax.set_title('CLUSTER V', fontsize=font_size)
ax.set_xlabel('METRICS', fontsize=font_size)
ax.set_ylabel('MODELS', fontsize=font_size)
ax.tick_params(axis='x', labelsize=font_size)
ax.tick_params(axis='y', labelsize=font_size)
ax.set_yticklabels([])
ax.set_ylabel('')

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_heatmap_stats_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
