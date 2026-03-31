# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot boxplot"

import os
import argparse
import numpy as np
import xarray as xr
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
			d_i = d_i.values
			mean_i.append(d_i)

			# Reading era5 
			d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
			d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.values 
			mean_ii.append(d_ii)
			
			# Reading regcm usp
			d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.values 
			mean_iii.append(d_iii)

			# Reading regcm ictp 
			d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.values 
			mean_iv.append(d_iv)
	
			# Reading regcm ictp pbl 1
			d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.values 
			mean_v.append(d_v)

			# Reading regcm ictp pbl 2
			d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.values 
			mean_vi.append(d_vi)

			# Reading wrf ncar
			d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.values 
			mean_vii.append(d_vii)
		
			# Reading wrf ucan
			d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_viii = d_viii.values 
			mean_viii.append(d_viii)
		
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii
	
					
def compute_stats(stations_data, threshold=None, valid_range=(-50, 60)):

    mean_list = []
    p95_list = []
    freq_list = []
    intensity_list = []

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
            intensity_list.append(float('nan'))
            continue

        # Annual mean
        mean_val = float(np.mean(da))
        mean_list.append(mean_val)

        # 95th percentile
        p95_val = float(np.percentile(da, 95))
        p95_list.append(p95_val)

        # Threshold for freq/intensity
        thr = threshold if threshold is not None else p95_val

        # Frequency (% of days above threshold)
        freq_val = float(np.sum(da > thr) / len(da) * 100)  # percentage	
        freq_list.append(freq_val)

        # Intensity (mean of values above threshold)
        intensity_val = float(np.mean(da[da > thr])) if np.any(da > thr) else 0.0
        intensity_list.append(intensity_val)

    return mean_list, p95_list, freq_list, intensity_list


def compute_relative_bias(obs_stats, sim_stats):

    mean_obs, perc_obs, freq_obs, int_obs =  [np.array(x, dtype=float) for x in obs_stats]
    mean_sim, perc_sim, freq_sim, int_sim =  [np.array(x, dtype=float) for x in sim_stats]
    
    # Avoid division by zero
    eps = 1e-10
    
    bias_mean = (mean_sim - mean_obs) / (mean_obs + eps) * 100
    bias_perc = (perc_sim - perc_obs) / (perc_obs + eps) * 100
    bias_freq = (freq_sim - freq_obs) / (freq_obs + eps) * 100
    bias_int = (int_sim - int_obs) / (int_obs + eps) * 100
    
    return bias_mean, bias_perc, bias_freq, bias_int
    

# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x, clim_viii_x = import_inmet()			

inmet_smn    = clim_i_x 
era5         = clim_ii_x 
reg_ictp     = clim_iii_x 
reg_ictp_i_  = clim_iv_x 
reg_ictp_ii_ = clim_v_x 
reg_usp      = clim_vi_x 
wrf_ncar     = clim_vii_x 
wrf_ucan     = clim_viii_x 

list_hc = [1, 2, 3, 2, 0, 1, 1, 0, 2, 2, 0, 3, 0, 2, 3, 0, 1, 2, 0, 3, 0, 4, 2, 4, 3, 1, 4, 2, 4, 2, 2, 2, 1, 2, 4, 2, 2, 3, 2, 4, 4, 4, 0, 2, 4, 3, 2, 0, 0, 0, 3, 2, 2, 2, 1, 2, 4, 1, 4, 3, 4, 3, 0, 2, 0, 3, 2, 3, 2, 4, 0, 1, 4, 2, 4, 4, 0, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 3, 2, 0, 0, 0, 4, 2, 3, 2, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 4, 0, 0, 4, 0, 4, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 2, 4, 3, 1, 4, 1, 2, 1, 1, 1, 4, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 2, 2, 2, 2, 2]
list_hc = list_hc[:len(inmet_smn)]

print(len(inmet_smn))
print(len(era5))
print(len(reg_usp))
print(len(reg_ictp))
print(len(reg_ictp_i_))
print(len(reg_ictp_ii_))
print(len(wrf_ncar))
print(len(wrf_ucan))
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

inmet_smn_i,   inmet_smn_ii,   inmet_smn_iii,   inmet_smn_iv,   inmet_smn_v   = [], [], [], [], []
era5_i,        era5_ii,        era5_iii,        era5_iv,        era5_v        = [], [], [], [], []
reg_usp_i,     reg_usp_ii,     reg_usp_iii,     reg_usp_iv,     reg_usp_v     = [], [], [], [], []
reg_ictp_i,    reg_ictp_ii,    reg_ictp_iii,    reg_ictp_iv,    reg_ictp_v    = [], [], [], [], []
reg_ictp_i_i,  reg_ictp_i_ii,  reg_ictp_i_iii,  reg_ictp_i_iv,  reg_ictp_i_v  = [], [], [], [], []
reg_ictp_ii_i, reg_ictp_ii_ii, reg_ictp_ii_iii, reg_ictp_ii_iv, reg_ictp_ii_v = [], [], [], [], []
wrf_ncar_i,    wrf_ncar_ii,    wrf_ncar_iii,    wrf_ncar_iv,    wrf_ncar_v    = [], [], [], [], []
wrf_ucan_i,    wrf_ucan_ii,    wrf_ucan_iii,    wrf_ucan_iv,    wrf_ucan_v    = [], [], [], [], []

for c_i in count_i:
	reg_usp_i.append(reg_usp[c_i])
	reg_ictp_i.append(reg_ictp[c_i])
	reg_ictp_i_i.append(reg_ictp_i_[c_i])
	reg_ictp_ii_i.append(reg_ictp_ii_[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])

for c_ii in count_ii:
	reg_usp_ii.append(reg_usp[c_ii])
	reg_ictp_ii.append(reg_ictp[c_ii])
	reg_ictp_i_ii.append(reg_ictp_i_[c_ii])
	reg_ictp_ii_ii.append(reg_ictp_ii_[c_ii])
	wrf_ncar_ii.append(wrf_ncar[c_ii])
	wrf_ucan_ii.append(wrf_ucan[c_ii])
	inmet_smn_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	
for c_iii in count_iii:
	reg_usp_iii.append(reg_usp[c_iii])
	reg_ictp_iii.append(reg_ictp[c_iii])
	reg_ictp_i_iii.append(reg_ictp_i_[c_iii])
	reg_ictp_ii_iii.append(reg_ictp_ii_[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	
for c_iv in count_iv:
	reg_usp_iv.append(reg_usp[c_iv])
	reg_ictp_iv.append(reg_ictp[c_iv])
	reg_ictp_i_iv.append(reg_ictp_i_[c_iv])
	reg_ictp_ii_iv.append(reg_ictp_ii_[c_iv])
	wrf_ncar_iv.append(wrf_ncar[c_iv])
	wrf_ucan_iv.append(wrf_ucan[c_iv])
	inmet_smn_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	
for c_v in count_v:
	reg_usp_v.append(reg_usp[c_v])
	reg_ictp_v.append(reg_ictp[c_v])
	reg_ictp_i_v.append(reg_ictp_i_[c_v])
	reg_ictp_ii_v.append(reg_ictp_ii_[c_v])
	wrf_ncar_v.append(wrf_ncar[c_v])
	wrf_ucan_v.append(wrf_ucan[c_v])
	inmet_smn_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])

# Group I
inmet_smn_c_i   = np.nanmean(inmet_smn_i, axis=0)
era5_c_i        = np.nanmean(era5_i, axis=0)
reg_usp_c_i     = np.nanmean(reg_usp_i, axis=0)
reg_ictp_c_i    = np.nanmean(reg_ictp_i, axis=0)
reg_ictp_i_c_i  = np.nanmean(reg_ictp_i_i, axis=0)
reg_ictp_ii_c_i = np.nanmean(reg_ictp_ii_i, axis=0)
wrf_ncar_c_i    = np.nanmean(wrf_ncar_i, axis=0)
wrf_ucan_c_i    = np.nanmean(wrf_ucan_i, axis=0)

mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i   = compute_stats(inmet_smn_c_i)
mean_era5_c_i,        perc_era5_c_i,        freq_era5_c_i,        int_era5_c_i        = compute_stats(era5_c_i)
mean_reg_usp_c_i,     perc_reg_usp_c_i,     freq_reg_usp_c_i,     int_reg_usp_c_i    = compute_stats(reg_usp_c_i)
mean_reg_ictp_c_i,    perc_reg_ictp_c_i,    freq_reg_ictp_c_i,    int_reg_ictp_c_i    = compute_stats(reg_ictp_c_i)
mean_reg_ictp_i_c_i,  perc_reg_ictp_i_c_i,  freq_reg_ictp_i_c_i,  int_reg_ictp_i_c_i  = compute_stats(reg_ictp_i_c_i)
mean_reg_ictp_ii_c_i, perc_reg_ictp_ii_c_i, freq_reg_ictp_ii_c_i, int_reg_ictp_ii_c_i = compute_stats(reg_ictp_ii_c_i)
mean_wrf_ncar_c_i,    perc_wrf_ncar_c_i,    freq_wrf_ncar_c_i,    int_wrf_ncar_c_i    = compute_stats(wrf_ncar_c_i)
mean_wrf_ucan_c_i,    perc_wrf_ucan_c_i,    freq_wrf_ucan_c_i,    int_wrf_ucan_c_i    = compute_stats(wrf_ucan_c_i) 

bias_era5_mean_c_i, bias_era5_perc_c_i, bias_era5_freq_c_i, bias_era5_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_era5_c_i,        perc_era5_c_i,        freq_era5_c_i,        int_era5_c_i))

bias_reg_usp_mean_c_i, bias_reg_usp_perc_c_i, bias_reg_usp_freq_c_i, bias_reg_usp_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_reg_usp_c_i,     perc_reg_usp_c_i,     freq_reg_usp_c_i,     int_reg_usp_c_i))

bias_reg_ictp_mean_c_i, bias_reg_ictp_perc_c_i, bias_reg_ictp_freq_c_i, bias_reg_ictp_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_reg_ictp_c_i,    perc_reg_ictp_c_i,    freq_reg_ictp_c_i,    int_reg_ictp_c_i))

bias_reg_ictp_i_mean_c_i, bias_reg_ictp_i_perc_c_i, bias_reg_ictp_i_freq_c_i, bias_reg_ictp_i_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_reg_ictp_i_c_i,  perc_reg_ictp_i_c_i,  freq_reg_ictp_i_c_i,  int_reg_ictp_i_c_i))

bias_reg_ictp_ii_mean_c_i, bias_reg_ictp_ii_perc_c_i, bias_reg_ictp_ii_freq_c_i, bias_reg_ictp_ii_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_reg_ictp_ii_c_i, perc_reg_ictp_ii_c_i, freq_reg_ictp_ii_c_i, int_reg_ictp_ii_c_i))

bias_wrf_ncar_mean_c_i, bias_wrf_ncar_perc_c_i, bias_wrf_ncar_freq_c_i, bias_wrf_ncar_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_wrf_ncar_c_i,    perc_wrf_ncar_c_i,    freq_wrf_ncar_c_i,    int_wrf_ncar_c_i))

bias_wrf_ucan_mean_c_i, bias_wrf_ucan_perc_c_i, bias_wrf_ucan_freq_c_i, bias_wrf_ucan_int_c_i = compute_relative_bias(
    (mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i),
    (mean_wrf_ucan_c_i,    perc_wrf_ucan_c_i,    freq_wrf_ucan_c_i,    int_wrf_ucan_c_i))

print(bias_era5_mean_c_i)

# Group II
inmet_smn_c_ii   = np.nanmean(inmet_smn_ii, axis=0)
era5_c_ii        = np.nanmean(era5_ii, axis=0)
reg_usp_c_ii     = np.nanmean(reg_usp_ii, axis=0)
reg_ictp_c_ii    = np.nanmean(reg_ictp_ii, axis=0)
reg_ictp_i_c_ii  = np.nanmean(reg_ictp_i_ii, axis=0)
reg_ictp_ii_c_ii = np.nanmean(reg_ictp_ii_ii, axis=0)
wrf_ncar_c_ii    = np.nanmean(wrf_ncar_ii, axis=0)
wrf_ucan_c_ii    = np.nanmean(wrf_ucan_ii, axis=0)

mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii   = compute_stats(inmet_smn_c_ii)
mean_era5_c_ii,        perc_era5_c_ii,        freq_era5_c_ii,        int_era5_c_ii        = compute_stats(era5_c_ii)
mean_reg_usp_c_ii,     perc_reg_usp_c_ii,     freq_reg_usp_c_ii,     int_reg_usp_c_ii    = compute_stats(reg_usp_c_ii)
mean_reg_ictp_c_ii,    perc_reg_ictp_c_ii,    freq_reg_ictp_c_ii,    int_reg_ictp_c_ii    = compute_stats(reg_ictp_c_ii)
mean_reg_ictp_i_c_ii,  perc_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii,  int_reg_ictp_i_c_ii  = compute_stats(reg_ictp_i_c_ii)
mean_reg_ictp_ii_c_ii, perc_reg_ictp_ii_c_ii, freq_reg_ictp_ii_c_ii, int_reg_ictp_ii_c_ii = compute_stats(reg_ictp_ii_c_ii)
mean_wrf_ncar_c_ii,    perc_wrf_ncar_c_ii,    freq_wrf_ncar_c_ii,    int_wrf_ncar_c_ii    = compute_stats(wrf_ncar_c_ii)
mean_wrf_ucan_c_ii,    perc_wrf_ucan_c_ii,    freq_wrf_ucan_c_ii,    int_wrf_ucan_c_ii    = compute_stats(wrf_ucan_c_ii) 
   
bias_era5_mean_c_ii, bias_era5_perc_c_ii, bias_era5_freq_c_ii, bias_era5_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_era5_c_ii,        perc_era5_c_ii,        freq_era5_c_ii,        int_era5_c_ii))

bias_reg_usp_mean_c_ii, bias_reg_usp_perc_c_ii, bias_reg_usp_freq_c_ii, bias_reg_usp_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_usp_c_ii,     perc_reg_usp_c_ii,     freq_reg_usp_c_ii,     int_reg_usp_c_ii))

bias_reg_ictp_mean_c_ii, bias_reg_ictp_perc_c_ii, bias_reg_ictp_freq_c_ii, bias_reg_ictp_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_ictp_c_ii,    perc_reg_ictp_c_ii,    freq_reg_ictp_c_ii,    int_reg_ictp_c_ii))

bias_reg_ictp_i_mean_c_ii, bias_reg_ictp_i_perc_c_ii, bias_reg_ictp_i_freq_c_ii, bias_reg_ictp_i_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_ictp_i_c_ii,  perc_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii,  int_reg_ictp_i_c_ii))

bias_reg_ictp_ii_mean_c_ii, bias_reg_ictp_ii_perc_c_ii, bias_reg_ictp_ii_freq_c_ii, bias_reg_ictp_ii_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_ictp_ii_c_ii, perc_reg_ictp_ii_c_ii, freq_reg_ictp_ii_c_ii, int_reg_ictp_ii_c_ii))

bias_wrf_ncar_mean_c_ii, bias_wrf_ncar_perc_c_ii, bias_wrf_ncar_freq_c_ii, bias_wrf_ncar_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_wrf_ncar_c_ii,    perc_wrf_ncar_c_ii,    freq_wrf_ncar_c_ii,    int_wrf_ncar_c_ii))

bias_wrf_ucan_mean_c_ii, bias_wrf_ucan_perc_c_ii, bias_wrf_ucan_freq_c_ii, bias_wrf_ucan_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_wrf_ucan_c_ii,    perc_wrf_ucan_c_ii,    freq_wrf_ucan_c_ii,    int_wrf_ucan_c_ii))
   
# Group III
inmet_smn_c_iii   = np.nanmean(inmet_smn_iii, axis=0)
era5_c_iii        = np.nanmean(era5_iii, axis=0)
reg_usp_c_iii     = np.nanmean(reg_usp_iii, axis=0)
reg_ictp_c_iii    = np.nanmean(reg_ictp_iii, axis=0)
reg_ictp_i_c_iii  = np.nanmean(reg_ictp_i_iii, axis=0)
reg_ictp_ii_c_iii = np.nanmean(reg_ictp_ii_iii, axis=0)
wrf_ncar_c_iii    = np.nanmean(wrf_ncar_iii, axis=0)
wrf_ucan_c_iii    = np.nanmean(wrf_ucan_iii, axis=0)

mean_inmet_smn_c_iii,   perc_inmet_smn_c_iii,   freq_inmet_smn_c_iii,   int_inmet_smn_c_iii   = compute_stats(inmet_smn_c_iii)
mean_era5_c_iii,        perc_era5_c_iii,        freq_era5_c_iii,        int_era5_c_iii        = compute_stats(era5_c_iii)
mean_reg_usp_c_iii,     perc_reg_usp_c_iii,     freq_reg_usp_c_iii,     int_reg_usp_c_iii    = compute_stats(reg_usp_c_iii)
mean_reg_ictp_c_iii,    perc_reg_ictp_c_iii,    freq_reg_ictp_c_iii,    int_reg_ictp_c_iii    = compute_stats(reg_ictp_c_iii)
mean_reg_ictp_i_c_iii,  perc_reg_ictp_i_c_iii,  freq_reg_ictp_i_c_iii,  int_reg_ictp_i_c_iii  = compute_stats(reg_ictp_i_c_iii)
mean_reg_ictp_ii_c_iii, perc_reg_ictp_ii_c_iii, freq_reg_ictp_ii_c_iii, int_reg_ictp_ii_c_iii = compute_stats(reg_ictp_ii_c_iii)
mean_wrf_ncar_c_iii,    perc_wrf_ncar_c_iii,    freq_wrf_ncar_c_iii,    int_wrf_ncar_c_iii    = compute_stats(wrf_ncar_c_iii)
mean_wrf_ucan_c_iii,    perc_wrf_ucan_c_iii,    freq_wrf_ucan_c_iii,    int_wrf_ucan_c_iii    = compute_stats(wrf_ucan_c_iii) 

bias_era5_mean_c_iii, bias_era5_perc_c_iii, bias_era5_freq_c_iii, bias_era5_int_c_iii = compute_relative_bias(
    (mean_inmet_smn_c_iii,   perc_inmet_smn_c_iii,   freq_inmet_smn_c_iii,   int_inmet_smn_c_iii),
    (mean_era5_c_iii,        perc_era5_c_iii,        freq_era5_c_iii,        int_era5_c_iii))

bias_reg_usp_mean_c_iii, bias_reg_usp_perc_c_ii, bias_reg_usp_freq_c_ii, bias_reg_usp_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_iii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_usp_c_iii,     perc_reg_usp_c_ii,     freq_reg_usp_c_ii,     int_reg_usp_c_ii))

bias_reg_ictp_mean_c_iii, bias_reg_ictp_perc_c_iii, bias_reg_ictp_freq_c_iii, bias_reg_ictp_int_c_iii = compute_relative_bias(
    (mean_inmet_smn_c_iii,   perc_inmet_smn_c_iii,   freq_inmet_smn_c_iii,   int_inmet_smn_c_iii),
    (mean_reg_ictp_c_iii,    perc_reg_ictp_c_iii,    freq_reg_ictp_c_iii,    int_reg_ictp_c_iii))

bias_reg_ictp_i_mean_c_ii, bias_reg_ictp_i_perc_c_ii, bias_reg_ictp_i_freq_c_ii, bias_reg_ictp_i_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_ictp_i_c_ii,  perc_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii,  int_reg_ictp_i_c_ii))

bias_reg_ictp_ii_mean_c_ii, bias_reg_ictp_ii_perc_c_ii, bias_reg_ictp_ii_freq_c_ii, bias_reg_ictp_ii_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_reg_ictp_ii_c_ii, perc_reg_ictp_ii_c_ii, freq_reg_ictp_ii_c_ii, int_reg_ictp_ii_c_ii))

bias_wrf_ncar_mean_c_ii, bias_wrf_ncar_perc_c_ii, bias_wrf_ncar_freq_c_ii, bias_wrf_ncar_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_wrf_ncar_c_ii,    perc_wrf_ncar_c_ii,    freq_wrf_ncar_c_ii,    int_wrf_ncar_c_ii))

bias_wrf_ucan_mean_c_ii, bias_wrf_ucan_perc_c_ii, bias_wrf_ucan_freq_c_ii, bias_wrf_ucan_int_c_ii = compute_relative_bias(
    (mean_inmet_smn_c_ii,   perc_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   int_inmet_smn_c_ii),
    (mean_wrf_ucan_c_ii,    perc_wrf_ucan_c_ii,    freq_wrf_ucan_c_ii,    int_wrf_ucan_c_ii))
    
# Group IV
inmet_smn_c_iv   = np.nanmean(inmet_smn_iv, axis=0)
era5_c_iv        = np.nanmean(era5_iv, axis=0)
reg_usp_c_iv     = np.nanmean(reg_usp_iv, axis=0)
reg_ictp_c_iv    = np.nanmean(reg_ictp_iv, axis=0)
reg_ictp_i_c_iv  = np.nanmean(reg_ictp_i_iv, axis=0)
reg_ictp_ii_c_iv = np.nanmean(reg_ictp_ii_iv, axis=0)
wrf_ncar_c_iv    = np.nanmean(wrf_ncar_iv, axis=0)
wrf_ucan_c_iv    = np.nanmean(wrf_ucan_iv, axis=0)

mean_inmet_smn_c_iv,   perc_inmet_smn_c_iv,   freq_inmet_smn_c_iv,   int_inmet_smn_c_iv   = compute_stats(inmet_smn_c_iv)
mean_era5_c_iv,        perc_era5_c_iv,        freq_era5_c_iv,        int_era5_c_iv        = compute_stats(era5_c_iv)
mean_reg_usp_c_iv,     perc_reg_usp_c_iv,     freq_reg_usp_c_iv,     int_reg_usp_c_iv    = compute_stats(reg_usp_c_iv)
mean_reg_ictp_c_iv,    perc_reg_ictp_c_iv,    freq_reg_ictp_c_iv,    int_reg_ictp_c_iv    = compute_stats(reg_ictp_c_iv)
mean_reg_ictp_i_c_iv,  perc_reg_ictp_i_c_iv,  freq_reg_ictp_i_c_iv,  int_reg_ictp_i_c_iv  = compute_stats(reg_ictp_i_c_iv)
mean_reg_ictp_ii_c_iv, perc_reg_ictp_ii_c_iv, freq_reg_ictp_ii_c_iv, int_reg_ictp_ii_c_iv = compute_stats(reg_ictp_ii_c_iv)
mean_wrf_ncar_c_iv,    perc_wrf_ncar_c_iv,    freq_wrf_ncar_c_iv,    int_wrf_ncar_c_iv    = compute_stats(wrf_ncar_c_iv)
mean_wrf_ucan_c_iv,    perc_wrf_ucan_c_iv,    freq_wrf_ucan_c_iv,    int_wrf_ucan_c_iv    = compute_stats(wrf_ucan_c_iv) 

# Group V
inmet_smn_c_v   = np.nanmean(inmet_smn_v, axis=0)
era5_c_v        = np.nanmean(era5_v, axis=0)
reg_usp_c_v     = np.nanmean(reg_usp_v, axis=0)
reg_ictp_c_v    = np.nanmean(reg_ictp_v, axis=0)
reg_ictp_i_c_v  = np.nanmean(reg_ictp_i_v, axis=0)
reg_ictp_ii_c_v = np.nanmean(reg_ictp_ii_v, axis=0)
wrf_ncar_c_v    = np.nanmean(wrf_ncar_v, axis=0)
wrf_ucan_c_v    = np.nanmean(wrf_ucan_v, axis=0)

mean_inmet_smn_c_v,   perc_inmet_smn_c_v,   freq_inmet_smn_c_v,   int_inmet_smn_c_v   = compute_stats(inmet_smn_c_v)
mean_era5_c_v,        perc_era5_c_v,        freq_era5_c_v,        int_era5_c_v        = compute_stats(era5_c_v)
mean_reg_usp_c_v,     perc_reg_usp_c_v,     freq_reg_usp_c_v,     int_reg_usp_c_v     = compute_stats(reg_usp_c_v)
mean_reg_ictp_c_v,    perc_reg_ictp_c_v,    freq_reg_ictp_c_v,    int_reg_ictp_c_v    = compute_stats(reg_ictp_c_v)
mean_reg_ictp_i_c_v,  perc_reg_ictp_i_c_v,  freq_reg_ictp_i_c_v,  int_reg_ictp_i_c_v  = compute_stats(reg_ictp_i_c_v)
mean_reg_ictp_ii_c_v, perc_reg_ictp_ii_c_v, freq_reg_ictp_ii_c_v, int_reg_ictp_ii_c_v = compute_stats(reg_ictp_ii_c_v)
mean_wrf_ncar_c_v,    perc_wrf_ncar_c_v,    freq_wrf_ncar_c_v,    int_wrf_ncar_c_v    = compute_stats(wrf_ncar_c_v)
mean_wrf_ucan_c_v,    perc_wrf_ucan_c_v,    freq_wrf_ucan_c_v,    int_wrf_ucan_c_v    = compute_stats(wrf_ucan_c_v) 
    
# Plot figure
fig = plt.figure(figsize=(8, 18))
font_size = 8

metrics = ['Mean', 'Percentile', 'Frequency', 'Intensity']

if var == 'tas':
	legend = 'Air emperature 2m (°C)'
	yvmin = 0
	yvmax = 40
	yvmax_ = 44
	yint_ = 4
else:
	legend = 'Wind speed 10m (m s⁻¹)'
	yvmin = 0
	yvmax =10
	yvmax_ = 11
	yint_ = 1
	
ax1 = fig.add_subplot(1, 5, 1)
plt.boxplot(bias_era5_mean_c_i)

ax2 = fig.add_subplot(1, 5, 2)
plt.boxplot(bias_era5_mean_c_i)

ax3 = fig.add_subplot(1, 5, 3)
plt.boxplot(bias_era5_mean_c_i)

ax4 = fig.add_subplot(1, 5, 4)
plt.boxplot(bias_era5_mean_c_i)

ax5 = fig.add_subplot(1, 5, 5)
plt.boxplot(bias_era5_mean_c_i)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_boxplot_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
