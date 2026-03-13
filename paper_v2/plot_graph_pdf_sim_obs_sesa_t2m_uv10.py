
# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot PDFs"

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

if var == 'tas':
	percentile = 75
else:
	percentile = 90
	
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
			
			# Reading regcm ictp 
			d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.values 
			mean_iii.append(d_iii)
	
			# Reading regcm ictp pbl 1
			d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.values 
			mean_iv.append(d_iv)

			# Reading regcm ictp pbl 2
			d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.values 
			mean_v.append(d_v)

			# Reading regcm usp
			d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
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
	
					
def mask_like(reference, target):

	reference = np.asarray(reference)
	masked_target = np.asarray(target, dtype=float).copy()
	masked_target[np.isnan(reference)] = np.nan
	
	return masked_target
	

def compute_pdf(ts_hourly, nbins=10000, max_kde_points=200000):

    ts = np.asarray(ts_hourly).ravel()

    # Remove NaN and inf
    ts[~np.isfinite(ts)] = np.nan
    ts[np.abs(ts) > 1e20] = np.nan

    # Filters
    if var == 'tas':        # Celsius
        ts = ts[(ts > -80) & (ts < 60)]

    elif var == 'sfcWind':  # m/s
        ts = ts[(ts >= 0) & (ts < 80)]

    if ts.size < 10:
        return np.array([]), np.array([]), np.nan, np.nan

    # Percentile and max
    perc = np.nanpercentile(ts, percentile)
    max_ = np.nanmax(ts)

    # KDE
    if ts.size > max_kde_points:
        idx = np.random.choice(ts.size, max_kde_points, replace=False)
        ts_kde = ts[idx]
    else:
        ts_kde = ts

    kde = gaussian_kde(ts_kde, bw_method="scott")

    x = np.linspace(ts.min(), ts.max(), nbins)
    pdf = kde(x)

    return x, pdf, max_, perc
    

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
inmet_smn_c_i   = np.concatenate(inmet_smn_i)
inmet_smn_c_ii  = np.concatenate(inmet_smn_ii)
inmet_smn_c_iii = np.concatenate(inmet_smn_iii)
inmet_smn_c_iv  = np.concatenate(inmet_smn_iv)
inmet_smn_c_v   = np.concatenate(inmet_smn_v)

# Group II
era5_c_i_   = np.concatenate(era5_i)
era5_c_ii_  = np.concatenate(era5_ii)
era5_c_iii_ = np.concatenate(era5_iii)
era5_c_iv_  = np.concatenate(era5_iv)
era5_c_v_   = np.concatenate(era5_v)

era5_c_i = mask_like(inmet_smn_c_i, era5_c_i_)
era5_c_ii = mask_like(inmet_smn_c_ii, era5_c_ii_)
era5_c_iii = mask_like(inmet_smn_c_iii, era5_c_iii_)
era5_c_iv = mask_like(inmet_smn_c_iv, era5_c_iv_)
era5_c_v = mask_like(inmet_smn_c_v, era5_c_v_)

# Group III
reg_usp_c_i_   = np.concatenate(reg_usp_i)
reg_usp_c_ii_  = np.concatenate(reg_usp_ii)
reg_usp_c_iii_ = np.concatenate(reg_usp_iii)
reg_usp_c_iv_  = np.concatenate(reg_usp_iv)
reg_usp_c_v_   = np.concatenate(reg_usp_v)

reg_usp_c_i = mask_like(inmet_smn_c_i, reg_usp_c_i_)
reg_usp_c_ii = mask_like(inmet_smn_c_ii, reg_usp_c_ii_)
reg_usp_c_iii = mask_like(inmet_smn_c_iii, reg_usp_c_iii_)
reg_usp_c_iv = mask_like(inmet_smn_c_iv, reg_usp_c_iv_)
reg_usp_c_v = mask_like(inmet_smn_c_v, reg_usp_c_v_)

# Group IV
reg_ictp_c_i_   = np.concatenate(reg_ictp_i)
reg_ictp_c_ii_  = np.concatenate(reg_ictp_ii)
reg_ictp_c_iii_ = np.concatenate(reg_ictp_iii)
reg_ictp_c_iv_  = np.concatenate(reg_ictp_iv)
reg_ictp_c_v_   = np.concatenate(reg_ictp_v)

reg_ictp_c_i = mask_like(inmet_smn_c_i, reg_ictp_c_i_)
reg_ictp_c_ii = mask_like(inmet_smn_c_ii, reg_ictp_c_ii_)
reg_ictp_c_iii = mask_like(inmet_smn_c_iii, reg_ictp_c_iii_)
reg_ictp_c_iv = mask_like(inmet_smn_c_iv, reg_ictp_c_iv_)
reg_ictp_c_v = mask_like(inmet_smn_c_v, reg_ictp_c_v_)

# Group V
reg_ictp_i_c_i_   = np.concatenate(reg_ictp_i_i)
reg_ictp_i_c_ii_  = np.concatenate(reg_ictp_i_ii)
reg_ictp_i_c_iii_ = np.concatenate(reg_ictp_i_iii)
reg_ictp_i_c_iv_  = np.concatenate(reg_ictp_i_iv)
reg_ictp_i_c_v_  = np.concatenate(reg_ictp_i_v)

reg_ictp_i_c_i = mask_like(inmet_smn_c_i, reg_ictp_i_c_i_)
reg_ictp_i_c_ii = mask_like(inmet_smn_c_ii, reg_ictp_i_c_ii_)
reg_ictp_i_c_iii = mask_like(inmet_smn_c_iii, reg_ictp_i_c_iii_)
reg_ictp_i_c_iv = mask_like(inmet_smn_c_iv, reg_ictp_i_c_iv_)
reg_ictp_i_c_v = mask_like(inmet_smn_c_v, reg_ictp_i_c_v_)

# Group VI
reg_ictp_ii_c_i_   = np.concatenate(reg_ictp_ii_i)
reg_ictp_ii_c_ii_  = np.concatenate(reg_ictp_ii_ii)
reg_ictp_ii_c_iii_ = np.concatenate(reg_ictp_ii_iii)
reg_ictp_ii_c_iv_  = np.concatenate(reg_ictp_ii_iv)
reg_ictp_ii_c_v_   = np.concatenate(reg_ictp_ii_v)

reg_ictp_ii_c_i = mask_like(inmet_smn_c_i, reg_ictp_ii_c_i_)
reg_ictp_ii_c_ii = mask_like(inmet_smn_c_ii, reg_ictp_ii_c_ii_)
reg_ictp_ii_c_iii = mask_like(inmet_smn_c_iii, reg_ictp_ii_c_iii_)
reg_ictp_ii_c_iv = mask_like(inmet_smn_c_iv, reg_ictp_ii_c_iv_)
reg_ictp_ii_c_v = mask_like(inmet_smn_c_v, reg_ictp_ii_c_v_)

# Group VII
wrf_ncar_c_i_   = np.concatenate(wrf_ncar_i)
wrf_ncar_c_ii_  = np.concatenate(wrf_ncar_ii)
wrf_ncar_c_iii_ = np.concatenate(wrf_ncar_iii)
wrf_ncar_c_iv_  = np.concatenate(wrf_ncar_iv)
wrf_ncar_c_v_   = np.concatenate(wrf_ncar_v)

wrf_ncar_c_i = mask_like(inmet_smn_c_i, wrf_ncar_c_i_)
wrf_ncar_c_ii = mask_like(inmet_smn_c_ii, wrf_ncar_c_ii_)
wrf_ncar_c_iii = mask_like(inmet_smn_c_iii, wrf_ncar_c_iii_)
wrf_ncar_c_iv = mask_like(inmet_smn_c_iv, wrf_ncar_c_iv_)
wrf_ncar_c_v = mask_like(inmet_smn_c_v, wrf_ncar_c_v_)

# Group VIII
wrf_ucan_c_i_   = np.concatenate(wrf_ucan_i)
wrf_ucan_c_ii_  = np.concatenate(wrf_ucan_ii)
wrf_ucan_c_iii_ = np.concatenate(wrf_ucan_iii)
wrf_ucan_c_iv_  = np.concatenate(wrf_ucan_iv)
wrf_ucan_c_v_   = np.concatenate(wrf_ucan_v)

wrf_ucan_c_i = mask_like(inmet_smn_c_i, wrf_ucan_c_i_)
wrf_ucan_c_ii = mask_like(inmet_smn_c_ii, wrf_ucan_c_ii_)
wrf_ucan_c_iii = mask_like(inmet_smn_c_iii, wrf_ucan_c_iii_)
wrf_ucan_c_iv = mask_like(inmet_smn_c_iv, wrf_ucan_c_iv_)
wrf_ucan_c_v = mask_like(inmet_smn_c_v,wrf_ucan_c_v_)

x_inmet_smn_c_i,   pdf_inmet_smn_c_i,   max_inmet_smn_c_i,   perc_inmet_smn_c_i   = compute_pdf(inmet_smn_c_i)
x_inmet_smn_c_ii,  pdf_inmet_smn_c_ii,  max_inmet_smn_c_ii,  perc_inmet_smn_c_ii  = compute_pdf(inmet_smn_c_ii)
x_inmet_smn_c_iii, pdf_inmet_smn_c_iii, max_inmet_smn_c_iii, perc_inmet_smn_c_iii = compute_pdf(inmet_smn_c_iii)
x_inmet_smn_c_iv,  pdf_inmet_smn_c_iv,  max_inmet_smn_c_iv,  perc_inmet_smn_c_iv  = compute_pdf(inmet_smn_c_iv)
x_inmet_smn_c_v,   pdf_inmet_smn_c_v,   max_inmet_smn_c_v,   perc_inmet_smn_c_v   = compute_pdf(inmet_smn_c_v)

x_era5_c_i,   pdf_era5_c_i,   max_era5_c_i,   perc_era5_c_i   = compute_pdf(era5_c_i)
x_era5_c_ii,  pdf_era5_c_ii,  max_era5_c_ii,  perc_era5_c_ii  = compute_pdf(era5_c_ii)
x_era5_c_iii, pdf_era5_c_iii, max_era5_c_iii, perc_era5_c_iii = compute_pdf(era5_c_iii)
x_era5_c_iv,  pdf_era5_c_iv,  max_era5_c_iv,  perc_era5_c_iv  = compute_pdf(era5_c_iv)
x_era5_c_v,   pdf_era5_c_v,   max_era5_c_v,   perc_era5_c_v   = compute_pdf(era5_c_v)

x_reg_usp_c_i,   pdf_reg_usp_c_i,   max_reg_usp_c_i,   perc_reg_usp_c_i   = compute_pdf(reg_usp_c_i)
x_reg_usp_c_ii,  pdf_reg_usp_c_ii,  max_reg_usp_c_ii,  perc_reg_usp_c_ii  = compute_pdf(reg_usp_c_ii)
x_reg_usp_c_iii, pdf_reg_usp_c_iii, max_reg_usp_c_iii, perc_reg_usp_c_iii = compute_pdf(reg_usp_c_iii)
x_reg_usp_c_iv,  pdf_reg_usp_c_iv,  max_reg_usp_c_iv,  perc_reg_usp_c_iv  = compute_pdf(reg_usp_c_iv)
x_reg_usp_c_v,   pdf_reg_usp_c_v,   max_reg_usp_c_v,   perc_reg_usp_c_v   = compute_pdf(reg_usp_c_v)

x_reg_ictp_c_i,   pdf_reg_ictp_c_i,   max_reg_ictp_c_i,   perc_reg_ictp_c_i   = compute_pdf(reg_ictp_c_i)
x_reg_ictp_c_ii,  pdf_reg_ictp_c_ii,  max_reg_ictp_c_ii,  perc_reg_ictp_c_ii  = compute_pdf(reg_ictp_c_ii)
x_reg_ictp_c_iii, pdf_reg_ictp_c_iii, max_reg_ictp_c_iii, perc_reg_ictp_c_iii = compute_pdf(reg_ictp_c_iii)
x_reg_ictp_c_iv,  pdf_reg_ictp_c_iv,  max_reg_ictp_c_iv,  perc_reg_ictp_c_iv  = compute_pdf(reg_ictp_c_iv)
x_reg_ictp_c_v,   pdf_reg_ictp_c_v,   max_reg_ictp_c_v,   perc_reg_ictp_c_v   = compute_pdf(reg_ictp_c_v)

x_reg_ictp_i_c_i,   pdf_reg_ictp_i_c_i,   max_reg_ictp_i_c_i,   perc_reg_ictp_i_c_i   = compute_pdf(reg_ictp_i_c_i)
x_reg_ictp_i_c_ii,  pdf_reg_ictp_i_c_ii,  max_reg_ictp_i_c_ii,  perc_reg_ictp_i_c_ii  = compute_pdf(reg_ictp_i_c_ii)
x_reg_ictp_i_c_iii, pdf_reg_ictp_i_c_iii, max_reg_ictp_i_c_iii, perc_reg_ictp_i_c_iii = compute_pdf(reg_ictp_i_c_iii)
x_reg_ictp_i_c_iv,  pdf_reg_ictp_i_c_iv,  max_reg_ictp_i_c_iv,  perc_reg_ictp_i_c_iv  = compute_pdf(reg_ictp_i_c_iv)
x_reg_ictp_i_c_v,   pdf_reg_ictp_i_c_v,   max_reg_ictp_i_c_v,   perc_reg_ictp_i_c_v   = compute_pdf(reg_ictp_i_c_v)

x_reg_ictp_ii_c_i,   pdf_reg_ictp_ii_c_i,   max_reg_ictp_ii_c_i,   perc_reg_ictp_ii_c_i   = compute_pdf(reg_ictp_ii_c_i)
x_reg_ictp_ii_c_ii,  pdf_reg_ictp_ii_c_ii,  max_reg_ictp_ii_c_ii,  perc_reg_ictp_ii_c_ii  = compute_pdf(reg_ictp_ii_c_ii)
x_reg_ictp_ii_c_iii, pdf_reg_ictp_ii_c_iii, max_reg_ictp_ii_c_iii, perc_reg_ictp_ii_c_iii = compute_pdf(reg_ictp_ii_c_iii)
x_reg_ictp_ii_c_iv,  pdf_reg_ictp_ii_c_iv,  max_reg_ictp_ii_c_iv,  perc_reg_ictp_ii_c_iv  = compute_pdf(reg_ictp_ii_c_iv)
x_reg_ictp_ii_c_v,   pdf_reg_ictp_ii_c_v,   max_reg_ictp_ii_c_v,   perc_reg_ictp_ii_c_v   = compute_pdf(reg_ictp_ii_c_v)

x_wrf_ncar_c_i,   pdf_wrf_ncar_c_i,   max_wrf_ncar_c_i, perc_wrf_ncar_c_i     = compute_pdf(wrf_ncar_c_i)
x_wrf_ncar_c_ii,  pdf_wrf_ncar_c_ii,  max_wrf_ncar_c_ii, perc_wrf_ncar_c_ii   = compute_pdf(wrf_ncar_c_ii)
x_wrf_ncar_c_iii, pdf_wrf_ncar_c_iii, max_wrf_ncar_c_iii, perc_wrf_ncar_c_iii = compute_pdf(wrf_ncar_c_iii)
x_wrf_ncar_c_iv,  pdf_wrf_ncar_c_iv,  max_wrf_ncar_c_iv, perc_wrf_ncar_c_iv   = compute_pdf(wrf_ncar_c_iv)
x_wrf_ncar_c_v,   pdf_wrf_ncar_c_v,   max_wrf_ncar_c_v, perc_wrf_ncar_c_v     = compute_pdf(wrf_ncar_c_v)

x_wrf_ucan_c_i,   pdf_wrf_ucan_c_i,   max_wrf_ucan_c_i,   perc_wrf_ucan_c_i   = compute_pdf(wrf_ucan_c_i)
x_wrf_ucan_c_ii,  pdf_wrf_ucan_c_ii,  max_wrf_ucan_c_ii,  perc_wrf_ucan_c_ii  = compute_pdf(wrf_ucan_c_ii)
x_wrf_ucan_c_iii, pdf_wrf_ucan_c_iii, max_wrf_ucan_c_iii, perc_wrf_ucan_c_iii = compute_pdf(wrf_ucan_c_iii)
x_wrf_ucan_c_iv,  pdf_wrf_ucan_c_iv,  max_wrf_ucan_c_iv,  perc_wrf_ucan_c_iv  = compute_pdf(wrf_ucan_c_iv)
x_wrf_ucan_c_v,   pdf_wrf_ucan_c_v,   max_wrf_ucan_c_v,   perc_wrf_ucan_c_v   = compute_pdf(wrf_ucan_c_v)

# Plot figure
fig = plt.figure(figsize=(8, 18))
font_size = 8

if var == 'tas':
	legend = 'Air emperature 2m (°C)'
	xvmin = 0
	xvmax = 40
	xvmax_ = 44
	xint_ = 4
	yvmin = 0
	yvmax = 0.1
	yvmax_ = 0.11
	yint_ = 0.01
	text_ = 1
	text_1 = 0.2
	text_2 = 0.3
	text_3 = 0.4
	text_4 = 0.5
	text_5 = 0.6
	text_6 = 0.7
	text_7 = 0.8
	text_8 = 0.9

else:
	legend = 'Wind speed 10m (m s⁻¹)'
	xvmin = 0
	xvmax =10
	xvmax_ = 11
	xint_ = 1
	yvmin = 0
	yvmax = 1
	yvmax_ = 1.1
	yint_ = 0.1
	text_ = 4
	text_1 = 0.2
	text_2 = 0.3
	text_3 = 0.4
	text_4 = 0.5
	text_5 = 0.6
	text_6 = 0.7
	text_7 = 0.8
	text_8 = 0.9

ax1 = fig.add_subplot(6, 2, 1)
plt.plot(x_inmet_smn_c_i,   pdf_inmet_smn_c_i,   linewidth=2, color='black',   label='INMET+SMN')
plt.plot(x_era5_c_i,        pdf_era5_c_i,        linewidth=1, color='red',     label='ERA5')
plt.plot(x_reg_usp_c_i,     pdf_reg_usp_c_i,     linewidth=1, color='blue',    label='Reg4')
plt.plot(x_reg_ictp_c_i,    pdf_reg_ictp_c_i,    linewidth=1, color='magenta', label='Reg5-Holt3')
plt.plot(x_reg_ictp_i_c_i,  pdf_reg_ictp_i_c_i,  linewidth=1, color='gray',    label='Reg5-Holt')
plt.plot(x_reg_ictp_ii_c_i, pdf_reg_ictp_ii_c_i, linewidth=1, color='brown',   label='Reg5-UW')
plt.plot(x_wrf_ncar_c_i,    pdf_wrf_ncar_c_i,    linewidth=1, color='green',   label='WRF-NCAR')
plt.plot(x_wrf_ucan_c_i,    pdf_wrf_ucan_c_i,    linewidth=1, color='orange',  label='WRF-UCAN')
plt.text(text_, text_8, 'INMET = {0} ({1})'.format(round(perc_inmet_smn_c_i, 1),      round(max_inmet_smn_c_i, 1)),   color='black',   fontsize=font_size)
plt.text(text_, text_7, 'ERA5 = {0} ({1})'.format(round(perc_era5_c_i, 1),            round(max_era5_c_i, 1)),        color='red',     fontsize=font_size)
plt.text(text_, text_6, 'Reg4 = {0} ({1})'.format(round(perc_reg_usp_c_i, 1),         round(max_reg_usp_c_i, 1)),     color='blue',    fontsize=font_size)
plt.text(text_, text_5, 'Reg5-Holt3 = {0} ({1})'.format(round(perc_reg_ictp_c_i, 1),  round(max_reg_ictp_c_i, 1)),    color='magenta', fontsize=font_size)
plt.text(text_, text_4, 'Reg5-Holt = {0} ({1})'.format(round(perc_reg_ictp_i_c_i, 1), round(max_reg_ictp_i_c_i, 1)),  color='gray',    fontsize=font_size)
plt.text(text_, text_3, 'Reg5-UW = {0} ({1})'.format(round(perc_reg_ictp_ii_c_i, 1),  round(max_reg_ictp_ii_c_i, 1)), color='brown',   fontsize=font_size)
plt.text(text_, text_2, 'WRF-NCAR = {0} ({1})'.format(round(perc_wrf_ncar_c_i, 1),    round(max_wrf_ncar_c_i, 1)),    color='green',   fontsize=font_size)
plt.text(text_, text_1, 'WRF-UCAN = {0} ({1})'.format(round(perc_wrf_ucan_c_i, 1),    round(max_wrf_ucan_c_i, 1)),    color='orange',  fontsize=font_size)
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlim(xvmin, xvmax)
plt.ylim(yvmin, yvmax)
plt.xticks(np.arange(xvmin, xvmax_, xint_), fontsize=font_size)
plt.yticks(np.arange(yvmin, yvmax_, yint_), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax2 = fig.add_subplot(6, 2, 2)
plt.plot(x_inmet_smn_c_ii,   pdf_inmet_smn_c_ii,   linewidth=2, color='black',   label='INMET+SMN')
plt.plot(x_era5_c_ii,        pdf_era5_c_ii,        linewidth=1, color='red',     label='ERA5')
plt.plot(x_reg_usp_c_ii,     pdf_reg_usp_c_ii,     linewidth=1, color='blue',    label='Reg4')
plt.plot(x_reg_ictp_c_ii,    pdf_reg_ictp_c_ii,    linewidth=1, color='magenta', label='Reg5-Holt3')
plt.plot(x_reg_ictp_i_c_ii,  pdf_reg_ictp_i_c_ii,  linewidth=1, color='gray',    label='Reg5-Holt')
plt.plot(x_reg_ictp_ii_c_ii, pdf_reg_ictp_ii_c_ii, linewidth=1, color='brown',   label='Reg5-UW')
plt.plot(x_wrf_ncar_c_ii,    pdf_wrf_ncar_c_ii,    linewidth=1, color='green',   label='WRF-NCAR')
plt.plot(x_wrf_ucan_c_ii,    pdf_wrf_ucan_c_ii,    linewidth=1, color='orange',  label='WRF-UCAN')
plt.text(text_, text_8, 'INMET = {0} ({1})'.format(round(perc_inmet_smn_c_ii, 1),      round(max_inmet_smn_c_ii, 1)),   color='black',   fontsize=font_size)
plt.text(text_, text_7, 'ERA5 = {0} ({1})'.format(round(perc_era5_c_ii, 1),            round(max_era5_c_ii, 1)),        color='red',     fontsize=font_size)
plt.text(text_, text_6, 'Reg4 = {0} ({1})'.format(round(perc_reg_usp_c_ii, 1),	     round(max_reg_usp_c_ii, 1)),     color='blue',    fontsize=font_size)
plt.text(text_, text_5, 'Reg5-Holt3 = {0} ({1})'.format(round(perc_reg_ictp_c_ii, 1),  round(max_reg_ictp_c_ii, 1)),    color='magenta', fontsize=font_size)
plt.text(text_, text_4, 'Reg5-Holt = {0} ({1})'.format(round(perc_reg_ictp_i_c_ii, 1), round(max_reg_ictp_i_c_ii, 1)),  color='gray',    fontsize=font_size)
plt.text(text_, text_3, 'Reg5-UW = {0} ({1})'.format(round(perc_reg_ictp_ii_c_ii, 1),  round(max_reg_ictp_ii_c_ii, 1)), color='brown',   fontsize=font_size)
plt.text(text_, text_2, 'WRF-NCAR = {0} ({1})'.format(round(perc_wrf_ncar_c_ii, 1),    round(max_wrf_ncar_c_ii, 1)),    color='green',   fontsize=font_size)
plt.text(text_, text_1, 'WRF-UCAN = {0} ({1})'.format(round(perc_wrf_ucan_c_ii, 1),    round(max_wrf_ucan_c_ii, 1)),    color='orange',  fontsize=font_size)
plt.title('(b) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlim(xvmin, xvmax)
plt.ylim(yvmin, yvmax)
plt.xticks(np.arange(xvmin, xvmax_, xint_), fontsize=font_size)
plt.yticks(np.arange(yvmin, yvmax_, yint_), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax3 = fig.add_subplot(6, 2, 3)
plt.plot(x_inmet_smn_c_iii,   pdf_inmet_smn_c_iii,   linewidth=2, color='black',   label='INMET+SMN')
plt.plot(x_era5_c_iii,        pdf_era5_c_iii,        linewidth=1, color='red',     label='ERA5')
plt.plot(x_reg_usp_c_iii,     pdf_reg_usp_c_iii,     linewidth=1, color='blue',    label='Reg4')
plt.plot(x_reg_ictp_c_iii,    pdf_reg_ictp_c_iii,    linewidth=1, color='magenta', label='Reg5-Holt3')
plt.plot(x_reg_ictp_i_c_iii,  pdf_reg_ictp_i_c_iii,  linewidth=1, color='gray',    label='Reg5-Holt')
plt.plot(x_reg_ictp_ii_c_iii, pdf_reg_ictp_ii_c_iii, linewidth=1, color='brown',   label='Reg5-UW')
plt.plot(x_wrf_ncar_c_iii,    pdf_wrf_ncar_c_iii,    linewidth=1, color='green',   label='WRF-NCAR')
plt.plot(x_wrf_ucan_c_iii,    pdf_wrf_ucan_c_iii,    linewidth=1, color='orange',  label='WRF-UCAN')
plt.text(text_, text_8, 'INMET = {0} ({1})'.format(round(perc_inmet_smn_c_iii, 1),      round(max_inmet_smn_c_iii, 1)),   color='black',   fontsize=font_size)
plt.text(text_, text_7, 'ERA5 = {0} ({1})'.format(round(perc_era5_c_iii, 1),            round(max_era5_c_iii, 1)),        color='red',     fontsize=font_size)
plt.text(text_, text_6, 'Reg4 = {0} ({1})'.format(round(perc_reg_usp_c_iii, 1),         round(max_reg_usp_c_iii, 1)),     color='blue',    fontsize=font_size)
plt.text(text_, text_5, 'Reg5-Holt3 = {0} ({1})'.format(round(perc_reg_ictp_c_iii, 1),  round(max_reg_ictp_c_iii, 1)),	 color='magenta', fontsize=font_size)
plt.text(text_, text_4, 'Reg5-Holt = {0} ({1})'.format(round(perc_reg_ictp_i_c_iii, 1), round(max_reg_ictp_i_c_iii, 1)),  color='gray',    fontsize=font_size)
plt.text(text_, text_3, 'Reg5-UW = {0} ({1})'.format(round(perc_reg_ictp_ii_c_iii, 1),  round(max_reg_ictp_ii_c_iii, 1)), color='brown',   fontsize=font_size)
plt.text(text_, text_2, 'WRF-NCAR = {0} ({1})'.format(round(perc_wrf_ncar_c_iii, 1),    round(max_wrf_ncar_c_iii, 1)),    color='green',   fontsize=font_size)
plt.text(text_, text_1, 'WRF-UCAN = {0} ({1})'.format(round(perc_wrf_ucan_c_iii, 1),    round(max_wrf_ucan_c_iii, 1)),    color='orange',  fontsize=font_size)
plt.title('(c) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlim(xvmin, xvmax)
plt.ylim(yvmin, yvmax)
plt.xticks(np.arange(xvmin, xvmax_, xint_), fontsize=font_size)
plt.yticks(np.arange(yvmin, yvmax_, yint_), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax4 = fig.add_subplot(6, 2, 4)
plt.plot(x_inmet_smn_c_iv,   pdf_inmet_smn_c_iv,   linewidth=2, color='black',   label='INMET+SMN')
plt.plot(x_era5_c_iv,        pdf_era5_c_iv,        linewidth=1, color='red',     label='ERA5')
plt.plot(x_reg_usp_c_iv,     pdf_reg_usp_c_iv,     linewidth=1, color='blue',    label='Reg4')
plt.plot(x_reg_ictp_c_iv,    pdf_reg_ictp_c_iv,    linewidth=1, color='magenta', label='Reg5-Holt3')
plt.plot(x_reg_ictp_i_c_iv,  pdf_reg_ictp_i_c_iv,  linewidth=1, color='gray',    label='Reg5-Holt')
plt.plot(x_reg_ictp_ii_c_iv, pdf_reg_ictp_ii_c_iv, linewidth=1, color='brown',   label='Reg5-UW')
plt.plot(x_wrf_ncar_c_iv,    pdf_wrf_ncar_c_iv,    linewidth=1, color='green',   label='WRF-NCAR')
plt.plot(x_wrf_ucan_c_iv,    pdf_wrf_ucan_c_iv,    linewidth=1, color='orange',  label='WRF-UCAN')
plt.text(text_, text_8, 'INMET = {0} ({1})'.format(round(perc_inmet_smn_c_iv, 1),      round(max_inmet_smn_c_iv, 1)),   color='black',   fontsize=font_size)
plt.text(text_, text_7, 'ERA5 = {0} ({1})'.format(round(perc_era5_c_iv, 1),	       round(max_era5_c_iv, 1)),	color='red',	 fontsize=font_size)
plt.text(text_, text_6, 'Reg4 = {0} ({1})'.format(round(perc_reg_usp_c_iv, 1),         round(max_reg_usp_c_iv, 1)),     color='blue',    fontsize=font_size)
plt.text(text_, text_5, 'Reg5-Holt3 = {0} ({1})'.format(round(perc_reg_ictp_c_iv, 1),  round(max_reg_ictp_c_iv, 1)),    color='magenta', fontsize=font_size)
plt.text(text_, text_4, 'Reg5-Holt = {0} ({1})'.format(round(perc_reg_ictp_i_c_iv, 1), round(max_reg_ictp_i_c_iv, 1)),  color='gray',	fontsize=font_size)
plt.text(text_, text_3, 'Reg5-UW = {0} ({1})'.format(round(perc_reg_ictp_ii_c_iv, 1),  round(max_reg_ictp_ii_c_iv, 1)), color='brown',   fontsize=font_size)
plt.text(text_, text_2, 'WRF-NCAR = {0} ({1})'.format(round(perc_wrf_ncar_c_iv, 1),    round(max_wrf_ncar_c_iv, 1)),    color='green',   fontsize=font_size)
plt.text(text_, text_1, 'WRF-UCAN = {0} ({1})'.format(round(perc_wrf_ucan_c_iv, 1),    round(max_wrf_ucan_c_iv, 1)),    color='orange',  fontsize=font_size)
plt.title('(d) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.xlim(xvmin, xvmax)
plt.ylim(yvmin, yvmax)
plt.xticks(np.arange(xvmin, xvmax_, xint_), fontsize=font_size)
plt.yticks(np.arange(yvmin, yvmax_, yint_), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax5 = fig.add_subplot(6, 2, 5)
plt.plot(x_inmet_smn_c_v,   pdf_inmet_smn_c_v,   linewidth=2, color='black',   label='INMET+SMN')
plt.plot(x_era5_c_v,        pdf_era5_c_v,        linewidth=1, color='red',     label='ERA5')
plt.plot(x_reg_usp_c_v,     pdf_reg_usp_c_v,     linewidth=1, color='blue',    label='Reg4')
plt.plot(x_reg_ictp_c_v,    pdf_reg_ictp_c_v,    linewidth=1, color='magenta', label='Reg5-Holt3')
plt.plot(x_reg_ictp_i_c_v,  pdf_reg_ictp_i_c_v,  linewidth=1, color='gray',    label='Reg5-Holt')
plt.plot(x_reg_ictp_ii_c_v, pdf_reg_ictp_ii_c_v, linewidth=1, color='brown',   label='Reg5-UW')
plt.plot(x_wrf_ncar_c_v,    pdf_wrf_ncar_c_v,    linewidth=1, color='green',   label='WRF-NCAR')
plt.plot(x_wrf_ucan_c_v,    pdf_wrf_ucan_c_v,    linewidth=1, color='orange',  label='WRF-UCAN')
plt.text(text_, text_8, 'INMET = {0} ({1})'.format(round(perc_inmet_smn_c_v, 1),      round(max_inmet_smn_c_v, 1)),   color='black',   fontsize=font_size)
plt.text(text_, text_7, 'ERA5 = {0} ({1})'.format(round(perc_era5_c_v, 1),            round(max_era5_c_v, 1)),        color='red',     fontsize=font_size)
plt.text(text_, text_6, 'Reg4 = {0} ({1})'.format(round(perc_reg_usp_c_v, 1),         round(max_reg_usp_c_v, 1)),     color='blue',    fontsize=font_size)
plt.text(text_, text_5, 'Reg5-Holt3 = {0} ({1})'.format(round(perc_reg_ictp_c_v, 1),  round(max_reg_ictp_c_v, 1)),    color='magenta', fontsize=font_size)
plt.text(text_, text_4, 'Reg5-Holt = {0} ({1})'.format(round(perc_reg_ictp_i_c_v, 1), round(max_reg_ictp_i_c_v, 1)),  color='gray',    fontsize=font_size)
plt.text(text_, text_3, 'Reg5-UW = {0} ({1})'.format(round(perc_reg_ictp_ii_c_v, 1),  round(max_reg_ictp_ii_c_v, 1)), color='brown',   fontsize=font_size)
plt.text(text_, text_2, 'WRF-NCAR = {0} ({1})'.format(round(perc_wrf_ncar_c_v, 1),    round(max_wrf_ncar_c_v, 1)),    color='green',   fontsize=font_size)
plt.text(text_, text_1, 'WRF-UCAN = {0} ({1})'.format(round(perc_wrf_ucan_c_v, 1),    round(max_wrf_ucan_c_v, 1)),    color='orange',  fontsize=font_size)
plt.title('(e) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.xlim(xvmin, xvmax)
plt.ylim(yvmin, yvmax)
plt.xticks(np.arange(xvmin, xvmax_, xint_), fontsize=font_size)
plt.yticks(np.arange(yvmin, yvmax_, yint_), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_pdf_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
