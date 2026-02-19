# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot PDFs"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from matplotlib.patches import Polygon
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

var = 't2m'
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
	
	mean = []
	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue
		if inmet[i][3] <= -48 and inmet[i][2] <= -16.5:
			yy=inmet[i][2]
			xx=inmet[i][3]

			if var == 't2m':
				# Reading inmet 
				df_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/tmp/'.format(path) + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
				df_i = df_i.tmp.sel(time=slice('2018-06-01','2021-05-31'))
				df_i = df_i.values
				mean.append(df_i)
			else:
				# Reading inmet 
				df_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/uv/'.format(path) + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
				df_i = df_i.uv.sel(time=slice('2018-06-01','2021-05-31'))
				df_i = df_i.values
				mean.append(df_i)
					
	return  mean



def compute_stats(dataset):

	if var == 't2m':
		perc = 97.5
	else:
		perc = 95
		
	var_ = np.asarray(dataset)
	horas = np.arange(len(var_)) % 24
	
	# Percentile (threshold)
	perc_global = np.nanpercentile(var_, perc)
	
	perc_hour = np.zeros(24)
	freq = np.zeros(24)
	intens = np.zeros(24)
	
	for h in range(24):
		var_h = var_[horas == h]
		valid = ~np.isnan(var_h)
	
		# Percentil diurno
		perc_hour[h] = np.nanpercentile(var_h, perc) if np.sum(valid) > 0 else np.nan
	
		# Extreme events 
		extreme = var_h[valid][var_h[valid] >= perc_global]

		# Frequency (%) 
		if np.sum(valid) > 0:
			freq[h] = 100.0 * len(extreme) / np.sum(valid)
		else:
			freq[h] = np.nan

		# Intensity 
		intens[h] = np.nanmean(extreme) if len(extreme) > 0 else np.nan

	return perc_hour, freq, intens
	
	
# Import dataset
inmet_smn = import_inmet()			

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

for c_i in count_i:
	inmet_smn_i.append(inmet_smn[c_i])

for c_ii in count_ii:
	inmet_smn_ii.append(inmet_smn[c_ii])
	
for c_iii in count_iii:
	inmet_smn_iii.append(inmet_smn[c_iii])
	
for c_iv in count_iv:
	inmet_smn_iv.append(inmet_smn[c_iv])
	
for c_v in count_v:
	inmet_smn_v.append(inmet_smn[c_v])

# Group I
inmet_smn_c_i = np.concatenate(inmet_smn_i)
inmet_smn_c_ii = np.concatenate(inmet_smn_ii)
inmet_smn_c_iii = np.concatenate(inmet_smn_iii)
inmet_smn_c_iv = np.concatenate(inmet_smn_iv)
inmet_smn_c_v = np.concatenate(inmet_smn_v)

perc_inmet_smn_c_i, freq_inmet_smn_c_i, int_inmet_smn_c_i = compute_stats(inmet_smn_c_i)
perc_inmet_smn_c_ii, freq_inmet_smn_c_ii, int_inmet_smn_c_ii = compute_stats(inmet_smn_c_ii)
perc_inmet_smn_c_iii, freq_inmet_smn_c_iii, int_inmet_smn_c_iii = compute_stats(inmet_smn_c_iii)
perc_inmet_smn_c_iv, freq_inmet_smn_c_iv, int_inmet_smn_c_iv = compute_stats(inmet_smn_c_iv)
perc_inmet_smn_c_v, freq_inmet_smn_c_v, int_inmet_smn_c_v = compute_stats(inmet_smn_c_v)

print(perc_inmet_smn_c_i)
print(freq_inmet_smn_c_i)
print(int_inmet_smn_c_i)

# Plot figure
fig = plt.figure(figsize=(14, 16))
time = np.arange(0, 24)
font_size = 8

pvmin = 20
pvmax = 40
pvmax_ = 42
pint_ = 2

fvmin = 0
fvmax = 2
fvmax_ = 2.2
fint_ = 0.2

ivmin = 20
ivmax = 40
ivmax_ = 42
iint_ = 2

ax = fig.add_subplot(5, 3, 1)
plt.plot(time, perc_inmet_smn_c_i, linewidth=1, color='black', label='INMET+SMN')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('P99 (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(pvmin, pvmax)
plt.yticks(np.arange(pvmin, pvmax_, pint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 2)
plt.plot(time, freq_inmet_smn_c_i, linewidth=1, color='black', label='INMET+SMN')
plt.title('(b) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (%)', fontsize=font_size, fontweight='bold')
plt.ylim(fvmin, fvmax)
plt.yticks(np.arange(fvmin, fvmax_, fint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 3)
plt.plot(time, int_inmet_smn_c_i, linewidth=1, color='black', label='INMET+SMN')
plt.title('(c) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Intensity (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(ivmin, ivmax)
plt.yticks(np.arange(ivmin, ivmax_, iint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 4)
plt.plot(time, perc_inmet_smn_c_ii, linewidth=1, color='black', label='INMET+SMN')
plt.title('(d) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('P99 (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(pvmin, pvmax)
plt.yticks(np.arange(pvmin, pvmax_, pint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 5)
plt.plot(time, freq_inmet_smn_c_ii, linewidth=1, color='black', label='INMET+SMN')
plt.ylabel('Frequency (%)', fontsize=font_size, fontweight='bold')
plt.title('(e) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(fvmin, fvmax)
plt.yticks(np.arange(fvmin, fvmax_, fint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 6)
plt.plot(time, int_inmet_smn_c_ii, linewidth=1, color='black', label='INMET+SMN')
plt.title('(f) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Intensity (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(ivmin, ivmax)
plt.yticks(np.arange(ivmin, ivmax_, iint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 7)
plt.plot(time, perc_inmet_smn_c_iii, linewidth=1, color='black', label='INMET+SMN')
plt.ylabel('P99 (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.title('(g) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(pvmin, pvmax)
plt.yticks(np.arange(pvmin, pvmax_, pint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 8)
plt.plot(time, freq_inmet_smn_c_iii, linewidth=1, color='black', label='INMET+SMN')
plt.title('(h) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (%)', fontsize=font_size, fontweight='bold')
plt.ylim(fvmin, fvmax)
plt.yticks(np.arange(fvmin, fvmax_, fint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 9)
plt.plot(time, int_inmet_smn_c_iii, linewidth=1, color='black', label='INMET+SMN')
plt.title('(i) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Intensity (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(ivmin, ivmax)
plt.yticks(np.arange(ivmin, ivmax_, iint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 10)
plt.plot(time, perc_inmet_smn_c_iv, linewidth=1, color='black', label='INMET+SMN')
plt.ylabel('P99 (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.title('(j) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(pvmin, pvmax)
plt.yticks(np.arange(pvmin, pvmax_, pint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 11)
plt.plot(time, freq_inmet_smn_c_iv, linewidth=1, color='black', label='INMET+SMN')
plt.title('(k) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (%)', fontsize=font_size, fontweight='bold')
plt.ylim(fvmin, fvmax)
plt.yticks(np.arange(fvmin, fvmax_, fint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 12)
plt.plot(time, int_inmet_smn_c_iv, linewidth=1, color='black', label='INMET+SMN')
plt.title('(l) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Intensity (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(ivmin, ivmax)
plt.yticks(np.arange(ivmin, ivmax_, iint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 13)
plt.plot(time, perc_inmet_smn_c_v, linewidth=1, color='black', label='INMET+SMN')
plt.ylabel('P99 (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.title('(m) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylim(pvmin, pvmax)
plt.yticks(np.arange(pvmin, pvmax_, pint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 14)
plt.plot(time, freq_inmet_smn_c_v, linewidth=1, color='black', label='INMET+SMN')
plt.title('(n) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (%)', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylim(fvmin, fvmax)
plt.yticks(np.arange(fvmin, fvmax_, fint_))
plt.grid(True, alpha=0.5, linestyle='--')

ax = fig.add_subplot(5, 3, 15)
plt.plot(time, int_inmet_smn_c_v, linewidth=1, color='black', label='INMET+SMN')
plt.title('(o) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Intensity (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.xlabel('Hours', fontsize=font_size, fontweight='bold')
plt.ylim(ivmin, ivmax)
plt.yticks(np.arange(ivmin, ivmax_, iint_))
plt.grid(True, alpha=0.5, linestyle='--')

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_diurnal_cycle_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
