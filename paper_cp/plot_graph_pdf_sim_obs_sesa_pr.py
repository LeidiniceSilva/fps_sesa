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

var = 'pr'
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
	
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii  = [], [], [], [], [], [], [], []
	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue
		if inmet[i][3] <= -48 and inmet[i][2] <= -16.5:
			yy=inmet[i][2]
			xx=inmet[i][3]

			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.values
			mean_vi.append(d_vi)
		
	return  mean_vi


def import_smn_i():

	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii  = [], [], [], [], [], [], [], []
	for i in range(1, 73):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
			
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.values
		mean_vi.append(d_vi)

	return  mean_vi


def compute_pdf(pr_hourly, nbins=200):

	p999 = np.nanpercentile(pr_hourly, 99.9)

	pr = np.asarray(pr_hourly)
	pr = pr[pr > 0]  

	kde = gaussian_kde(pr)
	x = np.linspace(pr.min(), pr.max(), nbins)
	pdf = kde(x)
	

	return x, pdf, p999
    
	
# Import dataset
clim_vi_x = import_inmet()			
clim_vi_y = import_smn_i()

inmet_smn = clim_vi_x + clim_vi_y

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
# Average
inmet_smn_c_i   = np.nanmean(inmet_smn_i, axis=0)
inmet_smn_c_ii   = np.nanmean(inmet_smn_ii, axis=0)
inmet_smn_c_iii   = np.nanmean(inmet_smn_iii, axis=0)
inmet_smn_c_iv   = np.nanmean(inmet_smn_iv, axis=0)
inmet_smn_c_v   = np.nanmean(inmet_smn_v, axis=0)

print(inmet_smn_c_i)
print(len(inmet_smn_c_i))
print()
print(inmet_smn_c_ii)
print(len(inmet_smn_c_ii))
print()
print(inmet_smn_c_iii)
print(len(inmet_smn_c_iii))
print()
print(inmet_smn_c_iv)
print(len(inmet_smn_c_iv))
print()
print(inmet_smn_c_v)
print(len(inmet_smn_c_v))

x_inmet_smn_c_i, pdf_inmet_smn_c_i, p999_inmet_smn_c_i = compute_pdf(inmet_smn_c_i)
x_inmet_smn_c_ii, pdf_inmet_smn_c_ii, p999_inmet_smn_c_ii = compute_pdf(inmet_smn_c_ii)
x_inmet_smn_c_iii, pdf_inmet_smn_c_iii, p999_inmet_smn_c_iii = compute_pdf(inmet_smn_c_iii)
x_inmet_smn_c_iv, pdf_inmet_smn_c_iv, p999_inmet_smn_c_iv = compute_pdf(inmet_smn_c_iv)
x_inmet_smn_c_v, pdf_inmet_smn_c_v, p999_inmet_smn_c_v = compute_pdf(inmet_smn_c_v)

# Plot figure
fig = plt.figure(figsize=(8, 18))
font_size = 8

ax1 = fig.add_subplot(6, 2, 1)
plt.plot(x_inmet_smn_c_i, pdf_inmet_smn_c_i, linewidth=1, color='black', label='INMET+SMN')
plt.axvline(p999_inmet_smn_c_i, linestyle='--', linewidth=0.75, color='black')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(0, 4)
plt.xticks(np.arange(0, 6.5, 0.5), fontsize=font_size)
plt.yticks(np.arange(0, 4.5, 0.5), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')
ax1.legend(loc=2, ncol=3, fontsize=8, frameon=False)

ax2 = fig.add_subplot(6, 2, 2)
ax2.plot(x_inmet_smn_c_ii, pdf_inmet_smn_c_ii, linewidth=1, color='black', label='INMET+SMN')
plt.axvline(p999_inmet_smn_c_ii, linestyle='--', linewidth=0.75, color='black')
plt.title('(b) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(0, 4)
plt.xticks(np.arange(0, 6.5, 0.5), fontsize=font_size)
plt.yticks(np.arange(0, 4.5, 0.5), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax3 = fig.add_subplot(6, 2, 3)
plt.plot(x_inmet_smn_c_iii, pdf_inmet_smn_c_iii, linewidth=1, color='black', label='INMET+SMN')
plt.axvline(p999_inmet_smn_c_iii, linestyle='--', linewidth=0.75, color='black')
plt.title('(c) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(0, 4)
plt.xticks(np.arange(0, 6.5, 0.5), fontsize=font_size)
plt.yticks(np.arange(0, 4.5, 0.5), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax4 = fig.add_subplot(6, 2, 4)
plt.plot(x_inmet_smn_c_iv, pdf_inmet_smn_c_iv, linewidth=1, color='black', label='INMET+SMN')
plt.axvline(p999_inmet_smn_c_iv, linestyle='--', linewidth=0.75, color='black')
plt.title('(d) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(0, 4)
plt.xticks(np.arange(0, 6.5, 0.5), fontsize=font_size)
plt.yticks(np.arange(0, 4.5, 0.5), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

ax5 = fig.add_subplot(6, 2, 5)
plt.plot(x_inmet_smn_c_v, pdf_inmet_smn_c_v, linewidth=1, color='black', label='INMET+SMN')
plt.axvline(p999_inmet_smn_c_v, linestyle='--', linewidth=0.75, color='black')
plt.title('(e) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=font_size, fontweight='bold')
plt.xlabel('Precipitation (mm h⁻¹)', fontsize=font_size, fontweight='bold')
plt.xlim(0, 6)
plt.ylim(0, 4)
plt.xticks(np.arange(0, 6.5, 0.5), fontsize=font_size)
plt.yticks(np.arange(0, 4.5, 0.5), fontsize=font_size)
plt.grid(True, alpha=0.5, linestyle='--')

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_pdf_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
