# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script check nan in automatic weather station"

import os
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

institution = 'smn_ii'

if institution == 'inmet':
	for i in range(1, 567):
		skip_list_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
		if i in skip_list_i:
			continue
		
		# Reading inmet 
		d_i = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/inmet/inmet_br/inmet_nc/daily/pre/' + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.values
		x_i = sum(math.isnan(x) for x in d_i)
		n_i = len(d_i)
		percent_missing_i = x_i * 100 / n_i
		if percent_missing_i > 30:
			print(
				f'Reading weather station: {i} {inmet[i][0]} {inmet[i][1]} | '
				f'Missing: {percent_missing_i:.1f}% ({x_i}/{n_i})'
			)
	
elif institution == 'smn_i':
	for j in range(1, 72):

		# Reading Uruguai weather stations
		d_j = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_i/smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[j][0]))
		d_j = d_j.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_j = d_j.values	
		y_j = sum(math.isnan(y) for y in d_j)
		n_j = len(d_j)
		percent_missing_j = y_j * 100 / n_j
		if percent_missing_j > 30:
			print(
				f'Reading weather station: {j} {smn_i[j][0]} {smn_i[j][1]} | '
				f'Missing: {percent_missing_j:.1f}% ({y_j}/{n_j})'
			)

else:
	for k in range(1, 110):

		# Reading Argentina weather stations
		d_k = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_ii/smn_nc/pre/' + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[k][0]))
		d_k = d_k.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_k = d_k.values
		z_k = sum(math.isnan(z) for z in d_k)
		n_k = len(d_k)
		percent_missing_z = z_k * 100 / n_k
		if percent_missing_z > 30:
			print(
				f'Reading weather station: {k} {smn_ii[k][0]} {smn_ii[k][1]} | '
				f'Missing: {percent_missing_z:.1f}% ({z_k}/{n_k})'
			)


