# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script check nan in automatic weather station"

import math
import xarray as xr

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

institution = 'smn_ii'
stations_30pcent_nan = []


def comp_nan(value):

	x_ = sum(math.isnan(x) for x in value)
	n_ = len(value)
	percent_missing = x_ * 100 / n_
	
	return percent_missing
	
		
if institution == 'inmet':
	for i in range(1, 567):
		skip_list_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
		if i in skip_list_i:
			continue
		d_i = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/inmet/inmet_br/inmet_nc/daily/pre/' + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.values
		pcent_i = comp_nan(d_i)
		if pcent_i > 30:
			stations_30pcent_nan.append(i)
	    
elif institution == 'smn_i':
	for j in range(1, 73):
		d_j = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_i/smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[j][0]))
		d_j = d_j.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_j = d_j.values	
		pcent_j = comp_nan(d_j)
		if pcent_j > 30:
			stations_30pcent_nan.append(j)

else:
	for k in range(1, 110):
		d_k = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_ii/smn_nc/pre/' + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[k][0]))
		d_k = d_k.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_k = d_k.values
		pcent_k = comp_nan(d_k)
		if pcent_k > 30:
			stations_30pcent_nan.append(k)

# Result
print("Stations with >30% NaNs:")
print(stations_30pcent_nan)




