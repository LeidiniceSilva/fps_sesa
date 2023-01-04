# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script compte table's stats of annual cycle to weather station"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import texttable as tt
import scipy.stats as st

from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn


def import_inmet_era5(dt):
	
	name_i = []
	ix = []		  
	iy = []
	corr_ii = []
	pvalue_ii = []
	rmse_ii = []

	# Select lat and lon 
	for i in range(1, 289):
		
		if i == 4:
			continue
		if i == 12:
			continue
		if i == 45:
			continue
		if i == 55:
			continue
		if i == 77:
			continue
		if i == 85:
			continue
		if i == 90:
			continue
		if i == 98:
			continue
		if i == 99:
			continue
		if i == 100:
			continue
		if i == 118:
			continue
		if i == 122:
			continue
		if i == 130:
			continue
		if i == 135:
			continue
		if i == 151:
			continue
		if i == 155:
			continue
		if i == 159:
			continue
		if i == 160:
			continue
		if i == 163:
			continue
		if i == 164:
			continue
		if i == 181:
			continue
		if i == 183:
			continue
		if i == 186:
			continue
		if i == 187:
			continue
		if i == 188:
			continue
		if i == 209:
			continue
		if i == 216:
			continue
		if i == 228:
			continue
		if i == 236:
			continue
		if i == 245:
			continue
		if i == 246:
			continue
		if i == 262:
			continue
		if i == 268:
			continue
		if i == 273:
			continue
		if i == 287:
			continue
				
		name_i.append(inmet[i][1])
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])

		print('Reading INMET weather station:', i, inmet[i][0], inmet[i][1])
		# Reading inmet weather station	
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_{1}.nc'.format(inmet[i][0], dt))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_i = d_i.groupby('time.hour').mean('time')
		values_i = d_i.values
		clim_i = values_i*24

		# reading era5 reanalisis
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
		d_ii = d_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_ii = d_ii.sel(latitude=inmet[i][2], longitude=inmet[i][3], method='nearest')
		d_ii = d_ii.groupby('time.hour').mean('time')
		values_ii = d_ii.values
		clim_ii = values_ii*24
			
		corr_i, pvalue_i = st.pearsonr(clim_i, clim_ii)
		corr_ii.append(corr_i)
		pvalue_ii.append(pvalue_i)

		rmse_i = np.sqrt(((np.array(clim_ii) - np.array(clim_i)) ** 2).mean()) 
		rmse_ii.append(rmse_i)

	return name_i, iy, ix, corr_ii, pvalue_ii, rmse_ii


def import_urug_smn_era5(dt):

	name_j = []
	jy = []
	jx = []
	corr_jj = []
	pvalue_jj = []
	rmse_jj = []
	
	# Select lat and lon 
	for j in range(1, 72):
		
		name_j.append(urug_smn[j][0])
		jy.append(urug_smn[j][1])
		jx.append(urug_smn[j][2])		

		print('Reading Uruguai weather station:', j, urug_smn[j][0])	
		# Reading Uruguai weather stations
		d_j = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
		d_j = d_j.pre.sel(time=slice('2018-01-01','2021-12-31'))
		d_j = d_j.groupby('time.hour').mean('time')
		values_j = d_j.values
		clim_j = values_j*24
		
		# Reading ERA5 reanalisis
		d_jj = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
		d_jj = d_jj.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_jj = d_jj.sel(latitude=urug_smn[j][1], longitude=urug_smn[j][2], method='nearest')
		d_jj = d_jj.groupby('time.hour').mean('time')
		values_jj = d_jj.values
		clim_jj = values_jj*24
		
		corr_j, pvalue_j = st.pearsonr(clim_j, clim_jj)
		corr_jj.append(corr_j)
		pvalue_jj.append(pvalue_j)

		rmse_j = np.sqrt(((np.array(clim_jj) - np.array(clim_j)) ** 2).mean()) 
		rmse_jj.append(rmse_j)
				
	return name_j, jy, jx, corr_jj, pvalue_jj, rmse_jj


def import_arg_emas_era5(dt):

	name_k = []
	ky = []
	kx = []
	corr_kk = []
	pvalue_kk = []
	rmse_kk = []
	
	# Select lat and lon 
	for k in range(1, 88):

		if k == 9:
			continue
		if k == 17:
			continue
		if k == 24:
			continue
		if k == 26:
			continue
		if k == 43:
			continue
		if k == 72:
			continue
		if k == 76:
			continue
			
		name_k.append(arg_emas[k][0])
		ky.append(arg_emas[k][2])
		kx.append(arg_emas[k][1])		
		
		print('Reading Argentina weather station:', k, arg_emas[k][0])	
		# Reading Argentina weather stations
		d_k = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[k][0], dt))
		d_k = d_k.precip.sel(time=slice('2018-01-01','2021-12-31'))
		d_k = d_k.groupby('time.hour').mean('time')
		values_k = d_k.values
		clim_k = values_k*24
		
		# Reading ERA5 reanalisis
		d_kk = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
		d_kk = d_kk.tp.sel(time=slice('2018-01-01','2021-12-31'))
		d_kk = d_kk.sel(latitude=arg_emas[k][2], longitude=arg_emas[k][1], method='nearest')
		d_kk = d_kk.groupby('time.hour').mean('time')
		values_kk = d_kk.values
		clim_kk = values_kk*24
		
		corr_k, pvalue_k = st.pearsonr(clim_k, clim_kk)
		corr_kk.append(corr_k)
		pvalue_kk.append(pvalue_k)

		rmse_k = np.sqrt(((np.array(clim_kk) - np.array(clim_k)) ** 2).mean()) 
		rmse_kk.append(rmse_k)
				
	return name_k, ky, kx, corr_kk, pvalue_kk, rmse_kk
	
	
dt = 'H_2018-01-01_2021-12-31'
tab = tt.Texttable(max_width=200)
tab_inform = [[]]

# Import name, latitude, longitude, correlation, pvalue and rmse
name_i, iy, ix, corr_ii, pvalue_ii, rmse_ii = import_inmet_era5(dt)			
name_j, jy, jx, corr_jj, pvalue_jj, rmse_jj = import_urug_smn_era5(dt)
name_k, ky, kx, corr_kk, pvalue_kk, rmse_kk = import_arg_emas_era5(dt)

name_tot = name_i+name_j+name_k
lat_yy = iy+jy+ky
lon_xx = ix+jx+kx
corr_tot = corr_ii+corr_jj+corr_kk
pvalue_tot = pvalue_ii+pvalue_jj+pvalue_kk
rmse_tot = rmse_ii+rmse_jj+rmse_kk

for idx in range(0, 403):
	
	print(name_tot[idx], lat_yy[idx], lon_xx[idx], corr_tot[idx], pvalue_tot[idx], rmse_tot[idx])
	
	name = name_tot[idx]
	latitude = lat_yy[idx]
	longitude = lon_xx[idx]
	correlation = corr_tot[idx]
	pvalue = pvalue_tot[idx]
	rmse = rmse_tot[idx]
	
	tab_inform.append([name, latitude, longitude, correlation, pvalue, rmse])
	
tab.add_rows(tab_inform)
tab.set_cols_align(['c', 'c', 'c', 'c', 'c', 'c'])
tab.header([u'Weather station', u'Latitude', u'Longitude', u'CORR', u'P_VALUE', U'RMSE'])
table = str(tab.draw())
	
file_name = 'table_stats_diurnal_cycle_stations_2018-2021.asc'
file_save = open(file_name, 'w')
file_save.write(table)
file_save.close()        
exit()

