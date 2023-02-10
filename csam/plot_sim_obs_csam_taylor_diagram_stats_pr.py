# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot annual cycle statistics in the csam  domain"

import os
import numpy as np
import xarray as xr
import scipy.stats as st
import matplotlib.pyplot as plt

from dict_csam_inmet_stations import inmet


def import_inmet(dt):
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []

	# Select lat and lon 
	for i in range(1, 155):
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		# reading regcm 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2019-01-01','2021-12-31'))
		d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		mean_i.append(values_i*86400)
		# Reading inmet 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2019-01-01','2021-12-31'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		mean_ii.append(values_ii*24)				
		# reading cmorph 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_4km_mon_20180101-20211231.nc')
		d_iii = d_iii.cmorph.sel(time=slice('2019-01-01','2021-12-31'))
		d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		mean_iii.append(values_iii)
		# reading era5 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_iv = d_iv.tp.sel(time=slice('2019-01-01','2021-12-31'))
		d_iv = d_iv.groupby('time.month').mean('time')
		d_iv = d_iv.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		values_iv = d_iv.values
		mean_iv.append(values_iv)
		
	return mean_i, mean_ii, mean_iii, mean_iv
	

dt = 'H_2018-01-01_2021-12-31'

print('Import dataset')
# Import dataset
mean_i, mean_ii, mean_iii, mean_iv = import_inmet(dt)			

list_hc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 0,
3, 1, 3, 1, 1, 1, 3, 1, 3, 4, 3, 1, 3, 1, 3, 3, 3, 3, 1, 3, 0, 3, 3, 3, 3,
3, 4, 3, 4, 0, 0, 0, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 3, 4,
3, 0, 4, 0, 4, 4, 4, 4, 4, 0, 4, 0, 0, 4, 0, 4, 4, 4, 4, 4, 4, 0, 4, 0, 2,
4, 4, 2, 4, 4, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 0, 0, 0, 2, 0, 2, 2, 2, 
0, 2, 2, 2, 2, 2, 4, 0]

clim_regcm = mean_i
clim_inmet = mean_ii
clim_cmorph = mean_iii
clim_era5 = mean_iv

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

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

regcm_i = []
regcm_ii = []
regcm_iii = []
regcm_iv = []
regcm_v = []

inmet_i = []
inmet_ii = []
inmet_iii = []
inmet_iv = []
inmet_v = []

cmorph_i = []
cmorph_ii = []
cmorph_iii = []
cmorph_iv = []
cmorph_v = []

era5_i = []
era5_ii = []
era5_iii = []
era5_iv = []
era5_v = []

for c_i in count_i:
	regcm_i.append(clim_regcm[c_i])
	inmet_i.append(clim_inmet[c_i])
	cmorph_i.append(clim_cmorph[c_i])
	era5_i.append(clim_era5[c_i])

for c_ii in count_ii:
	regcm_ii.append(clim_regcm[c_ii])
	inmet_ii.append(clim_inmet[c_ii])
	cmorph_ii.append(clim_cmorph[c_ii])
	era5_ii.append(clim_era5[c_ii])
	
for c_iii in count_iii:
	regcm_iii.append(clim_regcm[c_iii])
	inmet_iii.append(clim_inmet[c_iii])
	cmorph_iii.append(clim_cmorph[c_iii])
	era5_iii.append(clim_era5[c_iii])
	
for c_iv in count_iv:
	regcm_iv.append(clim_regcm[c_iv])
	inmet_iv.append(clim_inmet[c_iv])
	cmorph_iv.append(clim_cmorph[c_iv])
	era5_iv.append(clim_era5[c_iv])
	
for c_v in count_v:
	regcm_v.append(clim_regcm[c_v])
	inmet_v.append(clim_inmet[c_v])
	cmorph_v.append(clim_cmorph[c_v])
	era5_v.append(clim_era5[c_v])

regcm_cluster_i = np.nanmean(regcm_i, axis=0)
inmet_cluster_i = np.nanmean(inmet_i, axis=0)
cmorph_cluster_i = np.nanmean(cmorph_i, axis=0)
era5_cluster_i = np.nanmean(era5_i, axis=0)

regcm_cluster_ii = np.nanmean(regcm_ii, axis=0)
inmet_cluster_ii = np.nanmean(inmet_ii, axis=0)
cmorph_cluster_ii = np.nanmean(cmorph_ii, axis=0)
era5_cluster_ii = np.nanmean(era5_ii, axis=0)

regcm_cluster_iii = np.nanmean(regcm_iii, axis=0)
inmet_cluster_iii = np.nanmean(inmet_iii, axis=0)
cmorph_cluster_iii = np.nanmean(cmorph_iii, axis=0)
era5_cluster_iii = np.nanmean(era5_iii, axis=0)

regcm_cluster_iv = np.nanmean(regcm_iv, axis=0)
inmet_cluster_iv = np.nanmean(inmet_iv, axis=0)
cmorph_cluster_iv = np.nanmean(cmorph_iv, axis=0)
era5_cluster_iv = np.nanmean(era5_iv, axis=0)

regcm_cluster_v = np.nanmean(regcm_v, axis=0)
inmet_cluster_v = np.nanmean(inmet_v, axis=0)
cmorph_cluster_v = np.nanmean(cmorph_v, axis=0)
era5_cluster_v = np.nanmean(era5_v, axis=0)

# Compute standard deviation
std_regcm_cluster_i = regcm_cluster_i.std(ddof=0)
std_inmet_cluster_i = inmet_cluster_i.std(ddof=0)
std_cmorph_cluster_i = cmorph_cluster_i.std(ddof=0)
std_era5_cluster_i = era5_cluster_i.std(ddof=0)

std_regcm_cluster_ii = regcm_cluster_ii.std(ddof=0)
std_inmet_cluster_ii = inmet_cluster_ii.std(ddof=0)
std_cmorph_cluster_ii = cmorph_cluster_ii.std(ddof=0)
std_era5_cluster_ii = era5_cluster_ii.std(ddof=0)

std_regcm_cluster_iii = regcm_cluster_iii.std(ddof=0)
std_inmet_cluster_iii = inmet_cluster_iii.std(ddof=0)
std_cmorph_cluster_iii = cmorph_cluster_iii.std(ddof=0)
std_era5_cluster_iii = era5_cluster_iii.std(ddof=0)

std_regcm_cluster_iv = regcm_cluster_iv.std(ddof=0)
std_inmet_cluster_iv = inmet_cluster_iv.std(ddof=0)
std_cmorph_cluster_iv = cmorph_cluster_iv.std(ddof=0)
std_era5_cluster_iv = era5_cluster_iv.std(ddof=0)

std_regcm_cluster_v = regcm_cluster_v.std(ddof=0)
std_inmet_cluster_v = inmet_cluster_v.std(ddof=0)
std_cmorph_cluster_v = cmorph_cluster_v.std(ddof=0)
std_era5_cluster_v = era5_cluster_v.std(ddof=0)

# Compute correlation
corr_regcm_inmet_cluster_i = st.pearsonr(inmet_cluster_i, regcm_cluster_i)[0]
corr_regcm_cmorph_cluster_i = st.pearsonr(cmorph_cluster_i, regcm_cluster_i)[0]
corr_regcm_era5_cluster_i = st.pearsonr(era5_cluster_i, regcm_cluster_i)[0]

corr_regcm_inmet_cluster_ii = st.pearsonr(inmet_cluster_ii, regcm_cluster_ii)[0]
corr_regcm_cmorph_cluster_ii = st.pearsonr(cmorph_cluster_ii, regcm_cluster_ii)[0]
corr_regcm_era5_cluster_ii = st.pearsonr(era5_cluster_ii, regcm_cluster_ii)[0]

corr_regcm_inmet_cluster_iii = st.pearsonr(inmet_cluster_iii, regcm_cluster_iii)[0]
corr_regcm_cmorph_cluster_iii = st.pearsonr(cmorph_cluster_iii, regcm_cluster_iii)[0]
corr_regcm_era5_cluster_iii = st.pearsonr(era5_cluster_iii, regcm_cluster_iii)[0]

corr_regcm_inmet_cluster_iv = st.pearsonr(inmet_cluster_iv, regcm_cluster_iv)[0]
corr_regcm_cmorph_cluster_iv = st.pearsonr(cmorph_cluster_iv, regcm_cluster_iv)[0]
corr_regcm_era5_cluster_iv = st.pearsonr(era5_cluster_iv, regcm_cluster_iv)[0]

corr_regcm_inmet_cluster_v = st.pearsonr(inmet_cluster_v, regcm_cluster_v)[0]
corr_regcm_cmorph_cluster_v = st.pearsonr(cmorph_cluster_v, regcm_cluster_v)[0]
corr_regcm_era5_cluster_v = st.pearsonr(era5_cluster_v, regcm_cluster_v)[0]

print(std_regcm_cluster_i, std_inmet_cluster_i, std_cmorph_cluster_i, std_era5_cluster_i)
print(std_regcm_cluster_ii, std_inmet_cluster_ii, std_cmorph_cluster_ii, std_era5_cluster_ii)
print(std_regcm_cluster_iii, std_inmet_cluster_iii, std_cmorph_cluster_iii, std_era5_cluster_iii)
print(std_regcm_cluster_iv, std_inmet_cluster_iv, std_cmorph_cluster_iv, std_era5_cluster_iv)
print(std_regcm_cluster_v, std_inmet_cluster_v, std_cmorph_cluster_v, std_era5_cluster_v)

print(corr_regcm_inmet_cluster_i, corr_regcm_cmorph_cluster_i, corr_regcm_era5_cluster_i)
print(corr_regcm_inmet_cluster_ii, corr_regcm_cmorph_cluster_ii, corr_regcm_era5_cluster_ii)
print(corr_regcm_inmet_cluster_iii, corr_regcm_cmorph_cluster_iii, corr_regcm_era5_cluster_iii)
print(corr_regcm_inmet_cluster_iv, corr_regcm_cmorph_cluster_iv, corr_regcm_era5_cluster_iv)
print(corr_regcm_inmet_cluster_v, corr_regcm_cmorph_cluster_v, corr_regcm_era5_cluster_v)
exit()

