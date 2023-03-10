# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot boxplt of the INMET weather station"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_sesa_inmet_stations import inmet
from pylab import setp


def import_inmet(dt):
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []

	# Select lat and lon 
	for i in range(1, 155):
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])

		# reading regcm 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2019-01-01','2021-12-31'))
		d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_i = d_i.values
		mean_i.append(d_i*86400)

		# Reading inmet 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_ii = d_ii.pre.sel(time=slice('2019-01-01','2021-12-31'))
		d_ii = d_ii.resample(time='1D').sum()
		d_ii = d_ii.values
		mean_ii.append(d_ii)
				
		# reading cmorph 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/cmorph/' + 'CMORPH_V1.0_ADJ_CSAM_8km_day_20180101-20211231.nc')
		d_iii = d_iii.cmorph.sel(time=slice('2019-01-01','2021-12-31'))
		d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iii = d_iii.values
		mean_iii.append(d_iii)

		# reading era5 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'tp_era5_csam_8km_day_20180101-20211231.nc')
		d_iv = d_iv.tp.sel(time=slice('2019-01-01','2021-12-31'))
		d_iv = d_iv.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
		d_iv = d_iv.values
		mean_iv.append(d_iv)

	return mean_i, mean_ii, mean_iii, mean_iv


def setBoxColors(bp):
    setp(bp['boxes'][0], color='gray')
    setp(bp['medians'][0], color='gray')
    setp(bp['caps'][0], color='gray')
    setp(bp['caps'][1], color='gray')
    setp(bp['whiskers'][0], color='gray')
    setp(bp['whiskers'][1], color='gray')
        
    setp(bp['boxes'][1], color='blue')
    setp(bp['medians'][1], color='blue')
    setp(bp['caps'][2], color='blue')
    setp(bp['caps'][3], color='blue')
    setp(bp['whiskers'][2], color='blue')
    setp(bp['whiskers'][3], color='blue')
    
    setp(bp['boxes'][2], color='green')
    setp(bp['medians'][2], color='green')
    setp(bp['caps'][4], color='green')
    setp(bp['caps'][5], color='green')
    setp(bp['whiskers'][4], color='green')
    setp(bp['whiskers'][5], color='green')
       
    setp(bp['boxes'][3], color='red')
    setp(bp['medians'][3], color='red')
    setp(bp['caps'][6], color='red')
    setp(bp['caps'][7], color='red')
    setp(bp['whiskers'][6], color='red')
    setp(bp['whiskers'][7], color='red')


dt = 'H_2018-01-01_2021-12-31'

print('Import dataset')
# Import dataset
mean_i, mean_ii, mean_iii, mean_iv = import_inmet(dt)			
	
clim_regcm = mean_i
clim_inmet = mean_ii
clim_cmorph = mean_iii
clim_era5 = mean_iv

list_hc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 0,
3, 1, 3, 1, 1, 1, 3, 1, 3, 4, 3, 1, 3, 1, 3, 3, 3, 3, 1, 3, 0, 3, 3, 3, 3,
3, 4, 3, 4, 0, 0, 0, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 3, 4,
3, 0, 4, 0, 4, 4, 4, 4, 4, 0, 4, 0, 0, 4, 0, 4, 4, 4, 4, 4, 4, 0, 4, 0, 2,
4, 4, 2, 4, 4, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 0, 0, 0, 2, 0, 2, 2, 2, 
0, 2, 2, 2, 2, 2, 4, 0]

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

cluster_i = [regcm_cluster_i, inmet_cluster_i, cmorph_cluster_i, era5_cluster_i]
cluster_ii = [regcm_cluster_ii, inmet_cluster_ii, cmorph_cluster_ii, era5_cluster_ii]
cluster_iii = [regcm_cluster_iii, inmet_cluster_iii, cmorph_cluster_iii, era5_cluster_iii]
cluster_iv = [regcm_cluster_iv, inmet_cluster_iv, cmorph_cluster_iv, era5_cluster_iv]
cluster_v = [regcm_cluster_v, inmet_cluster_v, cmorph_cluster_v, era5_cluster_v]

print('Plot figure')
# Plot figure
fig = plt.figure(figsize=(10, 4))
x = np.arange(0.5, 25 + 0.5)

A = [[12, 4, 16, 7], [2, 14, 6, 7], [2, 4, 16, 7], [2, 4, 16, 7]]
B = [[12, 4, 16, 7], [2, 14, 6, 7], [2, 4, 16, 7], [2, 4, 16, 7]]
C = [[12, 4, 16, 7], [2, 14, 6, 7], [2, 4, 16, 7], [2, 4, 16, 7]]
D = [[12, 4, 16, 7], [2, 14, 6, 7], [2, 4, 16, 7], [2, 4, 16, 7]]
E = [[12, 4, 16, 7], [2, 14, 6, 7], [2, 4, 16, 7], [2, 4, 16, 7]]

ax = fig.add_subplot(1, 1, 1)
bp = plt.boxplot(cluster_i, positions=[1, 2, 3, 4], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_ii, positions=[6, 7, 8, 9], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_iii, positions=[11, 12, 13, 14], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_iv, positions=[16, 17, 18, 19], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_v, positions=[21, 22, 23, 24], sym='.')
setBoxColors(bp)

plt.xlim(0, 25)
plt.ylim(0, 40)
plt.yticks(np.arange(0, 44, 4))
plt.xticks(x, ('','','','','','','','','','','','','','','','','','','','','','','','',''))
plt.xlabel('Clusters', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.axhspan(36, 40, color='gray', alpha=0.5, lw=0)
plt.text(1.25, 37, 'Cluster I', fontweight='bold')
plt.text(6.25, 37, 'Cluster II', fontweight='bold')
plt.text(11.25, 37, 'Cluster III', fontweight='bold')
plt.text(16.25, 37, 'Cluster IV', fontweight='bold')
plt.text(21.25, 37, 'Cluster V', fontweight='bold')
plt.axhline(36., linewidth=1., linestyle='-',  color='black')
plt.axvline(5., linewidth=1., linestyle='--',  color='black')
plt.axvline(10, linewidth=1., linestyle='--',  color='black')
plt.axvline(15, linewidth=1., linestyle='--',  color='black')
plt.axvline(20, linewidth=1., linestyle='--',  color='black')

c1, = plt.plot([1,1],'gray')
c2, = plt.plot([1,1],'blue')
c3, = plt.plot([1,1],'green')
c4, = plt.plot([1,1],'red')
plt.legend((c1, c2, c3, c4),('RegCM4', 'INMET', 'CMORPH', 'ERA5'), bbox_to_anchor=(0.5, 1.09), loc=9, ncol=4, frameon=False)
c1.set_visible(False)
c2.set_visible(False)
c3.set_visible(False)
c4.set_visible(False)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_stations_cluster_boxplot_pr_sesa.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()
