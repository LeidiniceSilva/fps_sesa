# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot ts of the INMET weather station"

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

from dict_sesa_inmet_stations import inmet

dt = 'H_2018-01-01_2021-12-31'

# Getting the data 
for j in range(1, 101):

	print('Reading INMET weather station:', j, inmet[j][0], inmet[j][1])
	
	# Reading inmet weather station	"pr"
	pre_input = '/home/nice/Documentos/FPS_SESA/docs/inmet_nc_sesa/pre/pre_{0}_{1}.nc'.format(inmet[j][0], dt)
	pre_input = nc.Dataset(pre_input)
	pre = pre_input['pre'][:]

	# Reading inmet weather station	"tmp"
	tmp_input = '/home/nice/Documentos/FPS_SESA/docs/inmet_nc_sesa/tmp/tmp_{0}_{1}.nc'.format(inmet[j][0], dt)
	tmp_input = nc.Dataset(tmp_input)
	tmp = tmp_input['tmp'][:]

	# Reading inmet weather station	"dtmp"
	dtmp_input = '/home/nice/Documentos/FPS_SESA/docs/inmet_nc_sesa/dtmp/dtmp_{0}_{1}.nc'.format(inmet[j][0], dt)
	dtmp_input = nc.Dataset(dtmp_input)
	dtmp = dtmp_input['dtmp'][:]
	
	# Reading inmet weather station	"uv"
	uv_input = '/home/nice/Documentos/FPS_SESA/docs/inmet_nc_sesa/uv/uv_{0}_{1}.nc'.format(inmet[j][0], dt)
	uv_input = nc.Dataset(uv_input)
	uv = uv_input['uv'][:]
	
	print('Plot figure')
	# Plot figure
	fig = plt.figure(figsize=(10, 10))
	time = np.arange(0.2, 35064 + 0.2)
	
	ax = fig.add_subplot(4, 1, 1)
	plt.plot(time, pre, linewidth=0.5, color='blue', label = 'PRE')
	plt.title(u'(a) {0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], inmet[j][1], inmet[j][2], inmet[j][3]), loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm h⁻¹)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.ylim(-0.5, 20)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
		
	ax = fig.add_subplot(4, 1, 2)
	plt.plot(time, tmp, linewidth=0.5, color='red', label = 'TMP')
	plt.title(u'(b) {0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], inmet[j][1], inmet[j][2], inmet[j][3]), loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Temperature (°C)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 38, 3), fontsize=8)
	plt.ylim(-0.5, 35)
	plt.legend(loc=1, fontsize=8)
	plt.grid()

	ax = fig.add_subplot(4, 1, 3)
	plt.plot(time, dtmp, linewidth=0.5, color='red', label = 'DTMP')
	plt.title(u'(c) {0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], inmet[j][1], inmet[j][2], inmet[j][3]), loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Dewpoint temperature (°C)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 38, 3), fontsize=8)
	plt.ylim(-0.5, 35)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
	
	ax = fig.add_subplot(4, 1, 4)
	plt.plot(time, uv, linewidth=0.5, color='green', label = 'UV')
	plt.title(u'(d) {0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], inmet[j][1], inmet[j][2], inmet[j][3]), loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'00h 01/01/2018 - 23h 31/12/2021', fontsize=8, fontweight='bold')
	plt.ylabel(u'Wind speed (m s⁻¹)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.xticks(fontsize=8)
	plt.ylim(-0.5, 16)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
			
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_ts'
	name_out = 'pyplt_nc_vars_{0}_{1}.png'.format(inmet[j][0], inmet[j][1])
	plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

