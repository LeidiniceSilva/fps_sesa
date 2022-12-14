# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script plot time series to each INMET automatic station"

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

from dict_stations_inmet import inmet

dt = 'H_2018-01-01_2021-12-31'

# Getting the data 
for j in range(1, 289):
	
	if j == 4:
		continue
	if j == 12:
		continue
	if j == 45:
		continue
	if j == 55:
		continue
	if j == 77:
		continue
	if j == 98:
		continue
	if j == 99:
		continue
	if j == 118:
		continue
	if j == 122:
		continue
	if j == 130:
		continue
	if j == 135:
		continue
	if j == 151:
		continue
	if j == 155:
		continue
	if j == 159:
		continue
	if j == 160:
		continue
	if j == 163:
		continue
	if j == 164:
		continue
	if j == 181:
		continue
	if j == 183:
		continue
	if j == 186:
		continue
	if j == 187:
		continue
	if j == 188:
		continue
	if j == 209:
		continue
	if j == 216:
		continue
	if j == 228:
		continue
	if j == 236:
		continue
	if j == 246:
		continue
	if j == 268:
		continue
	if j == 287:
		continue
			
	yy = inmet[j][2]
	xx = inmet[j][3]
	name_i = inmet[j][1]
	name_ii = inmet[j][1]

	print('Reading INMET weather station:', j, inmet[j][0], inmet[j][1])
	# Reading inmet weather station	"pr"
	pr_input = '/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/pre_{0}_{1}.nc'.format(inmet[j][0], dt)
	pr_input = nc.Dataset(pr_input)
	pre = pr_input['pre'][:]

	# Reading inmet weather station	"tp"
	tp_input = '/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/tmp_{0}_{1}.nc'.format(inmet[j][0], dt)
	tp_input = nc.Dataset(tp_input)
	tmp = tp_input['tmp'][:]

	# Reading inmet weather station	"uv"
	uv_input = '/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/uv_{0}_{1}.nc'.format(inmet[j][0], dt)
	uv_input = nc.Dataset(uv_input)
	uv = uv_input['uv'][:]
	
	print('Plot figure')
	# Plot figure
	fig = plt.figure(figsize=(10, 6))
	time = np.arange(0.2, 35064 + 0.2)
	
	ax = fig.add_subplot(3, 1, 1)
	plt.plot(time, pre, linewidth=0.5, color='blue', label = 'INMET')
	plt.title(u'{0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], name_ii, yy, xx), fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm h$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.ylim(-0.5, 20)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
		
	ax = fig.add_subplot(3, 1, 2)
	plt.plot(time, tmp, linewidth=0.5, color='blue', label = 'INMET')
	plt.title(u'{0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], name_ii, yy, xx), fontsize=8, fontweight='bold')
	plt.ylabel(u'Temperature (°C)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 36, 3), fontsize=8)
	plt.ylim(-0.5, 33)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
	
	ax = fig.add_subplot(3, 1, 3)
	plt.plot(time, uv, linewidth=0.5, color='blue', label = 'INMET')
	plt.title(u'{0} - {1} lat: {2} lon: {3}'.format(inmet[j][0], name_ii, yy, xx), fontsize=8, fontweight='bold')
	plt.xlabel(u'00h 01/01/2018 - 23h 31/12/2021', fontsize=8, fontweight='bold')
	plt.ylabel(u'Wind speed (m s$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.xticks(fontsize=8)
	plt.ylim(-0.5, 20)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
			
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_nc_vars_{0}_{1}.png'.format(inmet[j][0], name_i)
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

