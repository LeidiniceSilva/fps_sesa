# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script plot time series to each inmet automatic station"

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

from dict_stations_urug_smn import urug_smn

dt = 'H_2018-01-01_2021-12-31'

for idx in range(1, 73):

	print()
	print('Reading Uruguai weather station:', idx, urug_smn[idx][0])	
	# Reading inmet weather station	"pr"
	pr_input = '/home/nice/Documentos/FPS_SESA/urug_smn/urug_smn_nc/pre_{0}_{1}.nc'.format(urug_smn[idx][0], dt)
	pr_input = nc.Dataset(pr_input)
	pre = pr_input['pre'][:]
	
	print('Plot figure')
	# Plot figure
	fig = plt.figure(figsize=(10, 6))
	time = np.arange(0.2, 35064 + 0.2)
	
	ax = fig.add_subplot(1, 1, 1)
	plt.plot(time, pre, linewidth=0.5, color='blue', label = 'INMET')
	plt.title(u'{0}'.format(urug_smn[idx][0]), fontsize=8, fontweight='bold')
	plt.xlabel(u'00h 01/01/2018 - 23h 31/12/2021', fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm h$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.ylim(-0.5, 20)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
		
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_nc_vars_{0}.png'.format(urug_smn[idx][0])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

