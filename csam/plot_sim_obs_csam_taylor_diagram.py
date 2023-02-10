# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot taylor diagram in the csam  domain"

import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

from taylor_diagram import TaylorDiagram

stddev = np.random.uniform(low=1, high=10, size=(10,))  # Generate 10 random numbers between 1 and 10
corrcoeff = np.random.uniform(low=0.5, high=1, size=(10,))  # Generate 10 random numbers between 0.5 and 1
refstd = 5
models = ['model '+str(i) for i in range(1,11)]

print('Plot figure')
# Plot figure
fig = plt.figure(figsize=(12,8))

titleprops_dict = dict(loc='left',
                       fontweight='bold',
                       x=0.0, 
                       y=1.05)

fig, ax1 = TaylorDiagram(stddev, corrcoeff, refstd, 
                         fig=fig, 
                         rect=231,  # 2 row, 2 column, 1st plot
                         title='(a) Cluster I', titleprops_dict=titleprops_dict,
                         normalize=True, 
                         labels=models, ref_label='Reference')

fig, ax2 = TaylorDiagram(stddev, corrcoeff, refstd, 
                         fig=fig, 
                         rect=232,  # 2 row, 2 column, 2nd plot
                         title='(b) Cluster II', titleprops_dict=titleprops_dict,
                         normalize=True, 
                         labels=models, ref_label='Reference')

fig, ax3 = TaylorDiagram(stddev, corrcoeff, refstd, 
                         fig=fig, 
                         rect=233,  # 2 row, 2 column, 3rd plot
                         title='(c) Cluster III', titleprops_dict=titleprops_dict,
                         normalize=True,  
                         labels=models, ref_label='Reference')

fig, ax4 = TaylorDiagram(stddev, corrcoeff, refstd, 
                         fig=fig, 
                         rect=234,  # 2 row, 2 column, 4th plot,
                         title='(d) Cluster IV', titleprops_dict=titleprops_dict,
                         normalize=True, 
                         labels=models, ref_label='Reference')

fig, ax5 = TaylorDiagram(stddev, corrcoeff, refstd, 
                         fig=fig, 
                         rect=235,  # 2 row, 2 column, 4th plot,
                         title='(e) Cluster V', titleprops_dict=titleprops_dict,
                         normalize=True, 
                         labels=models, ref_label='Reference')
                        
ax5.legend(bbox_to_anchor=(1.2, 0.0), loc='lower left', ncol=1)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/csam'
name_out = 'pyplt_stations_cluster_taylor_diagram_pr_csam.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
	
