# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sept 22, 2025"
__description__ = "This script plot annual variability"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii


era5_c_i_dt = [3.5352616, 3.1525583, 3.276934,  4.9777637, 5.799324,  4.9096255, 5.6049438,
 8.071165,  7.049877,  5.075211,  4.6197033, 6.0924354, 1.8022686, 2.2555478,
 1.7588763, 2.5753777, 6.274767,  5.1259375, 5.173312,  6.2439437, 5.956787,
 2.3435526, 2.225257,  2.6337051, 4.6608396, 3.2362761, 1.9720353, 3.6633186,
 3.1959972, 5.2666807, 7.793757,  9.040763,  5.233349,  5.553249,  2.737926,
 3.6671624]

era5_c_i_dt_iii = np.nanmean(era5_c_i_dt)

# Plot figure
fig = plt.figure(figsize=(12, 8))
font_size = 8

dt = pd.date_range(start="20180601", end="20210531", freq="ME")

ax = fig.add_subplot(1, 1, 1)

era5_c_i_dt        = pd.Series(data=era5_c_i_dt, index=dt)
plt.plot(era5_c_i_dt,        linewidth=1, color='red', label = 'ERA5')

ax.text(0.2, 0.95, 'Reg4 = {0}({1})'.format(round(era5_c_i_dt_iii, 2), round(era5_c_i_dt_iii, 2)), transform=ax.transAxes, ha='right', va='top', fontsize=8, color='black')
ax.text(0.2, 0.9, 'Reg4 = {0}({1})'.format(round(era5_c_i_dt_iii, 2), round(era5_c_i_dt_iii, 2)), transform=ax.transAxes, ha='right', va='top', fontsize=8, color='black')

	
plt.title('(a) Cluster I', loc='left', fontsize=8, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)

plt.show()
exit()


