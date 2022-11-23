# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot annual cycle for subregions on sesa"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.grid_finder as GF
import mpl_toolkits.axisartist.floating_axes as FA

from scipy import stats
from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn
from matplotlib.projections import PolarAxes


class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd, fig=None, rect=224, label='_', srange=(0., 3.), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        self.refstd = refstd           
        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.])
        
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
            
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0,self.tmax, # 1st quadrant
                                                     self.smin,self.smax),
                                           grid_locator1=gl1,
                                           tick_formatter1=tf1,
                                           )
		
        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text(u'Correlation')

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text(u'SD')
        #~ ax.axis["bottom"].set_tick_params(labelsize=4)
       
        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
			"bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)      # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*', ls='', ms=10, label=label)
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]


    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        'Figure.plot' command.
        """

        l, = self.ax.plot(np.arccos(corrcoef), stddev, *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)
   
        return l
        

    def add_grid(self, *args, **kwargs):
        """Add a grid."""
        
        self._ax.grid(*args, **kwargs)
        

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """
        
        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax), np.linspace(0, self.tmax))
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))
        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)
        
        return contours
  
        
if __name__=='__main__':
	       
	subregion_i = []
	subregion_ii = []
	subregion_iii = []
	subregion_iv = []
	subregion_v = []
	subregion_vi = []
	subregion_vii = []
	subregion_viii = []
	subregion_ix = []

	era5_i = []
	era5_ii = []
	era5_iii = []
	era5_iv = []
	era5_v = []
	era5_vi = []
	era5_vii = []
	era5_viii = []
	era5_ix = []

	# Getting the data subregion I
	for j in range(1, 4):
		
		if j == 4:
			continue
		if j == 12:
			continue

		# Reading inmet weather station	
		ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[j][0]))
		ds_i = ds_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
		ds_i = ds_i.groupby('time.month').mean('time')
		subregion_i.append(ds_i.values*24)

		# reading era5 reanalisis
		ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
		ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
		ds_ii = ds_ii.sel(latitude=inmet[j][2], longitude=inmet[j][3], method='nearest')
		ds_ii = ds_ii.groupby('time.month').mean('time')
		era5_i.append(ds_ii.values*24)

	inmet_subregion_i = np.nanmean(subregion_i, axis=0)
	inmet_subregion_ii = np.nanmean(subregion_i, axis=0)
	inmet_subregion_iii = np.nanmean(subregion_i, axis=0)
	inmet_subregion_iv = np.nanmean(subregion_i, axis=0)
	inmet_subregion_v = np.nanmean(subregion_i, axis=0)
	inmet_subregion_vi = np.nanmean(subregion_i, axis=0)
	inmet_subregion_vii = np.nanmean(subregion_i, axis=0)
	inmet_subregion_viii = np.nanmean(subregion_i, axis=0)
	inmet_subregion_ix = np.nanmean(subregion_i, axis=0)
		
	era5_subregion_i = np.nanmean(era5_i, axis=0)
	era5_subregion_ii = np.nanmean(era5_i, axis=0)
	era5_subregion_iii = np.nanmean(era5_i, axis=0)
	era5_subregion_iv = np.nanmean(era5_i, axis=0)
	era5_subregion_v = np.nanmean(era5_i, axis=0)
	era5_subregion_vi = np.nanmean(era5_i, axis=0)
	era5_subregion_vii = np.nanmean(era5_i, axis=0)
	era5_subregion_viii = np.nanmean(era5_i, axis=0)
	era5_subregion_ix = np.nanmean(era5_i, axis=0)

	print('Plot figure')
	# Plot figure
	fig = plt.figure(figsize=(9, 11))
	
	# Reference database standard desviation
	stdrefs = dict(region_i = era5_subregion_i.std(ddof=1),
	region_ii  = era5_subregion_i.std(ddof=1),
	region_iii = era5_subregion_i.std(ddof=1),
	region_iv  = era5_subregion_i.std(ddof=1),
	region_v   = era5_subregion_i.std(ddof=1))
	
	text1 = dict(region_i ='(a) Subregion I',
	region_ii  ='(b) Subregion II',
	region_iii ='(c) Subregion III',
	region_iv  ='(d) Subregion IV',
	region_v   ='(e) Subregion IV')  
				 
	# Compute stddev and correlation coefficient of models
	samples = dict(region_i = [[inmet_subregion_i.std(ddof=1), np.corrcoef(era5_subregion_i, inmet_subregion_i)[0,1], 'ERA5']],
	region_ii  = [[inmet_subregion_i.std(ddof=1), np.corrcoef(era5_subregion_i, inmet_subregion_i)[0,1], 'ERA5']],
	region_iii = [[inmet_subregion_i.std(ddof=1), np.corrcoef(era5_subregion_i, inmet_subregion_i)[0,1], 'ERA5']],
	region_iv  = [[inmet_subregion_i.std(ddof=1), np.corrcoef(era5_subregion_i, inmet_subregion_i)[0,1], 'ERA5']],
	region_v   = [[inmet_subregion_i.std(ddof=1), np.corrcoef(era5_subregion_i, inmet_subregion_i)[0,1], 'ERA5']])		          

	# Here set placement of the points marking 95th and 99th significance
	# levels. For more than 102 samples (degrees freedom > 100), critical
	# correlation levels are 0.195 and 0.254 for 95th and 99th
	# significance levels respectively. Set these by eyeball using the
	# standard deviation x and y axis.

	# x95 = [0.01, 0.55] # For Tair, this is for 95th level (r = 0.195)
	# y95 = [0.0, 3.0]
	# x99 = [0.01, 0.95] # For Tair, this is for 99th level (r = 0.254)
	# y99 = [0.0, 3.0]
	
	rects = dict(region_i = 321,
	region_ii  = 322,
	region_iii = 323,
	region_iv  = 324,
	region_v   = 325)
				 	
	for var in ['region_i', 'region_ii', 'region_iii', 'region_iv', 'region_v']:

		dia = TaylorDiagram(stdrefs[var], fig=fig, rect=rects[var], label=u'Reference', srange=(0., 3.), extend=False)
		dia.samplePoints[0].set_color('black')
		# ~ dia.ax.plot(x95,y95,color='black')
		# ~ dia.ax.plot(x99,y99,color='black')

		colors = plt.matplotlib.cm.jet(np.linspace(0.299, 1,len(samples['region_i'])))

		# Add samples to Taylor diagram
		for i, (stddev,corrcoef,name) in enumerate(samples[var]):
			dia.add_sample(stddev, corrcoef,
						   marker='$%d$' % (i+1), ms=10, ls='', c=colors[i], label=name)
						   
			plt.text(0.09, 3., text1[var], fontweight='bold')

		# Add RMS contours, and label them
		contours = dia.add_contours(colors='0.5')
		plt.clabel(contours, inline=2, fontsize=10, fmt='%.1f')

		# Tricky: ax is the polar ax (used for plots), _ax is the container (used for layout)
		dia.add_grid()                                  
		dia._ax.axis[:].major_ticks.set_tick_out(True) 

	# Add a figure legend and title. For loc option, place x,y tuple inside [ ].
	# Can also use special options here: http://matplotlib.sourceforge.net/users/legend_guide.html

	# ~ # Add a figure legend
	# ~ fig.legend(dia.samplePoints, 
			   # ~ [ p.get_label() for p in dia.samplePoints ], 
			   # ~ prop=dict(size=10), numpoints=1, loc=(0.73, 1.10))

	# ~ plt.subplots_adjust(left=0.10, bottom=0.10, right=0.70, top=0.90, wspace=0.55, hspace=0.20)
		
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_taylor_diagram_all_subregions_2018-2021.png'
	if not os.path.exists(path_out):
		create_path(path_out)
	plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
	plt.show()
	exit()
