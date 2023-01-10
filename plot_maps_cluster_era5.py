# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "01/02/2023"
__description__ = "This script plot cluster analysis from each weather station"

import os
import conda
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

type_cycle = 'annual_cycle'

longitude = [-53.37138888, -53.37638888, -52.1, -52.70111111, -54.01333333, -50.90555555, -54.61805554, -55.61277777, -51.83472221, -53.46749999, -56.43722221, -54.31083333, -51.16666666, -50.13527777, -53.82666666, -52.38249999, -57.08249999, -55.5261111, -53.70361111, -54.69416666, -51.06416666, -51.82416666, -50.82749999, -49.73305554, -54.88555555, -51.53472221, -50.14944444, -49.49805555, -52.55833333, -50.05833333, -51.87055554, -53.11194444, -56.01555555, -48.81333333, -53.6736111, -49.31555555, -50.88277777, -54.96222222, -49.93444444, -52.40388888, -51.51222222, -49.48333333, -53.31694443, -54.48055554, -53.7911111, -50.33555555, -49.04194444, -52.30638888, -48.61999999, -49.64666666, -53.41277777, -50.60416666, -51.57555554, -48.76194444, -52.415, -50.14638888, -49.26749999, -50.98555555, -53.50444444, -52.34888888, -52.85055555, -51.35444444, -50.36333333, -53.63277777, -49.58055554, -48.64166666, -50.3686111, -53.74805555, -53.09444444, -54.48305554, -51.08944444, -48.80861111, -49.26666666, -52.39194444, -49.1575, -50.8711111, -48.41638888, -49.99972221, -47.55, -54.01972221, -47.86416666, -51.96305555, -50.24555555, -48.88527777, -55.02416666, -46.53055554, -48.16444444, -46.14305555, -45.4025, -50.18055554, -49.94638888, -46.61666666, -54.18166666, -51.91666666, -50.57777777, -53.63472221, -52.93194443, -47.66666666, -45.41722221, -44.72666666, -48.94555555, -43.6, -45.52027777, -55.32944443, -43.19027777, -42.01666666, -44.30333333, -49.89416666, -43.40277777, -42.60916666, -46.04333333, -43.41111111, -54.60555555, -43.68333333, -45.60388888, -47.62305554, -45.00555555, -44.04083333, -42.41583333, -52.85888888, -43.28222221, -55.71638888, -52.13444444, -43.29138888, -44.45, -42.98722221, -46.80527777, -44.96166666, -41.81666667, -44.70305555, -50.97416666, -48.55722222, -49.02888888, -43.7, -42.67722221, -45.37277777, -53.82722221, -49.965, -54.91138888, -51.4, -56.53999999, -43.20861111, -41.05194444, -47.88333333, -50.49027777, -46.38305554, -48.66666666, -47.07972221, -54.52805554, -43.36416666, -52.47055554, -41.35, -57.88666666, -45.94444444, -49.73416666, -55.17777777, -41.95, -45.40416666, -43.26111111, -56.13694444, -51.55222222, -48.11388888, -50.93027777, -43.76694443, -44.97944443, -48.84027777, -42.37583332, -41.03888888, -49.92055554, -54.97194443, -48.48972221, -47.11416666, -51.71222222, -42.86388888, -41.48888888, -46.63388888, -44.86444444, -40.74138888, -47.37999999, -48.54472221, -43.76944443, -55.78388888, -40.40416666, -45.45388888, -54.6, -52.87555554, -49.96583333, -56.41666666, -47.77499999, -40.3, -42.18277777, -41.18999999, -44.87472221, -50.59499999, -41.10694444, -44.01111111, -46.00861111, -40.57944443, -48.15138888, -43.95861111, -44.41694444, -43.96944443, -47.43416666, -42.15361111, -47.96166666, -51.18138888, -46.94944444, -42.62222221, -49.51805554, -41.09083333, -45.59388888, -44.17333332, -54.55305555, -40.53972221, -40.0686111, -57.63749999, -46.98583333, -56.62305554, -50.61666666, -49.52499999, -48.25555555, -41.97666666, -53.26722221, -52.60277777, -42.94305555, -40.98583333, -44.45361111, -39.84833333, -40.4075, -46.44055555, -49.19194444, -40.73638888, -43.64805555, -45.45972222, -47.92638888, -38.6961111, -51.71777777, -41.51666666, -54.45361111, -40.24972221, -50.98138888, -46.11916666, -39.26666666, -49.1, -42.38916666, -47.19944444, -52.60111111, -53.22416666, -49.91472222, -48.28416666, -54.83722221, -44.83555555, -39.55138888, -56.542211, -56.224017, -55.866233, -56.555522, -57.4559, -57.485838, -57.097101, -56.618974, -56.456119, -58.165178, -57.325281, -57.78309, -56.972142, -57.429814, -56.214658, -56.326675, -58.295964, -56.657495, -56.629005, -58.141889, -57.021887, -57.85945, -57.680318, -56.230152, -56.905948, -57.421509, -58.070281, -57.903, -57.055016, -54.193437, -58.247317, -57.880716, -57.913847, -55.642269, -56.635879, -57.703535, -56.787211, -57.08769, -58.468999, -56.411658, -57.743836, -57.674446, -57.820722, -57.966907, -57.982836, -57.622519, -56.786308, -57.291694, -57.082856, -57.074737, -56.66754, -57.124992, -56.84021, -53.933283, -58.002561, -57.880984, 493837, -57.939, -56.392531, -57.977472, -57.916594, -58.317145, -55.129216, -57.565952, -57.92773, -56.028995, -56.215761, -56.879045, -58.193534, -57.440628, -57.480465, -60.2483, -59.2162, -59.5456, -59.436, -58.5787, -60.0885, -59.5572, -59.9712, -59.8411, -58.615, -59.209, -60.4211, -60.1331, -60.3219, -58.4125, -60.6044, -59.0099, -59.9312, -59.4519, -60.1917, -58.5287, -58.6634, -60.5174, -58.7863, -60.3981, -58.6356, -59.8476, -59.8144, -59.9075, -58.7729, -58.1467, -60.3631, -60.3077, -60.1638, -59.3818, -59.5882, -59.2999, -59.8915, -60.0818, -59.9056, -59.2894, -59.3278, -59.5185, -58.879, -58.9821, -59.0893, -57.9219, -58.083, -58.3011, -58.1921, -60.528, -58.4742, -60.1094, -59.524, -58.2346, -59.6821, -59.6263, -58.7438, -58.914, -59.3869, -59.161, -59.1903, -59.5917, -59.6215, -58.9641, -60.0252, -60.3073, -59.7595, -60.5324, -59.5292, -59.122]
latitude = [-33.74166666, -32.53749999, -32.03333333, -31.40583333, -31.34777777, -31.24777777, -31.0025, -30.84249999, -30.81055555, -30.54777777, -30.3686111, -30.34138888, -30.05, -30.00972222, -29.89388888, -29.87333332, -29.84249999, -29.71166666, -29.70833333, -29.70194444, -29.67444443, -29.45027777, -29.36888888, -29.35027777, -29.19138888, -29.16722221, -29.04888888, -28.93138888, -28.85361111, -28.75138888, -28.70472222, -28.65333333, -28.64944444, -28.60416666, -28.60361111, -28.53249999, -28.51361111, -28.41694443, -28.27555554, -28.22944443, -28.22194443, -28.13333333, -27.92166666, -27.89305555, -27.85416666, -27.80222222, -27.6786111, -27.66027777, -27.6025, -27.41833332, -27.39555555, -27.2886111, -27.16944443, -26.95083333, -26.9386111, -26.93749999, -26.91361111, -26.81944443, -26.77666666, -26.41722221, -26.40694444, -26.3986111, -26.3936111, -26.28638888, -26.2486111, -26.0811111, -25.83499999, -25.72194443, -25.69472221, -25.60194444, -25.56583333, -25.51277777, -25.43333333, -25.36888888, -25.32222221, -25.01333333, -24.96305555, -24.78944444, -24.71666666, -24.53583333, -24.53305554, -24.43722221, -24.23833333, -23.98138888, -23.96694443, -23.95194444, -23.85138888, -23.84499999, -23.81083333, -23.77305554, -23.50527777, -23.48333333, -23.44944444, -23.40777777, -23.40694444, -23.39027777, -23.37583332, -23.35, -23.22833332, -23.22333332, -23.09972221, -23.05, -23.04194444, -23.0025, -22.98833333, -22.98333333, -22.97583332, -22.9486111, -22.93972221, -22.8711111, -22.86138888, -22.86083333, -22.85722222, -22.8, -22.75027777, -22.70277777, -22.68888888, -22.65361111, -22.64583333, -22.63388888, -22.58972221, -22.5525, -22.49166666, -22.48166666, -22.45, -22.44888888, -22.415, -22.3961111, -22.38333333, -22.37388888, -22.37249999, -22.37083332, -22.35805555, -22.35, -22.33305554, -22.31444444, -22.30666666, -22.23527777, -22.19388888, -22.11999999, -22.10083333, -22.09833333, -22.04166666, -21.97972221, -21.92722221, -21.91805554, -21.85555555, -21.77972221, -21.77472221, -21.76999999, -21.75138888, -21.71666666, -21.70583333, -21.68083332, -21.66527777, -21.60916666, -21.56666666, -21.56638888, -21.54694444, -21.47833332, -21.45777777, -21.33833333, -21.31916666, -21.22888888, -21.22638888, -21.13305554, -21.10472222, -21.10111111, -21.08555555, -20.98166666, -20.94916666, -20.91, -20.78999999, -20.76277777, -20.75055555, -20.74499999, -20.71472222, -20.63638888, -20.57999999, -20.55888888, -20.54583333, -20.47555554, -20.46694443, -20.45472222, -20.45, -20.44444444, -20.40305555, -20.38333333, -20.35944444, -20.26666666, -20.26361111, -20.25222222, -20.17305554, -20.165, -20.10416666, -20.0311111, -20.0311111, -19.98833333, -19.98583333, -19.97999999, -19.88527778, -19.88416666, -19.87527777, -19.73583333, -19.71, -19.69527777, -19.60555555, -19.5736111, -19.5391666, -19.5327777, -19.48166666, -19.455, -19.42027777, -19.40722222, -19.35694444, -18.99666666, -18.99638888, -18.98888888, -18.96666666, -18.95277777, -18.91722221, -18.83027777, -18.82277777, -18.80222222, -18.78666666, -18.78027777, -18.76444444, -18.71388888, -18.69527777, -18.52055554, -18.40972222, -18.29166666, -18.2311111, -18.20055555, -18.15777777, -17.96416666, -17.92388888, -17.9, -17.89638888, -17.7986111, -17.78583333, -17.78472221, -17.73333333, -17.71666667, -17.70527777, -17.56166666, -17.45472222, -17.33944444, -17.33694444, -17.30416666, -17.29722221, -17.25777777, -17.00722222, -28.416682, -27.416682, -26.416682, -25.416682, -24.416682, -23.416682, -22.416682, -21.416682, -20.416682, -19.416682, -18.416682, -30.787252, -30.299302, -29.819517, -30.751196, -30.637259, -29.98349, -30.62427, -31.337833, -30.776241, -31.09242, -30.014261, -30.580333, -30.58639, -30.613837, -31.444407, -29.787069, -31.459, -30.755734, -27.298582, -30.971026, -30.971187, -30.992189, -28.177515, -31.037601, -31.286127, -30.432802, -29.05663, -30.590296, -30.865546, -29.660671, -29.845637, -30.690307, -30.214445, -30.626772, -30.249296, -30.787214, -30.289, -29.721745, -30.205281, -30.223953, -30.476397, -31.444635, -27.153304, -31.404212, -30.390258, -31, -32.135, -31.350331, -31.38978, -31.272706, -30.337372, -27.869343, -29.296744, -30.91468, -28.545073, -30.988557, -31.004877, -29.377358, -30.431886, -28.691829, -32.0706, -31.819, -32.409, -33.0357, -32.5229, -32.5708, -30.7405, -31.2033, -31.5181, -33.0328, -31.555, -31.9183, -31.9648, -31.7096, -32.1482, -32.054, -33.0512, -32.2483, -30.6953, -32.1525, -32.0017, -30.286, -31.8471, -30.9272, -32.3279, -31.7134, -31.9651, -32.4139, -32.7802, -30.6597, -32.2541, -32.2205, -32.4834, -32.589, -32.7225, -31.4774, -33.1011, -32.583, -31.5899, -31.69, -31.1092, -31.4762, -31.1867, -32.3825, -30.3711, -31.3108, -31.0585, -31.2972, -31.9419, -31.8041, -31.7282, -31.8546, -31.7531, -30.5624, -30.3611, -30.9942, -32.5596, -32.9224, -31.3839, -30.8536, -32.0626, -32.3259, -30.7977, -32.8519, -31.6438, -31.8577, -31.6222, -32.7458, -32.1278, -32.0201, -32.3807]

if type_cycle == 'diurnal_cycle':
	list_hc = [2, 4, 2, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 0, 4, 4, 0, 0, 4, 4, 0, 4, 4, 4, 2, 4, 0, 4,
	4, 1, 4, 4, 0, 4, 4, 4, 0, 0, 4, 4, 4, 4, 1, 4, 4, 4, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 4, 1, 0, 1, 0,
	1, 1, 0, 1, 4, 1, 4, 0, 1, 3, 1, 0, 3, 0, 0, 3, 3, 1, 1, 1, 1, 1, 1, 3, 1, 0, 3, 3, 1, 1, 3, 2, 0, 3, 3, 0, 0,
	1, 1, 0, 0, 3, 1, 0, 0, 1, 0, 1, 3, 0, 0, 0, 1, 1, 3, 0, 3, 3, 3, 0, 0, 0, 3, 3, 3, 3, 1, 1, 3, 3, 3, 0, 3, 3,
	3, 0, 3, 1, 2, 1, 3, 1, 3, 3, 0, 1, 3, 3, 3, 0, 3, 3, 0, 1, 3, 1, 3, 1, 3, 3, 1, 1, 1, 0, 1, 3, 0, 1, 4, 3, 1,
	3, 3, 1, 3, 2, 1, 1, 3, 3, 3, 1, 1, 0, 1, 1, 3, 1, 1, 3, 1, 3, 1, 2, 3, 3, 3, 3, 1, 3, 4, 3, 1, 3, 3, 3, 1, 2,
	1, 1, 3, 3, 3, 2, 3, 1, 3, 3, 1, 3, 3, 2, 1, 3, 1, 3, 1, 3, 2, 3, 1, 1, 1, 1, 3, 3, 1, 2, 3, 4, 4, 1, 4, 2, 2,
	4, 1, 1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
else:
	list_hc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1,
	1, 1, 1, 1, 2, 1, 1, 1, 1, 4, 1, 4, 4, 1, 4, 1, 4, 1, 2, 4, 4, 1, 4, 1, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4,
	4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4, 2, 0, 2, 2, 0, 0, 0, 4, 4, 0, 4, 4, 0, 0, 2, 0, 0, 0, 4, 0, 0, 2, 0, 0, 2, 2,
	0, 4, 2, 2, 0, 0, 2, 2, 0, 2, 4, 0, 2, 2, 2, 0, 2, 0, 2, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0,
	0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 4, 4,
	4, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3,
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3,
	1, 1, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
	 
print(len(longitude))
print(len(latitude))
print(len(list_hc))

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

lon_c_i = []
lon_c_ii = []
lon_c_iii = []
lon_c_iv = []
lon_c_v = []

lat_c_i = []
lat_c_ii = []
lat_c_iii = []
lat_c_iv = []
lat_c_v = []

for c_i in count_i:
	lon_c_i.append(longitude[c_i])
	lat_c_i.append(latitude[c_i])

for c_ii in count_ii:
	lon_c_ii.append(longitude[c_ii])
	lat_c_ii.append(latitude[c_ii])
	
for c_iii in count_iii:
	lon_c_iii.append(longitude[c_iii])
	lat_c_iii.append(latitude[c_iii])
	
for c_iv in count_iv:
	lon_c_iv.append(longitude[c_iv])
	lat_c_iv.append(latitude[c_iv])
	
for c_v in count_v:
	lon_c_v.append(longitude[c_v])
	lat_c_v.append(latitude[c_v])

print('Plot figure')
# Plot figure   
fig = plt.figure()

my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
my_map.drawmeridians(np.arange(-75.,-25.,5.), labels=[0,0,0,1], linewidth=0.5, color='black')
my_map.drawparallels(np.arange(-40.,5.,5.), labels=[1,0,0,0], linewidth=0.5, color='black') 
my_map.plot(lon_c_i, lat_c_i, 'o', color='blue', label='Cluster I', markersize=2)
my_map.plot(lon_c_ii, lat_c_ii, 'o', color='gray', label='Cluster II', markersize=2)
my_map.plot(lon_c_iii, lat_c_iii, 'o', color='green', label='Cluster III', markersize=2)
my_map.plot(lon_c_iv, lat_c_iv, 'o', color='red', label='Cluster IV', markersize=2)
my_map.plot(lon_c_v, lat_c_v, 'o', color='yellow', label='Cluster V', markersize=2)
plt.legend(loc=2, fontsize=10)

path = '/home/nice/Documentos/github_projects/shp'
my_map.readshapefile('{0}/lim_pais/lim_pais'.format(path), 'world', drawbounds=True, color='black', linewidth=0.5)

plt.title('Cluster Analysis')
plt.text(-68, -36, u'SESA', fontsize=10)
plt.text(-42, -36, u'\u25B2 \nN', fontsize=10, fontweight='bold')

x1,i1 = my_map(-70,-38)
x2,i2 = my_map(-70,-14)
x3,i3 = my_map(-38,-14)
x4,i4 = my_map(-38,-38)

poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_era5_cluster_analysis_{0}.png'.format(type_cycle)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.close('all')
plt.cla()
exit()

