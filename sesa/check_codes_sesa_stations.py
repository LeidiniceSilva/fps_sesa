# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot cluster analysis of the INMET weather station"

import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet
from dict_urug_smn_stations import urug_smn

longitude = [-53.37138888, -53.37638888, -52.1, -52.70111111, -54.01333333, -54.61805554, -55.61277777, -51.83472221, -53.46749999, -56.43722221, -54.31083333, -51.16666666, -53.82666666, -52.38249999, -57.08249999, -55.5261111, -53.70361111, -54.69416666, -51.06416666, -50.82749999, -49.73305554, -54.88555555, -51.53472221, -50.14944444, -49.49805555, -52.55833333, -50.05833333, -51.87055554, -53.11194444, -56.01555555, -48.81333333, -53.6736111, -49.31555555, -50.88277777, -54.96222222, -49.93444444, -52.40388888, -51.51222222, -49.48333333, -53.31694443, -54.48055554, -53.7911111, -50.33555555, -49.04194444, -52.30638888, -48.61999999, -49.64666666, -53.41277777, -50.60416666, -51.57555554, -48.76194444, -52.415, -49.26749999, -50.98555555, -52.34888888, -51.35444444, -53.63277777, -49.58055554, -50.3686111, -51.08944444, -48.80861111, -49.26666666, -52.39194444, -49.1575, -50.8711111, -54.01972221, -48.88527777, -48.16444444, -50.18055554, -49.94638888, -54.18166666, -51.91666666, -52.93194443, -52.13444444, -50.97416666, -49.02888888, -49.965, -51.4, -50.49027777, -54.52805554, -49.73416666, -51.55222222, -48.11388888, -50.93027777, -48.54472221, -54.6, -52.87555554, -49.96583333, -48.15138888, -49.51805554, -57.63749999, -48.25555555, -39.84833333, -49.19194444, -51.71777777, -49.1, -52.60111111, -53.22416666, -49.91472222, -48.28416666, -56.542211, -56.224017, -55.866233, -56.555522, -57.4559, -57.485838, -57.097101, -56.618974, -56.456119, -58.165178, -57.325281, -57.78309, -56.972142, -57.429814, -56.214658, -56.326675, -58.295964, -56.657495, -56.629005, -58.141889, -57.021887, -57.85945, -57.680318, -56.230152, -56.905948, -57.421509, -58.070281, -57.903, -57.055016, -54.193437, -58.247317, -57.880716, -57.913847, -55.642269, -56.635879, -57.703535, -56.787211, -57.08769, -58.468999, -56.411658, -57.743836, -57.674446, -57.820722, -57.966907, -57.982836, -57.622519, -56.786308, -57.291694, -57.082856, -57.074737, -56.66754, -57.124992, -56.84021, -53.933283, -58.002561, -57.880984, 493837, -57.939, -56.392531, -57.977472, -57.916594, -58.317145, -55.129216, -57.565952, -57.92773, -56.028995, -56.215761, -56.879045, -58.193534, -57.440628, -57.480465]
latitude = [-33.74166666, -32.53749999, -32.03333333, -31.40583333, -31.34777777, -31.0025, -30.84249999, -30.81055555, -30.54777777, -30.3686111, -30.34138888, -30.05, -29.89388888, -29.87333332, -29.84249999, -29.71166666, -29.70833333, -29.70194444, -29.67444443, -29.36888888, -29.35027777, -29.19138888, -29.16722221, -29.04888888, -28.93138888, -28.85361111, -28.75138888, -28.70472222, -28.65333333, -28.64944444, -28.60416666, -28.60361111, -28.53249999, -28.51361111, -28.41694443, -28.27555554, -28.22944443, -28.22194443, -28.13333333, -27.92166666, -27.89305555, -27.85416666, -27.80222222, -27.6786111, -27.66027777, -27.6025, -27.41833332, -27.39555555, -27.2886111, -27.16944443, -26.95083333, -26.9386111, -26.91361111, -26.81944443, -26.41722221, -26.3986111, -26.28638888, -26.2486111, -25.83499999, -25.56583333, -25.51277777, -25.43333333, -25.36888888, -25.32222221, -25.01333333, -24.53583333, -23.98138888, -23.85138888, -23.77305554, -23.50527777, -23.44944444, -23.40777777, -23.37583332, -22.49166666, -22.37249999, -22.35805555, -22.23527777, -22.11999999, -21.92722221, -21.77472221, -21.66527777, -21.45777777, -21.33833333, -21.31916666, -20.55888888, -20.45, -20.44444444, -20.40305555, -19.98583333, -19.5391666, -18.99666666, -18.91722221, -18.71388888, -18.40972222, -17.92388888, -17.71666667, -17.45472222, -17.33944444, -17.33694444, -17.30416666, -28.416682, -27.416682, -26.416682, -25.416682, -24.416682, -23.416682, -22.416682, -21.416682, -20.416682, -19.416682, -18.416682, -30.787252, -30.299302, -29.819517, -30.751196, -30.637259, -29.98349, -30.62427, -31.337833, -30.776241, -31.09242, -30.014261, -30.580333, -30.58639, -30.613837, -31.444407, -29.787069, -31.459, -30.755734, -27.298582, -30.971026, -30.971187, -30.992189, -28.177515, -31.037601, -31.286127, -30.432802, -29.05663, -30.590296, -30.865546, -29.660671, -29.845637, -30.690307, -30.214445, -30.626772, -30.249296, -30.787214, -30.289, -29.721745, -30.205281, -30.223953, -30.476397, -31.444635, -27.153304, -31.404212, -30.390258, -31, -32.135, -31.350331, -31.38978, -31.272706, -30.337372, -27.869343, -29.296744, -30.91468, -28.545073, -30.988557, -31.004877, -29.377358, -30.431886, -28.691829]

list_hc = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 0, 2, 3, 3, 2, 0,
3, 0, 0, 3, 3, 0, 2, 0, 0, 3, 3, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1,
1, 1, 1, 0, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 4, 2, 4, 2, 4, 2, 4, 4, 4, 4,
4, 4, 2, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 2,
3, 2, 2, 2, 2, 1, 2]
	
inmet_list = list_hc[:100]
smn_list = list_hc[100:]

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

for count, idx in enumerate(smn_list):
	
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

list_inmet = []
for i in range(1, 101):
	list_inmet.append(inmet[i][0])

list_smn = []
for j in range(1, 72):
	list_smn.append(urug_smn[j][0])

code_c_i = []
code_c_ii = []
code_c_iii = []
code_c_iv = []
code_c_v = []

for c_i in count_i:
	code_c_i.append(list_smn[c_i])

for c_ii in count_ii:
	code_c_ii.append(list_smn[c_ii])
	
for c_iii in count_iii:
	code_c_iii.append(list_smn[c_iii])
	
for c_iv in count_iv:
	code_c_iv.append(list_smn[c_iv])
	
for c_v in count_v:
	code_c_v.append(list_smn[c_v])

print()
print("List of C.I: ", code_c_i)
print()
print("List of C.II: ", code_c_ii)
print()
print("List of C.III: ", code_c_iii)
print()
print("List of C.IV: ", code_c_iv)
print()
print("List of C.V: ", code_c_v)
print()
exit()

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

print("Longitude x Latitude C.I: ", len(lon_c_i))
print("Longitude x Latitude C.II: ", len(lon_c_ii))
print("Longitude x Latitude C.III: ", len(lon_c_iii))
print("Longitude x Latitude C.IV: ", len(lon_c_iv))
print("Longitude x Latitude C.V: ", len(lon_c_v))
print()

print("List of longitude C.I: ", lon_c_i)
print()
print("List of longitude C.II: ", lon_c_ii)
print()
print("List of longitude C.III: ", lon_c_iii)
print()
print("List of longitude C.IV: ", lon_c_iv)
print()
print("List of longitude C.V: ", lon_c_v)
print()

print("List of latitude C.I: ", lat_c_i)
print()
print("List of latitude C.II: ", lat_c_ii)
print()
print("List of latitude C.III: ", lat_c_iii)
print()
print("List of latitude C.IV: ", lat_c_iv)
print()
print("List of latitude C.V: ", lat_c_v)

# ~ # InMet
# ~ List of C.I:  ['A845', 'A870', 'A806', 'A863', 'A841', 'A858', 'A817', 'A848', 'A862', 'A823', 'A873', 'A807', 'B804', 'B806', 'A818', 'A820', 'A835']
# ~ List of C.II:  ['A714', 'A715', 'A871', 'A821', 'A752', 'A869', 'A718', 'A705', 'A763', 'A707', 'A768', 'A727', 'A762', 'A747', 'A734', 'A748', 'A702', 'A756', 'A729', 'A520', 'A519', 'A724', 'A507', 'A616', 'A035', 'A016', 'A003', 'A026', 'A909', 'A029', 'A033']
# ~ List of C.III:  ['A830', 'A856', 'A865', 'A868', 'A875', 'A874', 'A850', 'A743']
# ~ List of C.IV:  ['A899', 'A836', 'A802', 'A811', 'A827', 'A881', 'A804', 'A838', 'A812', 'A831', 'A832', 'A801', 'A886', 'A813', 'A809', 'A826', 'A803', 'A889', 'A884', 'A879', 'A808', 'A833', 'A840', 'A897', 'A867', 'A837', 'A829', 'A894', 'A883', 'A866', 'A853', 'A814', 'A880', 'A852', 'A815', 'A839', 'A844', 'A810', 'A805', 'A828', 'A854', 'A860', 'A859', 'A876']
# ~ List of C.V:  []

# ~ # SMN
# ~ List of C.I:  ['Puntas_de_Valentín']
# ~ List of C.II:  ['Baibene', 'Diego_Lamas', 'Tomás_Gomensoro']
# ~ List of C.III:  ['Aguapey', 'Alto_Arapey', 'Alto_Cuareim', 'Alvear', 'Arapey_Chico', 'Arapey_Grande', 'Arapey_Grande_Ruta_4', 'Arerunguá', 'Artigas', 'Belén', 'Bernabé_Rivera', 'Bonpland', 'Catalán', 'Catalán_Grande_rio', 'Cazadores_Correntinos', 'Cerro_Amarillo', 'Cerro_Chato', 'Chajarí', 'Colonia_Lavalleja', 'Colonia_Libertad', 'Colonia_Palma', 'Cuareim_Río', 'Cuaró', 'Cuchilla_de_Salto', 'Curuzú_Cuatiá', 'Dayman_Ruta_3', 'El_Soberbio', 'Espinillar', 'Garruchos', 'Itapebí', 'Mocoretá_Lago', 'Monte_Caseros', 'Paso_Campamento', 'Paso_de_la_Cruz', 'Paso_de_los_Libres', 'Paso_del_León', 'Paso_del_Remanso', 'Paso_Farías', 'Paso_Potrero', 'Pepirí_Miní', 'Puerto_Concordia', 'Pujol', 'Queguay_Ruta_3', 'Quintana', 'Salto', 'Salto_Grande', 'San_Jaime', 'San_Roquito', 'Santo_Tomé', 'Sarandí_de_Arapey', 'Sequeira', 'Solari', 'Uguay']
# ~ List of C.IV:  ['Baltasar_Brum', 'San_Javier', 'Santa_Ana']
# ~ List of C.V:  ['El_Trompo', 'Federación', 'Guaviyú_de_Arapey', 'Javier_de_Viana', 'Loma_Alta', 'Los_Conquistadores', 'Meneses', 'Miriñay_Medio', 'Miriñay_Río', 'Mocoreta_Medio', 'Mocoretá_Río']
# ~ exit()
