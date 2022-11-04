# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/31/2022"
__description__ = "Dictionary of the smn automatic station"

"""
Dictionaries of the automatic weather station data

> {station}

    Name, latitude and longitude:

    where:

        1 - automatic station code.

    Examples:

         1: ['Aguapey',	-28,416682,	-56,542211]
         
"""

urug_smn = {1	:['Aguapey',	           -28.416682,	-56.542211],
			2	:['Alto_Arapey',	       -27.416682,	-56.224017],
			3	:['Alto_Cuareim',	       -26.416682,	-55.866233],
			4	:['Alvear',	               -25.416682,	-56.555522],
			5	:['Arapey_Chico',	       -24.416682,	-57.455900],
			6	:['Arapey_Grande',	       -23.416682,	-57.485838],
			7	:['Arapey_Grande_Ruta_4',  -22.416682,	-57.097101],
			8	:['Arerunguá',	           -21.416682,	-56.618974],
			9	:['Artigas',	           -20.416682,	-56.456119],
			10	:['Baibene',	           -19.416682,	-58.165178],
			11	:['Baltasar_Brum',	       -18.416682,	-57.325281],
			12	:['Belén',	               -30.787252,	-57.783090],
			13	:['Bernabé_Rivera',	       -30.299302,	-56.972142],
			14	:['Bonpland',	           -29.819517,	-57.429814],
			15	:['Catalán',	           -30.751196,	-56.214658],
			16	:['Catalán_Grande_(rio)',  -30.637259,	-56.326675],
			17	:['Cazadores_Correntinos', -29.983490,	-58.295964],
			18	:['Cerro_Amarillo',	       -30.624270,	-56.657495],
			19	:['Cerro_Chato',	       -31.337833,	-56.629005],
			20	:['Chajarí',	           -30.776241,	-58.141889],
			21	:['Colonia_Lavalleja',	   -31.092420,	-57.021887],
			22	:['Colonia_Libertad',	   -30.014261,	-57.859450],
			23	:['Colonia_Palma',	       -30.580333,	-57.680318],
			24	:['Cuareim_Río',	       -30.586390,	-56.230152],
			25	:['Cuaró',	               -30.613837,	-56.905948],
			26	:['Cuchilla_de_Salto',	   -31.444407,	-57.421509],
			27	:['Curuzú_Cuatiá',	       -29.787069,	-58.070281],
			28	:['Dayman_Ruta 3',	       -31.459000,	-57.903000],
			29	:['Diego_Lamas',	       -30.755734,	-57.055016],
			30	:['El_Soberbio',	       -27.298582,	-54.193437],
			31	:['El_Trompo',	           -30.971026,	-58.247317],
			32	:['Espinillar',            -30.971187,	-57.880716],
			33	:['Federación',	           -30.992189,	-57.913847],
			34	:['Garruchos',	           -28.177515,	-55.642269],
			35	:['Guaviyú_de_Arapey',	   -31.037601,	-56.635879],
			36	:['Itapebí',               -31.286127,	-57.703535],
			37	:['Javier_de_Viana',	   -30.432802,	-56.787211],
			38	:['Loma_Alta',	           -29.056630,	-57.087690],
			39	:['Los_Conquistadores',	   -30.590296,	-58.468999],
			40	:['Meneses',	           -30.865546,	-56.411658],
			41	:['Miriñay_Medio',	       -29.660671,	-57.743836],
			42	:['Miriñay_Río',	       -29.845637,	-57.674446],
			43	:['Mocoretá_Lago',	       -30.690307,	-57.820722],
			44	:['Mocoreta_Medio',	       -30.214445,	-57.966907],
			45	:['Mocoretá_Río',	       -30.626772,	-57.982836],
			46	:['Monte_Caseros',	       -30.249296,	-57.622519],
			47	:['Paso_Campamento',	   -30.787214,	-56.786308],
			48	:['Paso_de_la_Cruz',	   -30.289000,	-57.291694],
			49	:['Paso_de_los_Libres',	   -29.721745,	-57.082856],
			50	:['Paso_del_León',         -30.205281,	-57.074737],
			51	:['Paso_del_Remanso',      -30.223953,	-56.667540],
			52	:['Paso_Farías',           -30.476397,	-57.124992],
			53	:['Paso_Potrero',          -31.444635,	-56.840210],
			54	:['Pepirí_Miní',           -27.153304,	-53.933283],
			55	:['Puerto_Concordia',      -31.404212,	-58.002561],
			56	:['Pujol',                 -30.390258,	-57.880984],
			57	:['Puntas_de_Valentín',    -31,493837,	-57.120524],
			58	:['Queguay_Ruta_3',        -32.135000,	-57.939000],
			59	:['Quintana',              -31.350331,	-56.392531],
			60	:['Salto',                 -31.389780,	-57.977472],
			61	:['Salto_Grande',          -31.272706,	-57.916594],
			62	:['San_Jaime',             -30.337372,	-58.317145],
			63	:['San_Javier',            -27.869343,	-55.129216],
			64	:['San_Roquito',           -29.296744,	-57.565952],
			65	:['Santa_Ana',             -30.914680,	-57.927730],
			66	:['Santo_Tomé',            -28.545073,	-56.028995],
			67	:['Sarandí_de_Arapey',     -30.988557,	-56.215761],
			68	:['Sequeira',              -31.004877,	-56.879045],
			69	:['Solari',	               -29.377358,	-58.193534],
			70	:['Tomás_Gomensoro',       -30.431886,	-57.440628],
			71	:['Uguay',                 -28.691829,	-57.480465],
			72	:['Yapeyú',                -29.470444,	-56.811670]}

