# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "Dictionary of the inmet automatic station"

"""
Dictionaries of the station

> {station}

    Stations and locations:

    where:

        1 - automatic station name: code, name, latitude, longiitude

    Examples:

         1: ['A899', 'SANTA_VITORIA_DO_PALMAR', -33,74166666, -53,37138888']
"""

inmet = {1	:['A899', 'SANTA_VITORIA_DO_PALMAR',			 -33.74166666,	-53.37138888],
		2	:['A836', 'JAGUARAO',	                         -32.53749999,	-53.37638888],
		3	:['A802', 'RIO_GRANDE',	                         -32.03333333,	-52.10000000],
		4	:['A887', 'CAPAO_DO_LEÃO',	                     -31.80250000,	-52.40722222],
		5	:['A811', 'CANGUÇU',	                         -31.40583333,	-52.70111111],
		6	:['A827', 'BAGE',	                             -31.34777777,	-54.01333333],
		7	:['A878', 'MOSTARDAS',                           -31.24777777,	-50.90555555],
		8	:['A881', 'DOM_PEDRITO',	                     -31.00250000,	-54.61805554],
		9	:['A804', 'SANTANA_DO_LIVRAMENTO',   	         -30.84249999,	-55.61277777],
		10	:['A838', 'CAMAQUA',	                         -30.81055555,	-51.83472221],
		11	:['A812', 'CAÇAPAVA_DO_SUL',	                 -30.54777777,	-53.46749999],
		12	:['A893', 'ENCRUZILHADA_DO_SUL',	             -30.54305555,	-52.52472221],
		13	:['A831', 'QUARAI',	                             -30.36861110,	-56.43722221],
		14	:['A832', 'SãO_GABRIEL',	                     -30.34138888,	-54.31083333],
		15	:['A801', 'PORTO_ALEGRE',	                     -30.05000000,	-51.16666666],
		16	:['A834', 'TRAMANDAI',	                         -30.00972222,	-50.13527777],
		17	:['A886', 'TUPANCIRETA',	                     -29.89388888,	-53.82666666],
		18	:['A813', 'RIO_PARDO',	                         -29.87333332,	-52.38249999],
		19	:['A809', 'URUGUAIANA',	                         -29.84249999,	-57.08249999],
		20	:['A826', 'ALEGRETE',	                         -29.71166666,	-55.52611110],
		21	:['A803', 'SANTA_MARIA',	                     -29.70833333,	-53.70361111],
		22	:['A889', 'SAO_VICENTE_DO_SUL',	                 -29.70194444,	-54.69416666],
		23	:['A884', 'CAMPO_BOM',	                		 -29.67444443,	-51.06416666],
		24	:['A882', 'TEUTONIA',	               			 -29.45027777,	-51.82416666],
		25	:['A879', 'CANELA',	                             -29.36888888,	-50.82749999],
		26	:['A808', 'TORRES',	                             -29.35027777,	-49.73305554],
		27	:['A833', 'SANTIAGO',	                         -29.19138888,	-54.88555555],
		28	:['A840', 'BENTO_GONCALVES',	                 -29.16722221,	-51.53472221],
		29	:['A897', 'CAMBARA_DO_SUL',	                     -29.04888888,	-50.14944444],
		30	:['A867', 'ARARANGUA',	                         -28.93138888,	-49.49805555],
		31	:['A837', 'SOLEDADE',	                         -28.85361111,	-52.55833333],
		32	:['A829', 'SAO_JOSE_DOS_AUSENTES',	             -28.75138888,	-50.05833333],
		33	:['A894', 'SERAFINA_CORREA',	                 -28.70472222,	-51.87055554],
		34	:['A883', 'IBIRUBA',	                         -28.65333333,	-53.11194444],
		35	:['A830', 'SAO_BORJA',	                         -28.64944444,	-56.01555555],
		36	:['A866', 'LAGUNA_FAROL_DE_SANTA_MARTA',	     -28.60416666,	-48.81333333],
		37	:['A853', 'CRUZ_ALTA',	                         -28.60361111,	-53.67361110],
		38	:['A814', 'URUSSANGA',	                         -28.53249999,	-49.31555555],
		39	:['A880', 'VACARIA',	                         -28.51361111,	-50.88277777],
		40	:['A852', 'SÃO_LUIZ_GONZAGA',	                 -28.41694443,	-54.96222222],
		41	:['A815', 'SÃO_JOAQUIM',	                     -28.27555554,	-49.93444444],
		42	:['A839', 'PASSO_FUNDO',	                     -28.22944443,	-52.40388888],
		43	:['A844', 'LAGOA_VERMELHA',	                     -28.22194443,	-51.51222222],
		44	:['A845', 'BOM_JARDIM_DA_SERRA_MORRO_DA_IGREJA', -28.13333333,	-49.48333333],
		45	:['A895', 'CHAPECO',	                         -27.95527777,	-52.63555555],
		46	:['A856', 'PALMEIRA_DAS_MISSOES',	             -27.92166666,	-53.31694443],
		47	:['A810', 'SANTA_ROSA',	                         -27.89305555,	-54.48055554],
		48	:['A805', 'SANTO_AUGUSTO',	                     -27.85416666,	-53.79111110],
		49	:['A865', 'LAGES',	                             -27.80222222,	-50.33555555],
		50	:['A870', 'RANCHO_QUEIMADO',	                 -27.67861110,	-49.04194444],
		51	:['A828', 'ERECHIM',	                         -27.66027777,	-52.30638888],
		52	:['A806', 'FLORIANOPOLIS',	                     -27.60250000,	-48.61999999],
		53	:['A863', 'ITUPORANGA',	                         -27.41833332,	-49.64666666],
		54	:['A854', 'FREDERICO_WESTPHALEN',	             -27.39555555,	-53.41277777],
		55	:['A898', 'CAMPOS_NOVOS',	                     -27.38861110,	-51.21583333],
		56	:['A860', 'CURITIBANOS',	                     -27.28861110,	-50.60416666],
		57	:['A841', 'JOACABA',	                         -27.16944443,	-51.57555554],
		58	:['A868', 'ITAJAI',	                             -26.95083333,	-48.76194444],
		59	:['A858', 'XANXERE',	                         -26.93861110,	-52.41500000],
		60	:['A861', 'RIO_DO_CAMPO',	                     -26.93749999,	-50.14638888],
		61	:['A817', 'INDAIAL',	                         -26.91361111,	-49.26749999],
		62	:['A859', 'CACADOR',	                         -26.81944443,	-50.98555555],
		63	:['A857', 'SAO_MIGUEL_DO_OESTE',	             -26.77666666,	-53.50444444],
		64	:['A876', 'CLEVELANDIA',	                     -26.41722221,	-52.34888888],
		65	:['A816', 'NOVO_HORIZONTE',	                     -26.40694444,	-52.85055555],
		66	:['A875', 'GENERAL_CARNEIRO',	                 -26.39861110,	-51.35444444],
		67	:['A864', 'MAJOR_VIEIRA',	                     -26.39361110,	-50.36333333],
		68	:['A848', 'DIONISIO_CERQUEIRA',	                 -26.28638888,	-53.63277777],
		69	:['A862', 'RIO_NEGRINHO',	                     -26.24861110,	-49.58055554],
		70	:['A851', 'ITAPOA',	                             -26.08111110,	-48.64166666],
		71	:['A874', 'SAO_MATEUS_DO_SUL',	                 -25.83499999,	-50.36861110],
		72	:['A855', 'PLANALTO',	                         -25.72194443,	-53.74805555],
		73	:['A843', 'DOIS_VIZINHOS',	                     -25.69472221,	-53.09444444],
		74	:['A846', 'FOZ_DO_IGUACU',	                     -25.60194444,	-54.48305554],
		75	:['A823', 'INACIO_MARTINS',	                     -25.56583333,	-51.08944444],
		76	:['A873', 'MORRETES',	                         -25.51277777,	-48.80861111],
		77	:['A847', 'ILHA_DO_MEL',	                     -25.49444444,	-48.32583332],
		78	:['A807', 'CURITIBA',	                         -25.43333333,	-49.26666666],
		79	:['B804', 'LARANJEIRAS_DO_SUL',	                 -25.36888888,	-52.39194444],
		80	:['B806', 'COLOMBO',	                         -25.32222221,	-49.15750000],
		81	:['A818', 'IVAI',	                             -25.01333333,	-50.87111110],
		82	:['A746', 'BARRA_DO_TURVO',	                     -24.96305555,	-48.41638888],
		83	:['A819', 'CASTRO',	                             -24.78944444,	-49.99972221],
		84	:['A712', 'IGUAPE',	                             -24.71666666,	-47.55000000],
		85	:['B803', 'CAMPINA_DA_LAGOA',	                 -24.57083332,	-52.80250000],
		86	:['A820', 'MAL._CANDIDO_RONDON',	             -24.53583333,	-54.01972221],
		87	:['A766', 'REGISTRO',	                         -24.53305554,	-47.86416666],
		88	:['A822', 'NOVA_TEBAS',	                         -24.43722221,	-51.96305555],
		89	:['A872', 'VENTANIA',	                         -24.23833333,	-50.24555555],
		90	:['A825', 'GOIOERE',	                         -24.15833333,	-53.03055554],
		91	:['A714', 'ITAPEVA',	                         -23.98138888,	-48.88527777],
		92	:['A751', 'SETE_QUEDAS',	                     -23.96694443,	-55.02416666],
		93	:['A744', 'BRAGANCA_PAULISTA',	                 -23.95194444,	-46.53055554],
		94	:['A715', 'SAO_MIGUEL_ARCANJO',	                 -23.85138888,	-48.16444444],
		95	:['A765', 'BERTIOGA',	                         -23.84499999,	-46.14305555],
		96	:['A767', 'SAO_SEBASTIAO',	                     -23.81083333,	-45.40250000],
		97	:['A871', 'JAPIRA',	                             -23.77305554,	-50.18055554],
		98	:['A771', 'SAO_PAULO_INTERLAGOS',	             -23.72472221,	-46.67805554],
		99	:['S709', 'IGUATEMI',	                         -23.64444444,	-54.57027777],
		100	:['A755', 'BARUERI',	                         -23.52333332,	-46.86916666],
		101	:['A821', 'JOAQUIM_TAVORA',	                     -23.50527777,	-49.94638888],
		102	:['A701', 'SAO_PAULO_MIRANTE',	                 -23.48333333,	-46.61666666],
		103	:['A752', 'ITAQUIRAI',	                         -23.44944444,	-54.18166666],
		104	:['A835', 'MARINGA',	                         -23.40777777,	-51.91666666],
		105	:['A842', 'NOVA_FATIMA',	                     -23.40694444,	-50.57777777],
		106	:['A824', 'ICARAIMA',	                         -23.39027777,	-53.63472221],
		107	:['A869', 'CIDADE_GAUCHA',	                     -23.37583332,	-52.93194443],
		108	:['A713', 'SOROCABA',                            -23.35000000,	-47.66666666],
		109	:['A740', 'SAO_LUIS_DO_PARAITINGA',	             -23.22833332,	-45.41722221],
		110	:['A619', 'PARATI',	                             -23.22333332,	-44.72666666],
		111	:['A725', 'AVARE',	                             -23.09972221,	-48.94555555],
		112	:['A602', 'MARAMBAIA',	                         -23.05000000,	-43.60000000],
		113	:['A728', 'TAUBATE',	                         -23.04194444,	-45.52027777],
		114	:['A750', 'AMAMBAI',	                         -23.00250000,	-55.32944443],
		115	:['A652', 'FORTE_DE_COPACABANA',	             -22.98833333,	-43.19027777],
		116	:['A606', 'ARRAIAL_DO_CABO',	                 -22.98333333,	-42.01666666],
		117	:['A628', 'ANGRA_DOS_REIS',	                     -22.97583332,	-44.30333333],
		118	:['S702', 'ARAL_MOREIRA',	                     -22.95500000,	-55.62611110],
		119	:['A716', 'OURINHOS',	                         -22.94861110,	-49.89416666],
		120	:['A636', 'RIO_DE_JANEIRO_JACAREPAGUA',	         -22.93972221,	-43.40277777],
		121	:['A667', 'SAQUAREMA',	                         -22.87111110,	-42.60916666],
		122	:['A627', 'NITEROI',	                         -22.86749999,	-43.10194444],
		123	:['A509', 'MONTE_VERDE',	                     -22.86138888,	-46.04333333],
		124	:['A621', 'VILA_MILITAR',	                     -22.86083333,	-43.41111111],
		125	:['A749', 'JUTI',	                             -22.85722222,	-54.60555555],
		126	:['A601', 'ECOLOGIA_AGRICOLA',	                 -22.80000000,	-43.68333333],
		127	:['A706', 'CAMPOS_DO_JORDAO',	                 -22.75027777,	-45.60388888],
		128	:['A726', 'PIRACICABA',	                         -22.70277777,	-47.62305554],
		129	:['A769', 'CACHOEIRA_PAULISTA',	                 -22.68888888,	-45.00555555],
		130	:['S706', 'CAARAPO',	                         -22.65694444,	-54.81944443],
		131	:['A626', 'RIO_CLARO',	                         -22.65361111,	-44.04083333],
		132	:['A659', 'SILVA_JARDIM',	                     -22.64583333,	-42.41583333],
		133	:['A849', 'DIAMANTE_DO_NORTE',	                 -22.63388888,	-52.85888888],
		134	:['A603', 'XEREM',	                             -22.58972221,	-43.28222221],
		135	:['S711', 'LAGUNA CARAPA',	                     -22.57527777,	-55.16027777],
		136	:['A703', 'PONTA_PORA',	                         -22.55250000,	-55.71638888],
		137	:['A850', 'PARANAPOEMA',	                     -22.49166666,	-52.13444444],
		138	:['A610', 'PICO_DO_COUTO',	                     -22.48166666,	-43.29138888],
		139	:['A609', 'RESENDE',	                         -22.45000000,	-44.45000000],
		140	:['A618', 'TERESOPOLIS',	                     -22.44888888,	-42.98722221],
		141	:['A739', 'ITAPIRA',	                         -22.41500000,	-46.80527777],
		142	:['A529', 'PASSA_QUATRO',	                     -22.39611110,	-44.96166666],
		143	:['A608', 'MACAE',	                             -22.38333333,	-41.81666667],
		144	:['A635', 'ITATIAIA',	                         -22.37388888,	-44.70305555],
		145	:['A718', 'RANCHARIA',	                         -22.37249999,	-50.97416666],
		146	:['A741', 'BARRA_BONITA',	                     -22.37083332,	-48.55722222],
		147	:['A705', 'BAURU',	                             -22.35805555,	-49.02888888],
		148	:['A611', 'VALENCA',	                         -22.35000000,	-43.70000000],
		149	:['A624', 'NOVA_FRIBURGO',	                     -22.33305554,	-42.67722221],
		150	:['A531', 'MARIA_DA_FE',	                     -22.31444444,	-45.37277777],
		151	:['S708', 'FATIMA_DO_SUL',	                     -22.30861111,	-54.32583332],
		152	:['A709', 'IVINHEMA',	                         -22.30666666,	-53.82722221],
		153	:['A763', 'MARILIA',	                         -22.23527777,	-49.96500000],
		154	:['A721', 'DOURADOS',	                         -22.19388888,	-54.91138888],
		155	:['S701', 'ANGELICA',	                         -22.14805555,	-53.76361111],
		156	:['A707', 'PRESIDENTE_PRUDENTE',	             -22.11999999,	-51.40000000],
		157	:['A757', 'BELA_VISTA',	                         -22.10083333,	-56.53999999],
		158	:['A625', 'TRES_RIOS',	                         -22.09833333,	-43.20861111],
		159	:['S710', 'ITAPORA',	                         -22.09277777,	-54.79888888],
		160	:['S713', 'NOVA_ANDRADINA',	                     -22.07861110,	-53.46583333],
		161	:['A620', 'SAO_TOME',	                         -22.04166666,	-41.05194444],
		162	:['A711', 'SAO_CARLOS',	                         -21.97972221,	-47.88333333],
		163	:['A630', 'SANTA_MARIA_MADALENA',	             -21.95055555,	-42.01055555],
		164	:['A629', 'CARMO',	                             -21.93861110,	-42.60083333],
		165	:['A768', 'TUPA',	                             -21.92722221,	-50.49027777],
		166	:['A530', 'CALDAS',	                             -21.91805554,	-46.38305554],
		167	:['A737', 'IBITINGA',	                         -21.85555555,	-48.66666666],
		168	:['A738', 'CASA_BRANCA',	                     -21.77972221,	-47.07972221],
		169	:['A743', 'RIO_BRILHANTE',	                     -21.77472221,	-54.52805554],
		170	:['A518', 'JUIZ_DE_FORA',	                     -21.76999999,	-43.36416666],
		171	:['A759', 'BATAGUASSU',	                         -21.75138888,	-52.47055554],
		172	:['A607', 'CAMPOS',	                             -21.71666666,	-41.35000000],
		173	:['A723', 'PORTO_MURTINHO',	                     -21.70583333,	-57.88666666],
		174	:['A567', 'MACHADO',	                         -21.68083332,	-45.94444444],
		175	:['A727', 'LINS',	                             -21.66527777,	-49.73416666],
		176	:['A731', 'MARACAJU',	                         -21.60916666,	-55.17777777],
		177	:['A604', 'CAMBUCI',	                         -21.56666666,	-41.95000000],
		178	:['A515', 'VARGINHA',	                         -21.56638888,	-45.40416666],
		179	:['A557', 'CORONEL_PACHECO',	                 -21.54694444,	-43.26111111],
		180	:['A758', 'JARDIM',	                             -21.47833332,	-56.13694444],
		181	:['A770', 'SAO_SIMAO',	                         -21.46111111,	-47.57944443],
		182	:['A762', 'DRACENA',	                         -21.45777777,	-51.55222222],
		183	:['S712', 'NOVA_ALVORADA_DO_SUL',	             -21.45111111,	-54.34194444],
		184	:['A747', 'PRADOPOLIS',	                         -21.33833333,	-48.11388888],
		185	:['A734', 'VALPARAISO',	                         -21.31916666,	-50.93027777],
		186	:['S716', 'SANTA_RITA_DO_PARDO',	             -21.30583333,	-52.82027777],
		187	:['S705', 'BRASILANDIA',	                     -21.29833333,	-52.06888888],
		188	:['S704', 'BONITO',	                             -21.24666666,	-56.45055555],
		189	:['A502', 'BARBACENA',	                         -21.22888888,	-43.76694443],
		190	:['A514', 'SAO_JOAO_DEL_REI',	                 -21.22638888,	-44.97944443],
		191	:['A736', 'ARIRANHA',	                         -21.13305554,	-48.84027777],
		192	:['A517', 'MURIAE',	                             -21.10472222,	-42.37583332],
		193	:['A622', 'PRES._KENNEDY',	                     -21.10111111,	-41.03888888],
		194	:['A735', 'JOSE_BONIFACIO',	                     -21.08555555,	-49.92055554],
		195	:['A754', 'SIDROLANDIA',	                     -20.98166666,	-54.97194443],
		196	:['A764', 'BEBDOURO',	                         -20.94916666,	-48.48972221],
		197	:['A561', 'SAO_SEBASTIAO_DO_PARAISO',	         -20.91000000,	-47.11416666],
		198	:['A704', 'TRES_LAGOAS',	                     -20.78999999,	-51.71222222],
		199	:['A510', 'VICOSA',	                             -20.76277777,	-42.86388888],
		200	:['A617', 'ALEGRE',	                             -20.75055555,	-41.48888888],
		201	:['A516', 'PASSOS',	                             -20.74499999,	-46.63388888],
		202	:['A570', 'OLIVEIRA',	                         -20.71472222,	-44.86444444],
		203	:['A615', 'ALFREDO_CHAVES',	                     -20.63638888,	-40.74138888],
		204	:['A708', 'FRANCA',	                             -20.57999999,	-47.37999999],
		205	:['A748', 'BARRETOS',	                         -20.55888888,	-48.54472221],
		206	:['A513', 'OURO_BRANCO',	                     -20.54583333,	-43.76944443],
		207	:['A719', 'AQUIDAUANA',	                         -20.47555554,	-55.78388888],
		208	:['A634', 'VILA_VELHA',	                         -20.46694443,	-40.40416666],
		209	:['S715', 'RIBAS_DO_RIO_PARDO',	                 -20.46666666,	-53.76305555],
		210	:['A524', 'FORMIGA',	                         -20.45472222,	-45.45388888],
		211	:['A702', 'CAMPO_GRANDE',	                     -20.45000000,	-54.60000000],
		212	:['A756', 'AGUA_CLARA',	                         -20.44444444,	-52.87555554],
		213	:['A729', 'VOTUPORANGA',	                     -20.40305555,	-49.96583333],
		214	:['A722', 'MIRANDA',	                         -20.38333333,	-56.41666666],
		215	:['A753', 'ITUVERAVA',	                         -20.35944444,	-47.77499999],
		216	:['S717', 'SELVIRIA',	                         -20.35138888,	-51.43027777],
		217	:['A612', 'VITORIA',	                         -20.26666666,	-40.30000000],
		218	:['A556', 'MANHUACU',	                         -20.26361111,	-42.18277777],
		219	:['A633', 'VENDA_NOVA_DO_IMIGRANTE',	         -20.25222222,	-41.18999999],
		220	:['A564', 'DIVINOPOLIS',	                     -20.17305554,	-44.87472221],
		221	:['A733', 'JALES',	                             -20.16500000,	-50.59499999],
		222	:['A657', 'AFONSO_CLAUDIO',	                     -20.10416666,	-41.10694444],
		223	:['A555', 'IBIRITE_ROLA_MOCA',	                 -20.03111110,	-44.01111111],
		224	:['A565', 'BAMBUI',	                             -20.03111110,	-46.00861111],
		225	:['A613', 'SANTA_TERESA',	                     -19.98833333,	-40.57944443],
		226	:['A520', 'CONCEICAO_DAS_ALAGOAS',	             -19.98583333,	-48.15138888],
		227	:['F501', 'BELO_HORIZONTE_CERCADINHO',	         -19.97999999,	-43.95861111],
		228	:['S703', 'BANDEIRANTES',	                     -19.94555555,	-54.36861110],
		229	:['A535', 'FLORESTAL',	                         -19.88527778,	-44.41694444],
		230	:['A521', 'PAMPULHA',	                         -19.88416666,	-43.96944443],
		231	:['A525', 'SACRAMENTO',	                         -19.87527777,	-47.43416666],
		232	:['A554', 'CARATINGA',	                         -19.73583333,	-42.15361111],
		233	:['A568', 'UBERABA',	                         -19.71000000,	-47.96166666],
		234	:['A710', 'PARANAIBA',	                         -19.69527777,	-51.18138888],
		235	:['A505', 'ARAXA',	                             -19.60555555,	-46.94944444],
		236	:['S707', 'CAMAPUA',	                         -19.58749999,	-54.02999999],
		237	:['A511', 'TIMOTEO',	                         -19.57361110,	-42.62222221],
		238	:['A519', 'CAMPINA_VERDE',	                     -19.53916660,	-49.51805554],
		239	:['A534', 'AIMORES',	                         -19.53277770,	-41.09083333],
		240	:['A536', 'DORES_DO_INDAIA',	                 -19.48166666,	-45.59388888],
		241	:['A569', 'SETE_LAGOAS',	                     -19.45500000,	-44.17333332],
		242	:['A732', 'SAO_GABRIEL_DO_OESTE',	             -19.42027777,	-54.55305555],
		243	:['A632', 'MARILANDIA',	                         -19.40722222,	-40.53972221],
		244	:['A614', 'LINHARES',	                         -19.35694444,	-40.06861110],
		245	:['A560', 'POMPEU',	                             -19.23249999,	-44.96416666],
		246	:['A742', 'CASSILANDIA',	                     -19.12249999,	-51.72083332],
		247	:['A724', 'CORUMBA',	                         -18.99666666,	-57.63749999],
		248	:['A523', 'PATROCINIO',	                         -18.99638888,	-46.98583333],
		249	:['A717', 'NHUMIRIM',	                         -18.98888888,	-56.62305554],
		250	:['A011', 'SAO_SIMAO',	                         -18.96666666,	-50.61666666],
		251	:['A512', 'ITUIUTABA',	                         -18.95277777,	-49.52499999],
		252	:['A507', 'UBERLANDIA',	                         -18.91722221,	-48.25555555],
		253	:['A532', 'GOVERNADOR_VALADARES',	             -18.83027777,	-41.97666666],
		254	:['A760', 'COSTA_RICA',	                         -18.82277777,	-53.26722221],
		255	:['A730', 'CHAPADAO_DO_SUL',	                 -18.80222222,	-52.60277777],
		256	:['A533', 'GUANHAES',	                         -18.78666666,	-42.94305555],
		257	:['A540', 'MANTENA',	                         -18.78027777,	-40.98583333],
		258	:['A538', 'CURVELO',	                         -18.76444444,	-44.45361111],
		259	:['A616', 'SAO_MATEUS',	                         -18.71388888,	-39.84833333],
		260	:['A623', 'NOVA_VENECIA',	                     -18.69527777,	-40.40750000],
		261	:['A562', 'PATOS_DE_MINAS',	                     -18.52055554,	-46.44055555],
		262	:['A720', 'COXIM',	                             -18.51222222,	-54.73583333],
		263	:['A035', 'ITUMBIARA',	                         -18.40972222,	-49.19194444],
		264	:['A631', 'ECOPORANGA',	                         -18.29166666,	-40.73638888],
		265	:['A537', 'DIAMANTINA',	                         -18.23111110,	-43.64805555],
		266	:['A528', 'TRES_MARIAS',	                     -18.20055555,	-45.45972222],
		267	:['A034', 'CATALAO',	                         -18.15777777,	-47.92638888],
		268	:['S714', 'PEDRO_GOMES',	                     -18.07277777,	-54.54888888],
		269	:['A422', 'ABROLHOS',	                         -17.96416666,	-38.69611110],
		270	:['A016', 'JATAI',	                             -17.92388888,	-51.71777777],
		271	:['A527', 'TEOFILO_OTONI',	                     -17.90000000,	-41.51666666],
		272	:['A761', 'SONORA',	                             -17.89638888,	-54.45361111],
		273	:['A934', 'ALTO_TAQUARI',	                     -17.81555555,	-53.28749999],
		274	:['A522', 'SERRA_DOS_AIMORES',	                 -17.79861110,	-40.24972221],
		275	:['A025', 'RIO_VERDE',	                         -17.78583333,	-50.98138888],
		276	:['A553', 'JOAO_PINHEIRO',	                     -17.78472221,	-46.11916666],
		277	:['A405', 'CARAVELAS',	                         -17.73333333,	-39.26666666],
		278	:['A003', 'MORRINHOS',	                         -17.71666667,	-49.10000000],
		279	:['A541', 'CAPELINHA',	                         -17.70527777,	-42.38916666],
		280	:['A546', 'GUARDA-MOR',	                         -17.56166666,	-47.19944444],
		281	:['A026', 'MINEIROS',	                         -17.45472222,	-52.60111111],
		282	:['A909', 'ALTO_ARAGUAIA',	                     -17.33944444,	-53.22416666],
		283	:['A029', 'EDEIA',	                             -17.33694444,	-49.91472222],
		284	:['A033', 'PIRES_DO_RIO',	                     -17.30416666,	-48.28416666],
		285	:['A933', 'ITIQUIRA',	                         -17.29722221,  -54.83722221],
		286	:['A545', 'PIRAPORA',	                         -17.25777777,	-44.83555555],
		287	:['A571', 'PARACATU',	                         -17.24416666,	-46.88166666],
		288	:['A455', 'ITAMARAJU',	                         -17.00722222,	-39.55138888]}
		
