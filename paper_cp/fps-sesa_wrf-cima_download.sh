#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J Download 
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 10, 2026'
#__description__ = 'Download WRF-CIMA'

#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2018123123.nc
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2019010100-2019123123.nc      
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2020010100-2020123123.nc   
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2021010100-2021053123.nc 
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/sfcWind_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2018123123.nc 
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/sfcWind_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2019010100-2019123123.nc           
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/sfcWind_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2020010100-2020123123.nc
#wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/sfcWind_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2021010100-2021053123.nc
wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/tas_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2018123123.nc
wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/tas_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2019010100-2019123123.nc
wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/tas_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2020010100-2020123123.nc
wget --no-check-certificate https://meteo.unican.es/work/josipa/FPS_SESA/CIMA/tas_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2021010100-2021053123.nc
