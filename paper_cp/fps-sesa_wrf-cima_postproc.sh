#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J Postproc 
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 10, 2026'
#__description__ = 'Download WRF-CIMA'

{
set -eo pipefail

CDO(){
  cdo -O -L -f nc4 -z zip $@
}

DIR_IN="/home/mda_silv/users/FPS_SESA/database/rcm/wrf_cima"
DIR_OUT="/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_cima"

echo
cd ${DIR_OUT}
echo ${DIR_OUT}

CDO settaxis,2018-06-01,00:00:00,1hour ${DIR_IN}/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2018123123.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018.nc
CDO settaxis,2019-01-01,00:00:00,1hour ${DIR_IN}/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2019010100-2019123123.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2019.nc  
CDO settaxis,2020-01-01,00:00:00,1hour ${DIR_IN}/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2020010100-2020123123.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2020.nc 
CDO settaxis,2021-01-01,00:00:00,1hour ${DIR_IN}/pr_CSAM-4i_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2021010100-2021053123.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2021.nc
CDO mergetime pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_20*.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018-2021.nc
CDO mulc,3600 pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018-2021.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2021053123.nc  
CDO daysum pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2021053123.nc pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_day_2018060100-2021053123.nc

CDO mergetime ${DIR_IN}/tas_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_20*.nc tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018-2021.nc
CDO subc,273.15 tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018-2021.nc tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2021053123.nc  
CDO daymean tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2021053123.nc tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_day_2018060100-2021053123.nc

CDO mergetime ${DIR_IN}/sfcWind_CSAM-4_ECWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_20*.nc sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2021053123.nc
CDO daymean sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_1hr_2018060100-2021053123.nc sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_CIMA-WRF433_v1_daymean_2018060100-2021053123.nc

}

