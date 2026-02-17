#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J FPS-SESA
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Feb 16, 2026'
#__description__ = 'Posprocessing the RegCM5 with CDO'
 
{
set -eo pipefail

CDO(){
  cdo -O -L -f nc4 -z zip $@
}

VAR_LIST="pr tas sfcWind"
SEASON_LIST="DJF MAM JJA SON"
DATASET_LIST="wrf_ucan"
#DATASET_LIST="reg_ictp reg_ictp_pbl1 reg_ictp_pbl2 reg_usp wrf_ncar wrf_ucan"

echo
echo "--------------- INIT POSPROCESSING MODEL ----------------"

for DATASET in ${DATASET_LIST[@]}; do

    DIR="/home/mda_silv/users/FPS_SESA/database/rcm/${DATASET}"
    echo
    cd ${DIR}
    echo ${DIR}
    
    if [ ${DATASET} = 'reg_ictp'  ]; then
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_v1" 
    elif [ ${DATASET} = 'reg_ictp_pbl1'  ]; then
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0" 
    elif [ ${DATASET} = 'reg_ictp_pbl2'  ]; then
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0"
    elif [ ${DATASET} = 'reg_usp'  ]; then
    EXP="pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_USP-RegCM471_v2"
    elif [ ${DATASET} = 'wrf_ncar'  ]; then
    EXP="CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1"
    else
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1"
    fi
    
    for VAR in ${VAR_LIST[@]}; do
    
	CDO mergetime ${VAR}_${EXP}_1hr_20*.nc ${VAR}_${EXP}_1hr_2018-2021.nc
	CDO seldate,2018-06-01,2021-05-31 ${VAR}_${EXP}_1hr_2018-2021.nc ${VAR}_${EXP}_1hr_201806-202105.nc
	
	if [ ${VAR} = 'pr'  ]; then
	CDO -b f32 mulc,86400 ${VAR}_${EXP}_1hr_201806-202105.nc ${VAR}_${EXP}_1hr_2018060100-2021053123.nc
	CDO daysum ${VAR}_${EXP}_1hr_2018060100-2021053123.nc ${VAR}_${EXP}_day_2018060100-2021053123.nc
	CDO monmean ${VAR}_${EXP}_day_2018060100-2021053123.nc ${VAR}_${EXP}_mon_2018060100-2021053123.nc

	elif [ ${VAR} = 'tas'  ]; then
	CDO -b f32 subc,273.15 ${VAR}_${EXP}_1hr_201806-202105.nc ${VAR}_${EXP}_1hr_2018060100-2021053123.nc
	CDO daymean ${VAR}_${EXP}_1hr_2018060100-2021053123.nc ${VAR}_${EXP}_day_2018060100-2021053123.nc
	CDO monmean ${VAR}_${EXP}_day_2018060100-2021053123.nc ${VAR}_${EXP}_mon_2018060100-2021053123.nc
	
	else
	cp ${VAR}_${EXP}_1hr_201806-202105.nc ${VAR}_${EXP}_1hr_2018060100-2021053123.nc
	CDO daymean ${VAR}_${EXP}_1hr_2018060100-2021053123.nc ${VAR}_${EXP}_day_2018060100-2021053123.nc
	CDO monmean ${VAR}_${EXP}_day_2018060100-2021053123.nc ${VAR}_${EXP}_mon_2018060100-2021053123.nc
	fi
	
	for SEASON in ${SEASON_LIST[@]}; do
	    CDO selseas,${SEASON} ${VAR}_${EXP}_1hr_2018060100-2021053123.nc ${VAR}_${EXP}_${SEASON}_2018060100-2021053123.nc
	done
    done
done

echo
echo "--------------- THE END POSPROCESSING MODEL ----------------"

}
