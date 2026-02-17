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

VAR_LIST="sfcWind"
SEASON_LIST="DJF MAM JJA SON"
DATASET_LIST="wrf_ncar"

echo
echo "--------------- INIT POSPROCESSING MODEL ----------------"

for DATASET in ${DATASET_LIST[@]}; do

    DIR="/home/mda_silv/users/FPS_SESA/database/rcm/${DATASET}"
    echo
    cd ${DIR}
    echo ${DIR}
    
    if [ ${DATASET} = 'reg_ictp'  ]; then
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_v1" 
    elif [ ${DATASET} = 'reg_ictp'  ]; then
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0" 
    elif [ ${DATASET} = 'reg_ictp'  ]; then
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0"
    elif [ ${DATASET} = 'reg_usp'  ]; then
    EXP="pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_USP-RegCM471_v2"
    elif [ ${DATASET} = 'wrf_ncar'  ]; then
    EXP="CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1"
    else
    EXP="CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1"
    fi
    
    for VAR in ${VAR_LIST[@]}; do
    
	#CDO expr,'sfcWind=sqrt(uas*uas+vas*vas)' -merge uas_${EXP}_1hr_2018010100-2018123123.nc vas_${EXP}_1hr_2018010100-2018123123.nc ${VAR}_${EXP}_1hr_2018010100-2018123123.nc
	CDO expr,'sfcWind=sqrt(uas*uas+vas*vas)' -merge uas_${EXP}_1hr_2019010100-2019123123.nc vas_${EXP}_1hr_2019010100-2019123123.nc ${VAR}_${EXP}_1hr_2019010100-2019123123.nc
	CDO expr,'sfcWind=sqrt(uas*uas+vas*vas)' -merge uas_${EXP}_1hr_2020010100-2020123123.nc vas_${EXP}_1hr_2020010100-2020123123.nc ${VAR}_${EXP}_1hr_2020010100-2020123123.nc
	CDO expr,'sfcWind=sqrt(uas*uas+vas*vas)' -merge uas_${EXP}_1hr_2021010100-2021123123.nc vas_${EXP}_1hr_2021010100-2021123123.nc ${VAR}_${EXP}_1hr_2021010100-2021123123.nc

    done
done

echo
echo "--------------- THE END POSPROCESSING MODEL ----------------"

}
