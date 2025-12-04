#!/bin/bash

#SBATCH -A CMPNS_ictpclim
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH --ntasks-per-node=112
#SBATCH -t 1-00:00:00
#SBATCH -J sigma2p
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Dec 03, 2025'
#__description__ = 'Posprocessing the RCM output with CDO'
 
{
set -eo pipefail
CDO(){
  cdo -O -L -f nc4 -z zip $@
}

VAR_LIST="sfcWind"
YR="2018-2021"
IYR=$( echo $YR | cut -d- -f1 )
FYR=$( echo $YR | cut -d- -f2 )

DIR_IN="/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km/rcm/output"
DIR_OUT="/leonardo/home/userexternal/mdasilva"
BIN="/leonardo/home/userexternal/mdasilva/RegCM/bin"

echo
cd ${DIR_OUT}
echo ${DIR_OUT}

for VAR in ${VAR_LIST[@]}; do
    for YEAR in `seq -w ${IYR} ${FYR}`; do
        for MON in `seq -w 01 12`; do
	    if [ ${VAR} = 'sfcWind'  ]
	    then
            CDO selname,${VAR} ${DIR_IN}/SAM-3km_SRF.${YEAR}${MON}0100.nc ${VAR}_SAM-3km_${YEAR}${MON}0100.nc
	    else 
            CDO selname,${VAR} ${DIR_IN}/SAM-3km_STS.${YEAR}${MON}0100.nc ${VAR}_SAM-3km_${YEAR}${MON}0100.nc
	    fi
        done
    done

    CDO mergetime ${VAR}_SAM-3km_*0100.nc ${VAR}_SAM-3km_${YR}.nc
    CDO seldate,2018-06-01,2021-05-31 ${VAR}_SAM-3km_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_${YR}.nc
    
    if [ ${VAR} = 'pr'  ]
    then
        CDO -b f32 mulc,86400 ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc
	CDO monmean ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_${YR}.nc
	${BIN}/./regrid ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc -35.70235,-11.25009,0.0352 -78.66277,-35.48362,0.0352 bil
	${BIN}/./regrid ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_${YR}.nc -35.70235,-11.25009,0.0352 -78.66277,-35.48362,0.0352 bil
    elif [ ${VAR} = 'tas'  ]
    then
        CDO -b f32 subc,273.15 ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc
	CDO monmean ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_${YR}.nc
	${BIN}/./regrid ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc -35.70235,-11.25009,0.0352 -78.66277,-35.48362,0.0352 bil
	${BIN}/./regrid ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_${YR}.nc -35.70235,-11.25009,0.0352 -78.66277,-35.48362,0.0352 bil
    else
        CDO daymean ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc
	CDO monmean ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_${YR}.nc
	${BIN}/./regrid ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_day_${YR}.nc -35.70235,-11.25009,0.0352 -78.66277,-35.48362,0.0352 bil
	${BIN}/./regrid ${VAR}_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_${YR}.nc -35.70235,-11.25009,0.0352 -78.66277,-35.48362,0.0352 bil
    fi
done


}
