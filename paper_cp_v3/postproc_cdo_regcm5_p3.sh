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
DT="2018060100-2021053123"
DATASET_LIST="reg_ictp reg_ictp_pbl1 reg_ictp_pbl2 ictp_usp wrf_ncar wrf_ucan"

echo
echo "--------------- INIT POSPROCESSING MODEL ----------------"

for DATASET in ${DATASET_LIST[@]}; do

    DIR="/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/${DATASET}"
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
    
	if [ ${VAR} = 'pr'  ]; then
	P=99
	elif [ ${VAR} = 'tas'  ]; then
	P=97.5
	else
	P=95
	fi
	    
    	for HR in `seq -w 00 23`; do
	    CDO selhour,${HR} ${VAR}_${EXP}_1hr_${DT}.nc ${VAR}_${EXP}_${HR}_${DT}.nc
	    CDO timpctl,${P} ${VAR}_${EXP}_${HR}_${DT}.nc -timmin ${VAR}_${EXP}_${HR}_${DT}.nc -timmax ${VAR}_${EXP}_${HR}_${DT}.nc p${P}_${VAR}_${EXP}_${HR}_${DT}.nc
	    CDO ifthen -ge ${VAR}_${EXP}_${HR}_${DT}.nc p${P}_${VAR}_${EXP}_${HR}_${DT}.nc ext_${VAR}_${EXP}_${HR}_${DT}.nc
	    CDO div -timcount ext_${VAR}_${EXP}_${HR}_${DT}.nc -timcount ${VAR}_${EXP}_${HR}_${DT}.nc freq_${VAR}_${EXP}_${HR}_${DT}.nc
	    CDO timmean ext_${VAR}_${EXP}_${HR}_${DT}.nc int_${VAR}_${EXP}_${HR}_${DT}.nc
	    rm -f ${VAR}_${EXP}_${HR}_${DT}.nc ext_${VAR}_${EXP}_${HR}_${DT}.nc
	done
	
        CDO mergetime p${P}_${VAR}_${EXP}_*_${DT}.nc p${P}_${VAR}_diurnal-cycle_${EXP}_${DT}.nc
        CDO mergetime freq_${VAR}_${EXP}_*_${DT}.nc freq_${VAR}_diurnal-cycle_${EXP}_${DT}.nc
        CDO mergetime int_${VAR}_${EXP}_*_${DT}.nc int_${VAR}_diurnal-cycle_${EXP}_${DT}.nc
	rm -f p${P}_${VAR}_${EXP}_*_${DT}.nc freq_${VAR}_${EXP}_*_${DT}.nc int_${VAR}_${EXP}_*_${DT}.nc
    
    done
done

echo
echo "--------------- THE END POSPROCESSING MODEL ----------------"

}
