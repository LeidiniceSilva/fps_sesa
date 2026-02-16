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
#__description__ = 'Posprocessing the OBS with CDO'

{
set -eo pipefail

CDO(){
  cdo -O -L -f nc4 -z zip $@
}

EXP="CSAM-4i_ERA5"
DT="2018060100-2021053123"
VAR_LIST="tp t2m ws10"

DIR="/home/mda_silv/users/FPS_SESA/database/obs/era5"

echo
cd ${DIR}
echo ${DIR}

echo
echo "--------------- INIT POSPROCESSING MODEL ----------------"

for VAR in ${VAR_LIST[@]}; do

    if [ ${VAR} = 'tp'  ]; then
    P=99.9
    elif [ ${VAR} = 't2m'  ]; then
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
    done
	
    CDO mergetime p${P}_${VAR}_${EXP}_*_${DT}.nc p${P}_${VAR}_diurnal-cycle_${EXP}_${DT}.nc
    CDO mergetime freq_${VAR}_${EXP}_*_${DT}.nc freq_${VAR}_diurnal-cycle_${EXP}_${DT}.nc
    CDO mergetime int_${VAR}_${EXP}_*_${DT}.nc int_${VAR}_diurnal-cycle_${EXP}_${DT}.nc
    rm -f p${P}_${VAR}_${EXP}_*_${DT}.nc freq_${VAR}_${EXP}_*_${DT}.nc int_${VAR}_${EXP}_*_${DT}.nc
  
done

echo
echo "--------------- THE END POSPROCESSING MODEL ----------------"

}
