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
SEASON_LIST="DJF MAM JJA SON"
VAR_LIST="tp t2m ws10"

DIR="/home/mda_silv/users/FPS_SESA/database/obs/era5"

echo
cd ${DIR}
echo ${DIR}

echo
echo "--------------- INIT POSPROCESSING MODEL ----------------"

for VAR in ${VAR_LIST[@]}; do

    if [ ${VAR} = 'tp' ]; then
    CDO daysum ${VAR}_${EXP}_1hr_${DT}.nc ${VAR}_${EXP}_day_${DT}.nc

    else
    CDO daymean ${VAR}_${EXP}_1hr_${DT}.nc ${VAR}_${EXP}_day_${DT}.nc
    fi
    CDO monmean ${VAR}_${EXP}_day_${DT}.nc ${VAR}_${EXP}_mon_${DT}.nc

    for SEASON in ${SEASON_LIST[@]}; do
	CDO selseas,${SEASON} ${VAR}_${EXP}_1hr_${DT}.nc ${VAR}_${EXP}_${SEASON}_${DT}.nc
    done
	
done

echo
echo "--------------- THE END POSPROCESSING MODEL ----------------"

}
