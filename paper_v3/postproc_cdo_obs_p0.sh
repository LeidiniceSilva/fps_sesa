#!/bin/bash

#SBATCH -A CMPNS_ictpclim
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH --ntasks-per-node=112
#SBATCH -t 1-00:00:00
#SBATCH -J Postproc
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

YR="2018-2021"
EXP="CSAM-4i_ERA5"
VAR_LIST="t2m"

DIR_IN="/leonardo/home/userexternal/mdasilva/leonardo_work/OBS/ERA5"
DIR_OUT="/leonardo/home/userexternal/mdasilva/leonardo_work/SAM-3km/postproc/fps_sesa"

echo
cd ${DIR_OUT}
echo ${DIR_OUT}

echo
echo "--------------- INIT POSPROCESSING MODEL ----------------"

for VAR in ${VAR_LIST[@]}; do

    if [ ${VAR} = 'tp' ]; then
    CDO seldate,2018-06-01,2021-05-31 ${DIR_IN}/${VAR}_ERA5_1hr_2018-2021.nc ${VAR}_ERA5_1hr_2018-2021.nc
    CDO mulc,1000 ${VAR}_ERA5_1hr_2018-2021.nc ${VAR}_ERA5_1hr_2018060100-2021053100.nc
    CDO remapdis,/leonardo/home/userexternal/mdasilva/grid.txt ${VAR}_ERA5_1hr_2018060100-2021053100.nc ${VAR}_${EXP}_1hr_2018060100-2021053123.nc

    elif [ ${VAR} = 't2m' ]; then
    CDO mergetime ${DIR_IN}/${VAR}_ERA5_1hr_20*.nc ${VAR}_ERA5_1hr_2018010100-2021123100.nc
    CDO seldate,2018-06-01,2021-05-31 ${VAR}_ERA5_1hr_2018010100-2021123100.nc ${VAR}_ERA5_1hr_2018-2021.nc
    CDO subc,273.15 ${VAR}_ERA5_1hr_2018-2021.nc ${VAR}_ERA5_1hr_2018060100-2021053100.nc
    CDO remapdis,/leonardo/home/userexternal/mdasilva/grid.txt ${VAR}_ERA5_1hr_2018060100-2021053100.nc ${VAR}_${EXP}_1hr_2018060100-2021053123.nc
    
    else
    CDO expr,'ws10=sqrt(u10*u10+v10*v10)' -merge ${DIR_IN}/u10_ERA5_1hr_2018-2021.nc ${DIR_IN}/v10_ERA5_1hr_2018-2021.nc uv10_ERA5_1hr_2018-2021.nc
    CDO seldate,2018-06-01,2021-05-31 uv10_ERA5_1hr_2018-2021.nc uv10_ERA5_1hr_2018060100-2021053100.nc
    CDO remapdis,/leonardo/home/userexternal/mdasilva/grid.txt uv10_ERA5_1hr_2018060100-2021053100.nc uv10_${EXP}_1hr_2018060100-2021053123.nc
    fi
  
done

echo
echo "--------------- THE END POSPROCESSING MODEL ----------------"

}
