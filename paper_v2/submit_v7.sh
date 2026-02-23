#!/bin/bash

#SBATCH -N 2
#SBATCH -t 24:00:00
#SBATCH -J FPS-SESA
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Feb 16, 2026'
#__description__ = 'Extract timeseries from CP simulations'
 
{
set -eo pipefail

python3 extract_ts_inmet_smn_i_sim.py --var pr --mdl wrf_ucan --inst SMN
python3 extract_ts_inmet_smn_i_sim.py --var pr --mdl wrf_ucan --inst INMET

python3 extract_ts_inmet_smn_i_sim.py --var tas --mdl wrf_ucan --inst SMN 
python3 extract_ts_inmet_smn_i_sim.py --var tas --mdl wrf_ucan --inst INMET

python3 extract_ts_inmet_smn_i_sim.py --var sfcWind --mdl wrf_ucan --inst SMN
python3 extract_ts_inmet_smn_i_sim.py --var sfcWind --mdl wrf_ucan --inst INMET

}
