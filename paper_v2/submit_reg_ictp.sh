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

python3 test_reg_ictp.py --var sfcWind --mdl reg_ictp --inst INMET

}
