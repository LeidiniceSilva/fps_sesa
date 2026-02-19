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

python3 plot_graph_pdf_sim_obs_sesa_t2m_uv10.py

}
