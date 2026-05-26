#!/bin/bash

#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J Plot
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Feb 16, 2026'
#__description__ = 'Extract timeseries from CP simulations'
 
{
set -eo pipefail

#python3 plot_graph_pdf_sim_obs_sesa_pr.py
python3 plot_graph_pdf_sim_obs_sesa_t2m_uv10.py --var tas
#python3 plot_graph_pdf_sim_obs_sesa_t2m_uv10.py --var sfcWind


}
