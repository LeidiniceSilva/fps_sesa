#!/bin/bash

#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH -J Plot
#SBATCH -p testing
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mda_silv@ictp.it

python3 plot_maps_bias_sim_obs_sesa_pr.py
python3 plot_maps_bias_sim_obs_sesa_t2m_uv10.py
python3 plot_maps_clim_sim_obs_sesa_pr.py
python3 plot_maps_clim_sim_obs_sesa_t2m_uv10.py
python3 plot_maps_kge_sim_obs_sesa_pr.py
python3 plot_maps_kge_sim_obs_sesa_t2m_uv10.py
