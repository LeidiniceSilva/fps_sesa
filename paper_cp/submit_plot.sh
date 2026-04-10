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

python3 plot_maps_bias_sim_obs_sesa_pr_mmd.py --var pr
python3 plot_maps_bias_sim_obs_sesa_pr_mmh.py --var pr
python3 plot_maps_bias_sim_obs_sesa_pr_seas.py --var pr --seas 12 1 2
python3 plot_maps_bias_sim_obs_sesa_pr_seas.py --var pr --seas 6 7 8

python3 plot_maps_bias_sim_obs_sesa_t2m_uv10.py --var tas
python3 plot_maps_bias_sim_obs_sesa_t2m_uv10_seas.py --var tas --seas 12 1 2
python3 plot_maps_bias_sim_obs_sesa_t2m_uv10_seas.py --var tas --seas 6 7 8

python3 plot_maps_bias_sim_obs_sesa_t2m_uv10.py --var sfcWind
python3 plot_maps_bias_sim_obs_sesa_t2m_uv10_seas.py --var sfcWind --seas 12 1 2
python3 plot_maps_bias_sim_obs_sesa_t2m_uv10_seas.py --var sfcWind --seas 6 7 8

python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_pr.py --var pr
python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_pr_seas.py --var pr --seas 12 1 2
python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_pr_seas.py --var pr --seas 6 7 8

python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_t2m_uv10.py --var tas
python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_t2m_uv10_seas.py --var tas --seas 12 1 2
python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_t2m_uv10_seas.py --var tas --seas 6 7 8

python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_t2m_uv10.py --var sfcWind
python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_t2m_uv10_seas.py --var sfcWind --seas 12 1 2
python3 plot_heatmap_diurnal_cycle_sim_obs_sesa_t2m_uv10_seas.py --var sfcWind --seas 6 7 8

python3 plot_graph_pdf_stats_sim_obs_sesa_pr_mmd.py --var pr
python3 plot_graph_pdf_stats_sim_obs_sesa_pr_mmh.py --var pr
python3 plot_graph_pdf_stats_sim_obs_sesa_t2m_uv10.py --var tas
python3 plot_graph_pdf_stats_sim_obs_sesa_t2m_uv10.py --var sfcWind

}
