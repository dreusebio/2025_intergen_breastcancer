#!/bin/bash
#SBATCH --job-name=Figures_Intergen_CHDS
#SBATCH --mail-user=gekuodza@health.ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=23:00:00
#SBATCH --partition=low
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=250G
#SBATCH --account=publicgrp #genome-center-grp, publicgrp ,lasallegrp
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err
#SBATCH --chdir=/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George

set -euo pipefail

###################
# Run Information #
###################
start=$(date +%s)
hostname
date

####################
# Load Environment #
####################
source ~/.profile

# Pixi settings for cluster stability
export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME

mkdir -p "$PIXI_CACHE_DIR"

# Keep R isolated from user libraries/configs
unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Optional threading hygiene
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# # Move into Pixi workspace
# cd /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/reproducibility/pixi

conda activate Comethyl_v1.3.0

# # Sanity check
# Rscript --vanilla -e 'cat("R version:", R.version.string, "\n"); library(comethyl); cat("comethyl OK\n")'

##############
# Run Script #
##############

# ##############
# # Figures Scripts#
# ##############
# # # Figure1a

Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/figures/figure1_panelA_sv_black_race.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --pc_meta /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/protected_model/PC_Metadata.tsv \
  --bicor_stats /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/protected_model/PC_Trait_Correlation_Stats_Bicor.tsv \
  --p_thresh 0.05 \
  --fig_width_mm 168 \
  --fig_height_mm 65 \
  --dpi 600

# # # Figure1b

Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/figures/figure1_panelB_heatmap.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --trait_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/00_figure1_heatmap_primary_traits.txt \
  --fig_width_mm 168 \
  --fig_height_mm 140 \
  --dpi 600

# # # Figure1c
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/figures/figure1_panelC_lollipop.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --annotated_regions /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Annotated_Regions.xlsx \
  --modules yellow,red,saddlebrown,brown \
  --top_n 10 \
  --fig_width_mm 168 \
  --fig_height_mm 90 \
  --dpi 600 \
  --color_overrides "yellow=darkgoldenrod3"

#figure 2-5
#indivdual module-trait plots
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/figures/figure2-5.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/07_module_detection/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules.rds \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
  --set_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/00_figure1_heatmap_primary_traits.txt \
  --module_figure_map /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/figures/module_figure_map.tsv \
  --output_mode individual \
  --p_thresh 0.05 \
  --dpi 600

#multipanel module-trait plots
Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/figures/figure2-5.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/07_module_detection/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules.rds \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026_full_dataset_priority_numeric.xlsx \
  --set_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/00_figure1_heatmap_primary_traits.txt \
  --trait_config /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/figures/trait_config.tsv \
  --module_figure_map /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/figures/module_figure_map.tsv \
  --modules_filter yellow,red,saddlebrown,brown \
  --output_mode multipanel \
  --p_thresh 0.05 \
  --fig_width_mm 168 \
  --fig_height_mm 200 \
  --base_size 8 \
  --dpi 600

end=$(date +%s)
echo "Run completed in $((end - start)) seconds"
date