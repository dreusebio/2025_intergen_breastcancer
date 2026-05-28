#!/bin/bash
#SBATCH --job-name=Cometh_Intergen_CHDS
#SBATCH --mail-user=gekuodza@ucdavis.edu,smholesinger@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=23:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH --account=lasallegrp #genome-center-grp, publicgrp
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

# Move into Pixi workspace
cd /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/reproducibility/pixi

# Sanity check
pixi run Rscript --vanilla -e 'cat("R version:", R.version.string, "\n"); library(comethyl); cat("comethyl OK\n")'

##############
# Run Script #
##############

##############
# Comethyl Scripts#
##############
# Script 00
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/00_import_cpg_reports.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --meta_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
  --cyto_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/processed/cytosine_reports \
  --workers 4

# Script 01
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/01_filter_cpgs_call_regions.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/00_cpg_extraction/Unfiltered_BSseq.rds \
  --cov 3 \
  --per_sample 0.75

#Script 02
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/02_filter_regions.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_regions /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/01_filter_regions/cov3_75pct/Regions.txt \
  --cov_min 4 \
  --meth_sd 0.08


# Script 03
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/03_region_methylation.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_bs /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/01_filter_regions/cov3_75pct/Filtered_BSseq.rds \
  --input_regions /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/02_filter_regions/cov3_75pct/covMin10_methSD0p07/Filtered_Regions.txt


# Script 04
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/04_pcs_and_trait_diagnostics.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
  --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
  --threads 4



# Script 05
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/05_adjust_region_methylation.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/03_region_methylation/cov3_75pct/covMin10_methSD0p07/Region_Methylation.rds \
  --input_pcs /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/PCs.rds \
  --pc_trait_stats /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin10_methSD0p07/PC_Trait_Correlation_Stats_Bicor.tsv \
  --protected_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/technical_traits.txt \
  --run_v1 TRUE \
  --run_v2 TRUE \
  --run_v3 TRUE \
  --v2_cor_method bicor \
  --v2_p_thresh 0.05 \
  --v3_p_thresh 0.05


# Script 06
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/06_soft_power.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --input_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
  --input_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --input_v3 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v3_technical_pcs_only/Adjusted_Region_Methylation_technicalPCs_only.rds \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 20 \
  --threads 4

# Script 07
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/07_module_detection.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --regions_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/02_filter_regions/cov3_75pct/covMin10_methSD0p07/Filtered_Regions.txt \
  --input_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
  --softpower_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/06_soft_power/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/SoftPower_pearson.rds \
  --input_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --softpower_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/06_soft_power/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/SoftPower_pearson.rds \
  --module_cor pearson \
  --deep_split 4 \
  --min_module_size 10 \
  --merge_cut_height 0.1 \
  --threads 4

# Script 08
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/08_module__sample_diagnostics.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/07_module_detection/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules.rds \
  --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/07_module_detection/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules.rds \
  --module_cor pearson \
  --sample_dendro_distance euclidean \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx

 
#script 09a
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/09a_me_trait_analysis.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173_05222026.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --run_bicor TRUE \
  --p_thresh 0.05 \
  --top_n 25 \
  --module_dendro_distance bicor \

#script 09b
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/09b_me_trait_presentation_heatmaps.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --modules_rds_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --set_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --p_thresh 0.05 \
  --top_n 25 \
  --min_sig_traits 1 \
  --focused_width 14 \
  --focused_height 9 \
  --module_order_mode exposure_outcome \
  --module_order_trait_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/module_ordering/exposure_outcome_traits.tsv \
  --module_order_p_thresh 0.05 \
  --module_order_score p \
  --trait_order_mode exposure_outcome \
  --drop_grey_modules TRUE

# #Script 10
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/10_circos_and_module_summaries.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --genome hg38 \
  --single_module_color red \
  --all_modules_width 4.0 \
  --all_modules_height 4.0 \
  --single_module_width 5.5 \
  --single_module_height 5.5 \
  --counts_width 7 \
  --counts_height 6


# #Script 11
#v1
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/11_me_trait_pair_plots.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
  --set_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --p_thresh 0.05 \
  --only_significant TRUE \
  --cat_plot_style violin \
  --top_n_pairs 50

#v2
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/11_me_trait_pair_plots.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --modules_rds /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --stats_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
  --set_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --p_thresh 0.05 \
  --only_significant TRUE \
  --cat_plot_style violin \
  --top_n_pairs 50

#script 12a
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12a_annotate_modules.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --sample_info /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/Modules_for_downstream.rds \
  --modules_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/08_module_and_sample_diagnostics/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Modules_for_downstream.rds \
  --genome hg38 \
  --annotation_mode auto

#script 12b
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12b_enrichment.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v1_all_pcs \
  --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs \
  --trait_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets \
  --stats_file_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --p_thresh 0.05 \
  --do_kegg false \
  --do_go false \
  --do_reactome false \
  --do_enrichr true \
  --enrichr_dbs GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,MSigDB_Hallmark_2020,MSigDB_Oncogenic_Signatures,ClinVar_2025,WikiPathways_2024_Human \
  --use_enrichr_background false \
  --enrichr_retries 5 \
  --enrichr_sleep 2 \
  --plot_p_col pvalue \
  --plot_top_n 10 \
  --plot_sig_cutoff 0.05

 #full list
  #  --enrichr_dbs ClinVar_2025,WikiPathways_2024_Human,DisGeNET,GTEx_Tissue_Expression_Up,Human_Phenotype_Ontology \
#  --trait_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/presentation_sets/12_heatmap_primary_traits.txt \
#script 12b just one text file only
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/12b_enrichment.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v1_all_pcs \
  --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs \
  --stats_file_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --stats_file_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --module_list_file /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/config/enrichrment_modules/modules.txt \
  --p_thresh 0.05 \
  --do_kegg false \
  --do_go false \
  --do_reactome false \
  --do_enrichr true \
  --enrichr_dbs GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,MSigDB_Hallmark_2020,MSigDB_Oncogenic_Signatures,ClinVar_2025,WikiPathways_2024_Human,GTEx_Tissue_Expression_Up,Reactome_Pathways_2024,Human_Phenotype_Ontology,Epigenomics_Roadmap_HM_ChIP-seq,Elsevier_Pathway_Collection,Panther_2016,UK_Biobank_GWAS_v1,GWAS_Catalog_2025 \
  --use_enrichr_background false \
  --enrichr_retries 5 \
  --enrichr_sleep 2 \
  --plot_p_col pvalue \
  --plot_top_n 10 \
  --plot_sig_cutoff 0.05 \
  --plot_width 10 \
  --plot_height 7

#p.adjust", "pvalue
# Script 13a
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/comethyl_scripts/13A_core_regulatory_state_overlap.R \
  --project_root /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George \
  --region_version covMin10_methSD0p07 \
  --annotation_dir_v2 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/covMin10_methSD0p07/v2_exclude_protected_pcs \
  --modules turquoise \
  --source roadmap 

  --annotation_dir_v1 /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/comethyl_output/12_annotation/covMin10_methSD0p07/v1_all_pcs \

# Script 14
pixi run Rscript /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/scripts/tables/00_table1.R \
  --input /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/metadata/Intergen173.xlsx \
  --output /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/00_table1/table1_intergen.xlsx \
  --output_pdf /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/results/00_table1/table1_intergen.pdf \
  --verify

end=$(date +%s)
echo "Run completed in $((end - start)) seconds"
date