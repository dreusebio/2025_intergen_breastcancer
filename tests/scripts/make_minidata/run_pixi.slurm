#!/bin/bash
#SBATCH --job-name=Make_mini
#SBATCH --mail-user=gekuodza@ucdavis.edu,smholesinger@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=23:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=700G
#SBATCH --account=lasallegrp
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err
#SBATCH --chdir=/quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/tests/scripts/make_minidata

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

##############
# Run Script #
##############

python /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/tests/scripts/make_minidata/00_make_comethyl_minidata.py \
  --input_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/processed/cytosine_reports \
  --output_dir /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/tests/data/processed/job \
  --sample_list /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/data/processed/sample_ids.txt \
  --keep_chr chr22 
