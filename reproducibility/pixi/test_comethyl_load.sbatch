#!/bin/bash
#SBATCH --job-name=test_comethyl
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=test_comethyl_%j.out
#SBATCH --error=test_comethyl_%j.err

set -euo pipefail

source ~/.profile

export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME

unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

mkdir -p "$PIXI_CACHE_DIR"

cd /quobyte/lasallegrp/projects/CHDS/WGBS/2025_intergen_BrCa_comethyl_George/reproducibility/pixi

hostname
date

time timeout 600 pixi run Rscript --vanilla -e '
cat("Before Rsamtools\n")
library(Rsamtools)
cat("After Rsamtools\n")

cat("Before comethyl\n")
library(comethyl)
cat("After comethyl\n")
'

date
