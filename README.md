# intergen_BrCa_wgbs_comethyl

Bioinformatics analysis of blood whole-genome bisulfite sequencing data from the CHDS Disparities Study to identify co-methylation signatures associated with in utero chemical exposures, neighborhood stressors, cumulative lifetime stress, and breast cancer-related risk factors.

---

## Project summary

**Data type:** Whole Genome Bisulfite Sequencing (WGBS)

**Study / cohort:** Child Health and Development Studies (CHDS) — Disparities Study

**Genome build:** hg38

**Biological sample:** Blood

**Generation profiled:** F1 daughters at midlife

**Sample size:** 173 F1 participants

**Starting input for this repository:** Bismark cytosine reports generated upstream

**Downstream analysis framework:** Comethyl / WGCNA-based co-methylation analysis

**Primary outputs:** Co-methylation modules, module eigengene–trait associations, annotated regions, figures, and enrichment results

---

## Project overview

This project investigates whether in utero exposure to environmental chemicals, neighborhood disadvantage, and cumulative lifetime stress are associated with midlife blood DNA methylation signatures related to breast cancer risk.

The study uses existing blood specimens from F1 daughters in the Child Health and Development Studies (CHDS) Disparities Study. In utero chemical exposures were measured in mothers' pregnancy blood using high-resolution metabolomics and targeted biomonitoring measures. Neighborhood context was captured using Index of Concentration at the Extremes (ICE) measures across education, race, income, and race-by-income dimensions. Cumulative stress burden was assessed using allostatic load and related cardiometabolic, inflammatory, and neurodegenerative biomarkers.

This repository focuses on the downstream WGBS cytosine report to co-methylation analysis workflow. Raw WGBS preprocessing was completed upstream, and the final Bismark cytosine reports were used as input for this Comethyl analysis.

The central hypothesis is that in utero exposures to DDT/DDE, PCBs, PFAS, and other organochlorine compounds, together with neighborhood disadvantage and cumulative lifetime stress, are associated with midlife DNA methylation signatures in pathways relevant to breast cancer risk. These effects may be more pronounced or biologically distinct in Black women.

---

## Study design

| Feature               | Detail                                                                          |
| --------------------- | ------------------------------------------------------------------------------- |
| Design                | Cross-sectional methylome analysis using midlife F1 blood samples               |
| Parent cohort         | Child Health and Development Studies                                            |
| Substudy              | CHDS Disparities Study                                                          |
| Generation profiled   | F1 daughters                                                                    |
| Sample size           | 173 participants                                                                |
| Biological sample     | Blood                                                                           |
| Exposure window       | In utero / prenatal period                                                      |
| Main exposure domains | DDT/DDE, PCBs, PFAS, organochlorines, ICE neighborhood indices, allostatic load |
| Main outcome          | Genome-wide DNA methylation                                                     |
| Statistical framework | Comethyl / WGCNA-based co-methylation analysis                                  |

---

## Trait variable domains

| Domain                                   | Description                                                                        |
| ---------------------------------------- | ---------------------------------------------------------------------------------- |
| Demographics                             | Race/ethnicity, age, batch, technical variables                                    |
| Breast cancer risk factors               | Obesity, overweight status, age at menarche                                        |
| Breast cancer outcomes                   | F1 breast cancer status and family breast cancer history where available           |
| DDT/DDE exposures                        | Continuous and categorical measures of DDT, DDE, and o,p'-DDT                      |
| PCB exposures                            | Individual PCB congeners                                                           |
| PFAS exposures                           | Individual PFAS compounds                                                          |
| Other organochlorines                    | HCB, oxychlordane, β-BHC, and related compounds                                    |
| Neighborhood context                     | ICE measures for education, race, income, and race-by-income dimensions            |
| Cumulative stress                        | Allostatic load summary score and component biomarkers                             |
| Lipid markers                            | Total cholesterol, triglycerides, total lipids                                     |
| Cardiometabolic and inflammatory markers | hsCRP, DHEAS, HDL, IL-6, HbA1c, blood pressure, BMI, body fat, waist circumference |
| Neurodegeneration / AD biomarkers        | GFAP, NfL, Aβ40/42, Tau, pTau-217                                                  |
| Kidney function                          | Creatinine, cystatin C                                                             |
| Genetic risk markers                     | APOE genotype and ε4 carrier status                                                |

---

## Analysis workflow

This repository implements the downstream WGBS cytosine report to co-methylation analysis workflow.

### Upstream preprocessing

Raw WGBS preprocessing was completed outside this repository.

The upstream preprocessing included:

* FASTQ processing
* Adapter trimming
* Bisulfite read alignment
* Methylation extraction
* Generation of Bismark cytosine reports
* Quality-control summaries

The final Bismark cytosine report files were used as the starting input for this repository.

### 1. Import cytosine reports

* Import Bismark cytosine reports
* Match cytosine reports to sample metadata
* Construct methylation objects for downstream Comethyl analysis

### 2. CpG and region filtering

* CpG-level coverage filtering
* CpG clustering into genomic regions
* Filtering regions by coverage and variability thresholds

### 3. Region-level methylation

* Aggregation of CpG methylation values into region-level methylation
* Construction of region-by-sample methylation matrices
* Preparation of methylation matrices for downstream network analysis

### 4. Adjustment

* Principal component analysis
* Identification of technical or unwanted variation
* PC-based and/or surrogate-variable-informed correction
* Sensitivity analyses excluding biologically relevant exposure or outcome variables from adjustment where appropriate

### 5. Network analysis

* Soft-thresholding power selection
* Co-methylation network construction
* Module detection using a WGCNA-like framework
* Calculation of module eigengenes

### 6. Association testing

* Module eigengene associations with in utero chemical exposures
* Module eigengene associations with ICE neighborhood measures
* Module eigengene associations with allostatic load and biomarker traits
* Module eigengene associations with obesity, age at menarche, and breast cancer-related traits
* Race-stratified or interaction analyses where applicable
* Correlation analyses using bicor and/or Pearson correlation
* Regression models for selected traits
* Multiple-testing correction

### 7. Downstream interpretation

* Annotation of genomic regions
* Gene mapping
* Functional enrichment analysis
* Pathway interpretation for breast cancer risk, endocrine disruption, inflammation, stress biology, metabolism, and environmental exposure pathways
* Prioritization of modules linked to exposures and breast cancer-related risk factors

---

## Repository structure

```text
intergen_BrCa_wgbs_comethyl/
├── data/            # raw + processed data; not committed
├── scripts/         # analysis scripts + SLURM scripts
├── analysis/        # downstream analyses + configuration files
├── docs/            # methods, notes, and documentation
├── results/         # analysis outputs and figures
├── tests/           # test data and example workflows
├── reproducibility/ # environment setup files
└── README.md
```

---

## Quick start: reproducible environment

This project uses Pixi for reproducible software environments.

### 1. Clone repository

```bash
git clone https://github.com/dreusebio/intergen_BrCa_wgbs_comethyl.git
cd intergen_BrCa_wgbs_comethyl
```

### 2. Configure Pixi environment

For general use:

```bash
export PIXI_HOME=/path/to/your/pixi_home
export PIXI_CACHE_DIR=/scratch/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME
```

For UC Davis Hive users:

```bash
export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME
```

### 3. Install environment

```bash
cd reproducibility/pixi
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

### 4. Enable Bioconductor post-link scripts if needed

Some Bioconductor annotation database packages require post-link scripts. If packages such as `GenomeInfoDbData`, `GO.db`, `org.Hs.eg.db`, or TxDb packages fail to load, run:

```bash
pixi config set --local run-post-link-scripts insecure
pixi install
```

### 5. Install Comethyl

```bash
pixi run install-comethyl
```

### 6. Test installation

```bash
pixi run Rscript -e 'library(comethyl); packageVersion("comethyl")'
```

Expected:

```r
[1] "1.3.0"
```

---

## Running analysis

All scripts should be executed through Pixi.

Example:

```bash
pixi run Rscript scripts/<your_script>.R
```

For cluster jobs, use the SLURM scripts provided in:

```text
scripts/slurm/
```

---

## Run test data

You can validate that the pipeline is working by running the included test dataset.

The test data and example workflows are located in:

```text
tests/
```

Example:

```bash
pixi run Rscript tests/scripts/comethyl_scripts/00_import_cpg_reports.R
```

If no error occurs, the expected test results will be written to:

```text
tests/comethyl_output/
```

---

## Reproducibility

This project is designed to be reproducible using the files in:

```text
reproducibility/
├── env/
│   ├── environment.portable.yml
│   └── install_comethyl.R
└── pixi/
    ├── pixi.toml
    └── pixi.lock
```

### Key reproducibility features

* `pixi.lock` guarantees identical software environments
* `install_comethyl.R` installs Bioconductor and GitHub dependencies
* Comethyl is pinned to a specific version or commit
* Analysis scripts are intended to be run through the locked Pixi environment
* Large raw and intermediate data files are excluded from GitHub

---

## Platform support

| Platform            | Support                     |
| ------------------- | --------------------------- |
| Linux / cluster     | Full                        |
| macOS Intel         | Full                        |
| macOS Apple Silicon | Partial; Docker recommended |
| Windows             | Use WSL2 or Docker          |

---

## Troubleshooting

### Pixi install fails because of network or cache issues

Set the Pixi cache to a local temporary directory:

```bash
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
```

Then rerun:

```bash
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

### Reinstall environment

```bash
rm -rf /tmp/$USER/pixi-cache
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

### Missing R or Bioconductor packages

Run:

```bash
pixi run install-comethyl
```

### Bioconductor database packages fail to load

Some Bioconductor annotation database packages require Pixi post-link scripts. Run:

```bash
pixi config set --local run-post-link-scripts insecure
pixi install
```

### Comethyl does not load

Test directly:

```bash
pixi run Rscript -e 'library(comethyl); sessionInfo()'
```

**Note:** Bioconductor package installation steps that require internet access should be run from a login node, not a compute node.

---

## Data availability

Raw WGBS data, cytosine reports, exposure files, phenotype files, biomarker files, and participant-level metadata are not committed to this repository because they may contain sensitive or restricted research information.

This repository contains scripts, documentation, configuration files, and reproducibility resources needed to rerun the analysis when appropriate data access permissions are available.

---

## License

MIT License. See `LICENSE`.

---

## Environment details

**R version:** 4.3.3
**Comethyl version:** 1.3.0
**Installed from:** GitHub
**Environment manager:** Pixi
**Locked environment:** `pixi.lock`

---

## Citation

If using this workflow, please cite or acknowledge:

* Comethyl for co-methylation network analysis
* Bismark for bisulfite read alignment and methylation extraction
* WGCNA
* Bioconductor
* Relevant WGBS and DNA methylation analysis tools
* Child Health and Development Studies resources where appropriate

---

## Contact

**George Eusebio Kuodza**
Postdoctoral Scholar
University of California, Davis
[gekuodza@health.ucdavis.edu](mailto:gekuodza@health.ucdavis.edu)

---

## Future directions

* Add Docker container for full portability
* Add Snakemake or Nextflow workflow integration
* Add automated test workflows
