## Genotype comparisons of histopathology scores: Dp16 vs. WT mice

## Objective

This analysis project is part of the [Trisomy 21 Model Atlas](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html).

This workflow compares **histopathology** scores between Dp(16)1Yey (Dp16) and wild-type (WT) mice across tissues at 1 month and 4 months of age (adipose was graded at 4 months only). Tissues were graded by blinded veterinary pathologists using tissue-specific component criteria, and component scores were summed within each tissue to give a composite score. Composite and component scores are compared between genotypes separately for each tissue-timepoint using one-sided, unpaired Wilcoxon rank-sum tests (Dp16 > WT), with rank-based effect sizes (r). Effect sizes are reported both as magnitude (r) and with direction (signed r; positive = higher scores in Dp16).

p values are nominal (not adjusted across tissues or components), and group sizes are small (12 mice per genotype per timepoint), so results are best treated as exploratory.

Please refer to the top-level `README.md` in the `t21-model-atlas-analysis/` repository for a full overview of all analyses and general data access instructions.

------------------------------------------------------------------------

## Repository contents  
```         
histo_mouse_Dp16_tissue-timepoint/ 
  ├── histo_mouse_Dp16_tissue-timepoint.Rproj   # RStudio project file
  ├── histo_mouse_Dp16_tissue-timepoint.R       # Main analysis script 
  ├── helper_functions_histo.R                  # Custom R functions used in analysis 
  ├── data/                                     # Input datasets 
  ├── results/                                  # Test results and summary tables 
  ├── plots/                                    # Generated plots 
  ├── rdata/                                    # Workspace images and RDS objects 
  ├── renv.lock                                 # Reproducible package versions 
  └── README.md                                 # This README file
```

------------------------------------------------------------------------

## System Requirements 

The R packages used in this analysis can be run on any standard computer with enough RAM to support the operations.

This analysis was originally run on a MacBook Pro running macOS 14.4+ and R version 4.4+.

The `renv` package can be used to manage the R environment.

Exact versions of all R packages can be found in the `renv.lock` file.

------------------------------------------------------------------------

## Data Sources 

The following files are provided in the `data/` directory:

* Histopathology component and composite scores for 1-month (`T21_atlas_Dp16_histo_scores_1mo_v1.0.tsv`) and 4-month (`T21_atlas_Dp16_histo_scores_4mo_v1.0.tsv`) mice, and data dictionary (`T21_atlas_Dp16_histo_scores_v1.0_dictionary.tsv`); scores are recorded per animal (`Animal_No`) within each age group
* Sample metadata (`T21_atlas_Dp16_metadata_v1.0.tsv`) and data dictionary (`T21_atlas_Dp16_metadata_v1.0_dictionary.tsv`), used for genotype and sex; animals are matched to scores by age and `Animal_No`, built from `Mouse_ID` and `Cage_ID`

Grading criteria and scales for each tissue are given in Table S14 of the manuscript.

Atlas datasets, including histopathology scores, can also be explored in the [Trisomy 21 Model Atlas collection](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html) on the INCLUDE Experimental Models of Down Syndrome portal.

No pre-processing beyond what is described in the manuscript is required prior to running the script.

------------------------------------------------------------------------

## Outputs

Output file names start with the prefix set in the script (`out_file_prefix`, `histo_mouse_Dp16_tissue-timepoint_`).

* `results/`
  * Composite score Wilcoxon tests with effect sizes (`*_composite_wilcox.txt`)
  * Component score Wilcoxon tests with effect sizes (`*_component_wilcox.txt`)
  * Both tables in one Excel workbook (`*_wilcox_results.xlsx`)
  * R session information (`*_session_info.txt`)
* `plots/`: individual composite and component score plots by genotype, and composite and component effect size heatmaps (magnitude and signed)
* `rdata/`: saved workspace (`.RData`)

------------------------------------------------------------------------

## R Environment Setup and Running Analyses  

1.  Clone the repository.

    ```         
    git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/t21-model-atlas-analysis.git
    ```

2.  Change to the `histo_mouse_Dp16_tissue-timepoint/` directory and open the R project via the `.Rproj` file.

3.  Set up reproducible R environment (requires `renv` package to be installed).

    Option A. Restore the R environment.\
    This will install the exact versions of all R packages but requires matching R version.

    ```         
    install.packages("renv")
    renv::restore()
    ```

    Option B. Initialize the R environment.\
    This will install all R packages but will not ensure identical versions.

    ```         
    install.packages("renv")
    renv::init(bioconductor = TRUE)
    ```

4.  Follow the workflow in `histo_mouse_Dp16_tissue-timepoint.R`.

### Helper functions

The `helper_functions_histo.R` script contains functions used throughout the analysis, including:

-   Custom ggplot theme setup for consistent figure formatting
-   A function for individual score plots by genotype
-   A function to export results tables as Excel workbooks

These functions are customized for this project and require no modification for standard execution of the workflow.
