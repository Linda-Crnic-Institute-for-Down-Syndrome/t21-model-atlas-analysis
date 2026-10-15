## DESeq2 analysis of bulk RNA-seq against histopathology score: Dp16 mice

## Objective

This analysis project is part of the [Trisomy 21 Model Atlas](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html).

This workflow models **bulk RNA-seq** gene expression in Dp(16)1Yey (Dp16) mice as a function of the composite histopathology score measured in the same animals, using DESeq2. Composite score is the sum of tissue-specific component scores from blinded histopathological grading, and is treated as a continuous predictor within the Dp16 genotype. Log2 fold-changes are therefore reported per one-unit increase in composite score. The preferred model adjusts for surrogate variables estimated with SVA; simple (score only) and multivariable (sex, RIN, age in days) models are included for comparison.

This script as written will run the analysis workflow for **4-month lung**. The same workflow was applied to 4-month kidney; section 0.3 of the script describes what to change to analyze another tissue-timepoint.

Because the analysis is restricted to Dp16 animals, it has limited statistical power and is best suited to identifying candidate genes rather than precise effect-size estimation.

Please refer to the top-level `README.md` in the `t21-model-atlas-analysis/` repository for a full overview of all analyses and general data access instructions.

------------------------------------------------------------------------

## Repository contents  
```         
integrated_mouse_Dp16_bulk-histo/ 
  ├── integrated_mouse_Dp16_bulk-histo.Rproj   # RStudio project file
  ├── integrated_mouse_Dp16_bulk-histo.R       # Main analysis script 
  ├── helper_functions_DESeq.R                 # Custom R functions used in analysis 
  ├── data/                                    # Input datasets (counts and RPKMs not included in repository) 
  ├── results/                                 # Model outputs and summary tables 
  ├── plots/                                   # Generated plots 
  ├── rdata/                                   # Workspace images and RDS objects 
  ├── renv.lock                                # Reproducible package versions 
  └── README.md                                # This README file
```

------------------------------------------------------------------------

## System Requirements 

The R packages used in this analysis can be run on any standard computer with enough RAM to support the operations.

This analysis was originally run on a MacBook Pro running macOS 14.4+ and R version 4.4+.

The `renv` package can be used to manage the R environment.

Exact versions of all R packages can be found in the `renv.lock` file.

Parallel processing uses `BiocParallel::MulticoreParam()`, which is not available on Windows; Windows users can replace this with `SnowParam()` in section 0.1 of the script.

------------------------------------------------------------------------

## Data Sources 

Bulk RNA-seq data for each tissue-timepoint can be obtained from the Gene Expression Omnibus (GEO):

* Gene-level counts and RPKMs: GEO SuperSeries [GSE347972](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347972), with one SubSeries per tissue-timepoint (4-month lung: [GSE347804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347804); 4-month kidney: [GSE347802](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347802))

The following files are provided in the `data/` directory:

* Sample metadata (`T21_atlas_Dp16_metadata_v1.0.tsv`) and data dictionary (`T21_atlas_Dp16_metadata_v1.0_dictionary.tsv`); the `GEO_SubSeries` column gives the SubSeries for each sample
* Composite histopathology scores for 4-month lung and kidney (`bulk_mouse_Dp16_histo_4mo_lung_kidney.txt`), matched to RNA-seq samples by `Sampleid`
* Gene annotation derived from [GENCODE mouse release M24](https://www.gencodegenes.org/mouse/release_M24.html) basic annotation GTF (`gene_annotation_Gencode.vM24.basic.txt.gz`)
* [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) mouse Hallmark gene sets v2022.1.Mm (`mh.all.v2022.1.Mm.symbols.gmt`)

Section 1.1 of the script downloads the counts and RPKMs directly from GEO using the `GEOquery` package.

Atlas datasets, including complete differential expression results and histopathology scores, can also be explored in the [Trisomy 21 Model Atlas collection](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html) on the INCLUDE Experimental Models of Down Syndrome portal.

No pre-processing beyond what is described in the manuscript is required prior to running the script.

------------------------------------------------------------------------

## Outputs

Output file names start with the prefix set in the script (`out_file_prefix`, e.g. `integrated_mouse_Dp16_bulk-histo_4mo-lung_`).

* `results/`
  * DESeq2 results table for the preferred SVA model (`*_results_sva_Score_0.5cpm.txt` and `.xlsx`); log2 fold-changes are per one-unit increase in composite score
  * GSEA Hallmark results (`*_GSEA_Hallmarks.xlsx`)
  * R session information (`*_session_info.txt`)
* `plots/`: QC, PCA, model comparison, volcano, individual gene scatter (expression vs. composite score), heatmap (samples ordered by composite score), and GSEA plots
* `rdata/`: saved workspace (`.RData`)

------------------------------------------------------------------------

## R Environment Setup and Running Analyses  

1.  Clone the repository.

    ```         
    git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/t21-model-atlas-analysis.git
    ```

2.  Change to the `integrated_mouse_Dp16_bulk-histo/` directory and open the R project via the `.Rproj` file.

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

4.  Follow the workflow in `integrated_mouse_Dp16_bulk-histo.R`.

### Helper functions

The `helper_functions_DESeq.R` script is shared with the other bulk RNA-seq projects and contains functions used throughout the analysis, including:

-   Custom ggplot theme setup for consistent figure formatting
-   Functions to generate and export DESeq2 results tables with log2 fold-change shrinkage, for both group comparisons and continuous predictors (`get_results_tbl_names()`)
-   Functions for QC and visualization: PCA, dispersion, and volcano plots
-   Functions to run and plot GSEA (combined positive and negative enrichment)

These functions are customized for this project and require no modification for standard execution of the workflow.

Note: identification of triplicated genes uses the coordinates of the Dp(16)1Yey duplicated region on mouse chromosome 16 (set in section 0.2 of the script). These are specific to the Dp16 model and must be changed for other mouse models.
