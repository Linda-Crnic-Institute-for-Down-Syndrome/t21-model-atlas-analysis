## DESeq2 analysis of bulk RNA-seq for individual iPSC-derived cell types: T21 vs. D21

## Objective

This analysis project is part of the [Trisomy 21 Model Atlas](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html).

This workflow uses DESeq2 for differential expression analysis of **bulk RNA-seq** data from human induced pluripotent stem cell (iPSC) lines with trisomy 21 (T21) and disomic (D21) controls. The lines were derived from Human Trisome Project participants. Each cell type (undifferentiated iPSCs or an iPSC-derived cell type) is analyzed separately, comparing T21 with D21. The preferred model adjusts for surrogate variables estimated with SVA. Simple (genotype only) and multivariable (sex, RIN) models are included for comparison.

As written, this script runs the analysis workflow for **undifferentiated iPSCs**. The same workflow was applied to each iPSC-derived cell type in the Atlas. Section 0.3 of the script describes what to change to analyze another cell type.

Please refer to the top-level `README.md` in the `t21-model-atlas-analysis/` repository for a full overview of all analyses and general data access instructions.

------------------------------------------------------------------------

## Repository contents  
```         
bulk_human_iPSC_cell-type/ 
  ├── bulk_human_iPSC_cell-type.Rproj   # RStudio project file
  ├── bulk_human_iPSC_cell-type.R       # Main analysis script 
  ├── helper_functions_DESeq.R          # Custom R functions used in analysis 
  ├── data/                             # Input datasets (counts and RPKMs not included in repository) 
  ├── results/                          # Model outputs and summary tables 
  ├── plots/                            # Generated plots 
  ├── rdata/                            # Workspace images and RDS objects 
  ├── renv.lock                         # Reproducible package versions 
  └── README.md                         # This README file
```

------------------------------------------------------------------------

## System Requirements 

The R packages used in this analysis can be run on any standard computer with enough RAM to support the operations.

This analysis was originally run on a MacBook Pro running macOS 14.4+ and R version 4.4+.

The `renv` package can be used to manage the R environment.

Exact versions of all R packages can be found in the `renv.lock` file.

Parallel processing uses `BiocParallel::MulticoreParam()`, which is not available on Windows. Windows users can replace it with `SnowParam()` in section 0.1 of the script.

------------------------------------------------------------------------

## Data Sources 

Bulk RNA-seq data for each cell type can be obtained from the Gene Expression Omnibus (GEO):

* Gene-level counts and RPKMs: GEO SuperSeries [GSE347972](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347972), with one SubSeries per cell type (undifferentiated iPSCs: [GSE347268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347268); see `GEO_SubSeries` in the metadata for other cell types)

The following files are provided in the `data/` directory:

* Sample metadata for all iPSC datasets (`T21_atlas_iPSC_metadata_v1.2.tsv`) and data dictionary (`T21_atlas_iPSC_metadata_v1.2_dictionary.tsv`)
* Gene annotation derived from the [GENCODE human release 33](https://www.gencodegenes.org/human/release_33.html) basic annotation GTF (`gene_annotation_Gencode.v33.basic.txt.gz`)
* [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) human Hallmark gene sets v7.4 (`h.all.v7.4.symbols.gmt`)

Download the required files and place them in the `data/` directory before running the analysis. Section 1.1 of the script includes optional code to download the counts and RPKMs directly from GEO using the `GEOquery` package.

Atlas datasets, including complete differential expression results for every cell type, can also be explored in the [Trisomy 21 Model Atlas collection](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html) on the INCLUDE Experimental Models of Down Syndrome portal.

No pre-processing beyond what is described in the manuscript is required before running the script.

------------------------------------------------------------------------

## Outputs

Output file names start with the cell-type prefix set in the script (`out_file_prefix`, e.g. `bulk_human_iPSC_undifferentiated_`).

* `results/`
  * DESeq2 results tables for each model (`*_results_<model>_T21_vs_D21_0.5cpm.txt` and `.xlsx`); the SVA model (`results_sva`) is the preferred model
  * GSEA Hallmark results (`*_GSEA_Hallmarks.xlsx`)
  * R session information (`*_session_info.txt`)
* `plots/`: QC, PCA, model comparison, MA, volcano, Manhattan, individual gene, heatmap, and GSEA plots
* `rdata/`: saved workspace (`.RData`)

------------------------------------------------------------------------

## R Environment Setup and Running Analyses  

1.  Clone the repository.

    ```         
    git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/t21-model-atlas-analysis.git
    ```

2.  Change to the `bulk_human_iPSC_cell-type/` directory and open the R project via the `.Rproj` file.

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

4.  Follow the workflow in `bulk_human_iPSC_cell-type.R`.

### Helper functions

The `helper_functions_DESeq.R` script is shared with the Dp16 bulk RNA-seq project and contains functions used throughout the analysis, including:

-   Custom ggplot theme setup for consistent figure formatting
-   Functions to generate and export DESeq2 results tables with log2 fold-change shrinkage
-   Functions for QC and visualization: PCA, sample clustering, dispersion, MA, and volcano plots
-   Functions to run and plot GSEA (combined positive and negative enrichment)

These functions are customized for this project and require no modification for standard execution of the workflow.

Note: triplicated genes are identified as all genes on human chromosome 21 (set in section 0.2 of the script, and in the helper `volcano_plot_chr21()`).
