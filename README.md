# Trisomy 21 Model Atlas: Analysis Code

Analysis code for the **Trisomy 21 Model Atlas**, an openly accessible resource of transcriptomic and histopathological data from mouse and human iPSC-derived models of Down syndrome.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23198762.svg)](https://doi.org/10.5281/zenodo.23198762)

> **Reproducing published results:** the `main` branch reflects the current state of the Atlas and may include updates made after a given publication. To reproduce the analyses in a specific preprint or paper, use the corresponding tagged release listed under [Releases](#releases). For example, code as used in [Rossmassler, Timkovich, Niemeyer et al. (bioRxiv 2026)](https://doi.org/10.64898/2026.09.20.751055) is release **v1.0.0** ([VERSION_DOI]).

------------------------------------------------------------------------

## Overview
* [About the Atlas](#about-the-atlas)
* [Releases](#releases)
* [Repository Structure](#repository-structure-and-table-of-contents)
* [Data Sources](#data-sources)
* [Upstream Processing](#upstream-processing)
* [Software & Dependencies](#software--dependencies)
* [R Environment Setup and Running Analyses](#r-environment-setup-and-running-analyses)
* [Citation & License](#citation--license)

------------------------------------------------------------------------

## About the Atlas

The Trisomy 21 Model Atlas maps the transcriptional and pathological consequences of trisomy 21 across tissues, developmental stages, sexes, and species. The initial release comprises:

* **Mouse:** bulk RNA-seq and blinded histopathology from Dp(16)1Yey (Dp16) and wild-type mice across seven tissues (blood, brain, heart, kidney, intestine, liver, lung) at three timepoints (E18.5, 1 month, 4 months).
* **Human:** bulk RNA-seq from a panel of 12 [Human Trisome Project](https://www.trisome.org/) iPSC lines (6 trisomic, 6 disomic), undifferentiated and differentiated into five cell types (iHSPCs, iMonocytes, iHepatocytes, iNeurons, iAstrocytes).

Future releases will add single-cell transcriptomics, multiplexed imaging, and additional mouse models and iPSC-derived systems. Atlas data can be explored interactively in the [Trisomy 21 Model Atlas collection](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html) on the INCLUDE Experimental Models of Down Syndrome portal.

This repository contains the R-based secondary analyses (statistics, integration, and visualization). Primary processing of sequencing data is performed separately; see [Upstream Processing](#upstream-processing).

------------------------------------------------------------------------

## Releases

| Version | Date | Atlas data release | Associated publication | Code DOI |
|---|---|---|---|---|
| v1.0.0 | 2026-10-09 | Release 1 | Rossmassler, Timkovich, Niemeyer et al. *bioRxiv* 2026. [doi:10.64898/2026.09.20.751055](https://doi.org/10.64898/2026.09.20.751055) | ()[VERSION_DOI] |

Versioning: major versions correspond to Atlas data releases; minor versions to analysis updates within a release (e.g., manuscript revisions).

------------------------------------------------------------------------

## Repository Structure and Table of Contents

Each analysis workflow is a self-contained R Project. Project folders are named `<datatype>_<species>_<model>_<analysis>`, so projects sort by data type; projects combining data types use the `integrated_` prefix.
<!--
Project naming extends to additional data modalities and models:
- sc_mouse_Dp16_tissue-timepoint, spatial_mouse_Dp16_lung
- bulk_mouse_Ts65Dn_tissue-timepoint, bulk_human_CO_timepoint
-->

```
t21-model-atlas-analysis/
│
├── bulk_mouse_Dp16_tissue-timepoint/             # Self-contained R Project directory
│    ├── bulk_mouse_Dp16_tissue-timepoint.Rproj     # RStudio project file; double click to open in RStudio
│    ├── bulk_mouse_Dp16_tissue-timepoint.R         # Main analysis script
│    ├── helper_functions_DESeq.R                   # Associated R functions
│    ├── data/                                      # Input data (meta data and annotation provided; counts/RPKMs downloaded from GEO)
│    ├── results/                                   # Results tables, processed data, model outputs
│    ├── plots/                                     # Visualizations and plots
│    ├── rdata/                                     # Workspace images and RDS objects
│    ├── renv.lock                                  # R package versions for reproducibility
│    └── README.md                                  # Analysis-specific README
├── bulk_mouse_Dp16_timepoint-comparison/         # Self-contained R Project directory
├── bulk_mouse_Dp16_sex-stratified/               # Self-contained R Project directory
├── bulk_human_iPSC_cell-type/                    # Self-contained R Project directory
├── histo_mouse_Dp16_tissue-timepoint/            # Self-contained R Project directory
├── integrated_mouse_Dp16_bulk-histo/             # Self-contained R Project directory
├── .zenodo.json                                  # Metadata for Zenodo DOI registration
├── LICENSE.md                                    # Software license
└── README.md                                     # This README file
```

### Analysis R Projects

| Project | Description | Inputs |
|---|---|---|
| `bulk_mouse_Dp16_tissue-timepoint` | Dp16 vs. WT differential expression within each tissue–timepoint; cross-tissue comparisons; gene set enrichment | gene-level count and RPKM data |
| `bulk_mouse_Dp16_timepoint-comparison` | Comparisons between timepoints (E18.5, 1 month, 4 months) within each tissue and genotype | gene-level count and RPKM data |
| `bulk_mouse_Dp16_sex-stratified` | Dp16 vs. WT differential expression within each sex at 4 months | gene-level count and RPKM data |
| `bulk_human_iPSC_cell-type` | T21 vs. D21 differential expression in undifferentiated and iPSC-derived cell types, each analyzed separately | gene-level count and RPKM data |
| `histo_mouse_Dp16_tissue-timepoint` | Composite and component histopathology scores; genotype comparisons | Histopathology scores |
| `integrated_mouse_Dp16_bulk-histo` | Gene expression modeled against composite pathology score in 4-month Dp16 kidney and lung | gene-level count and RPKM data; Histopathology scores |

Where a project uses outputs of another, copy the listed files from the upstream project's `results/` into the downstream project's `data/` directory, as described in each project's README.

------------------------------------------------------------------------

## Data Sources

Download each dataset to the appropriate `data/` directory within each R project.

### RNA sequencing
* **Mouse tissues, undifferentiated iPSCs, iHSPCs, iMonocytes:** GEO SuperSeries [GSE347972](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347972), with SubSeries for each mouse tissue–timepoint and iPSC cell type.
* **iHepatocytes:** GEO [GSE296449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296449) (Dunn et al. 2026).
* **iNeurons and iAstrocytes:** GEO [GSE344531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344531) and [GSE344533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344533) (Dooling et al. 2026).

Per-series accessions and sequencing metrics are provided in Table S1 of the preprint.

### Histopathology
* Composite and component histopathology scores: Table S2 and [Zenodo LOCATION / DOI]

### Interactive access
Atlas datasets, including complete differential expression results, can be queried and visualized in the [Trisomy 21 Model Atlas collection](https://experimentalmodels.includedcc.org/trisomy-21-model-atlas.html) on the [INCLUDE Experimental Models of Down Syndrome portal](https://experimentalmodels.includedcc.org), part of the [INCLUDE Data Hub](https://portal.includedcc.org).

------------------------------------------------------------------------

## Upstream Processing

The R projects in this repository start from gene-level count and RPKM matrices. Raw sequencing reads were processed with [BaSH_seq](https://github.com/mattgalbraith/BaSH_seq):

* Quality filtering and trimming: bbduk (BBTools v37.99) and fastq-mcf (ea-utils v1.05)
* Alignment: HISAT2 v2.1.0 to GRCm38 with Gencode M24 annotation (mouse) or GRCh38 with Gencode v33 annotation (human)
* Filtering: Samtools v1.5 (MAPQ > 10)
* Quantification: HTSeq-count v0.6.1 (`--stranded=reverse --minaqual=10 --type=exon --mode=intersection-nonempty`)

------------------------------------------------------------------------

## Software & Dependencies

* [R](https://cran.r-project.org/) and [Bioconductor](https://bioconductor.org/)
* [RStudio](https://posit.co/download/rstudio-desktop)

Key packages include:
* renv
* DESeq2
* sva
* fgsea
* tidyverse
* ggplot2

The `renv.lock` file within each analysis project directory contains a full list of packages and versions.

------------------------------------------------------------------------

## R Environment Setup and Running Analyses

1. Clone the repository.
    ```
    git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/t21-model-atlas-analysis.git
    ```
    To reproduce a specific publication, check out its release (see [Releases](#releases)):
    ```
    git checkout v1.0.0
    ```

2. Change to the desired R Project directory and open the R project via its `.Rproj` file.

3. Set up the reproducible R environment (requires the `renv` package).

    Option A. Restore the R environment.
    This will install the exact versions of all R packages but requires a matching R version.
    ```
    install.packages("renv")
    renv::restore()
    ```

    Option B. Initialize the R environment.
    This will install all R packages but will not ensure identical versions.
    ```
    install.packages("renv")
    renv::init(bioconductor = TRUE)
    ```

4. Follow the workflow in the project's analysis script.

------------------------------------------------------------------------

## Citation & License

If you use this code or the Atlas, please cite:

**Preprint**
Rossmassler K\*, Timkovich A\*, Niemeyer B\*, et al. A Trisomy 21 Model Atlas reveals conserved dosage effects and context-specific transcriptional responses across mouse and human models of Down syndrome. *bioRxiv* 2026. [doi:10.64898/2026.09.20.751055](https://doi.org/10.64898/2026.09.20.751055)

**Code**
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23198762.svg)](https://doi.org/10.5281/zenodo.23198762)  
The DOI above always resolves to the latest release. To cite the exact code used in a publication, use the version DOI listed under [Releases](#releases).

This project is licensed under the MIT License – see the LICENSE file for details.

## Funding

This work was supported primarily by NIH grant R24OD035579 as part of the INCLUDE Project at the NIH Office of the Director, with additional support from NIH grants U2CHL156291 and U24AG092191, the Linda Crnic Institute for Down Syndrome, the Global Down Syndrome Foundation, the Anna and John J. Sie Foundation, and the Boettcher Foundation.

## Contact

Questions and bug reports: please open an [issue](https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/t21-model-atlas-analysis/issues), or contact the corresponding authors of the associated publication.
