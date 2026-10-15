################################################
# Title: Differential expression analysis of Dp16 mouse bulk RNA-seq versus histopathology scores
# Author(s):
#   - Karen Rossmassler
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# DESeq2 modelling for differential expression in Dp16 mice as a function of composite histopathology score
# within 4-month lung and 4-month kidney of the Trisomy 21 Model Atlas.
# This script uses 4-month lung as an example dataset; the same workflow was applied
# to 4-month kidney (see section 0.3 to adapt to other tissue-timepoints).
# The preferred model is the SVA model (surrogate variables + Composite score); simple and
# multivariable models are included for comparison.
# See README.md for more details
#
### Data type(s):
#   A. T21 Model Atlas mouse meta data
#      - T21_atlas_Dp16_metadata_v1.0.tsv (provided in data/)
#   B. Bulk RNA-seq data (GEO SuperSeries GSE347972; 4-month lung SubSeries)
#      - GSE347804_Lung_4mo_counts_GEO.txt.gz (obtain from GEO)
#      - GSE347804_Lung_4mo_RPKMs_GEO.txt.gz (obtain from GEO)
#   C. Gene annotation derived from Gencode vM24 basic GTF
#      - gene_annotation_Gencode.vM24.basic.txt.gz (provided in data/)
#   D. MSigDB mouse Hallmark gene sets v2022.1.Mm
#      - mh.all.v2022.1.Mm.symbols.gmt (provided in data/)
#   E. Dp16 histopathology scores for 4-month lung and kidney, matched to RNA-seq samples
#      - T21_atlas_Dp16_histo_scores_4mo_lung-kidney_v1.0.tsv (provided in data/)

### Workflow:
#   Step 1  - Read in meta, counts, rpkms data for all samples and QC
#   Step 2  - Filter by minimum cpm
#   Step 3  - Comparison groups and covariates setup
#   Step 4  - Generate DESeqDataSet(s) and assess models
#   Step 5  - Run DESeq2 analysis
#   Step 6  - QC checks and overall sample groupings
#   Step 7  - Get results for comparison(s) of interest
#   Step 8  - Export results
#   Step 9  - Extra plots
#   Step 10 - Volcano plot(s)
#   Step 11 - Individual gene plots
#   Step 12 - Heatmap
#   Step 13 - GSEA

### Change Log:
# v1.0
# Initial version
#


# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
# if you get 'project is out-of-sync' warning from renv, try running renv::install() which will often correct the issue.
# OR to install the exact package versions used (requires matching R version):
# renv::restore()
#

## 0.1 Load required libraries ----
library("GEOquery") # for direct GEO downloads
library("DESeq2") # differential expression analysis
library("edgeR") # for cpm() function (can also be used for differential expression analysis)
library("limma") # for removeBatchEffect() function (can also be used for differential expression analysis)
library("BiocParallel") # enables multi-cpu for some of DEseq2 functions
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores)) # enables multi-cpu for some of DEseq2 functions; MulticoreParam is not available on Windows (use SnowParam)
library("apeglm") # used with DESeq2 to 'moderate' fold-changes
library("sva") # for Surrogate Vector Analysis
library("biobroom") # for tidy_sva()
library("openxlsx") # for exporting results as Excel workbooks
library("dendextend") # used for coloring dendograms
library("tidyverse") # required for ggplot2, dplyr etc
library("ggforce") # used for sina plots
library("ggrastr") # required for rasterizing some layers of plots
library("ggrepel") # required for using geom_text and geom_text_repel() to make sample labels for PCA plot
library("RColorBrewer") # color palettes
library("circlize") # color scale generation
library("tidyHeatmap") # tidy heatmaps
library("patchwork") # combining plots
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::select,
  dplyr::count,
  dplyr::rename,
  base::paste,
  matrixStats::rowVars
)
library("here")
source(here("helper_functions_DESeq.R")) # load helper functions
#

## 0.2 Set required parameters ----
# Input data files
gene_anno_file <- here("data", "gene_annotation_Gencode.vM24.basic.txt.gz") # MOUSE; provided with repository
counts_file <- here("data", "GSE347804_Lung_4mo_counts_GEO.txt.gz") # obtain from GEO
rpkms_file <- here("data", "GSE347804_Lung_4mo_RPKMs_GEO.txt.gz") # obtain from GEO
meta_data_file <- here("data", "T21_atlas_Dp16_metadata_v1.0.tsv") # provided with repository
hallmarks_file <- here("data", "mh.all.v2022.1.Mm.symbols.gmt") # MSigDB mouse Hallmarks v2022.1.Mm; provided with repository
histo_file <- here("data", "T21_atlas_Dp16_histo_scores_4mo_lung-kidney_v1.0.tsv") # provided with repository
#
# Settings
min_cpm <- 0.5 # used for low count filtering; default is 0.5
min_samples <- "auto" # used for low count filtering; use a number, "all", or "auto" (sets to half number of samples)
standard_colors <- c("WT" = "grey60", "Dp16" = "#009b4e") # these should be named
#
# Triplicated region - used to identify triplicated genes (chromosome AND coordinates, not chr16 alone)
# Dp16-SPECIFIC: Dp(16)1Yey segmental duplication of MMU16 (GRCm38 coordinates).
# Must be changed for other mouse models (e.g. Ts65Dn, TcMAC21); for human data use chr21.
trip_chr <- "chr16"
trip_start <- 75540514
trip_end <- 97962622
#
out_file_prefix <- "integrated_mouse_Dp16_bulk-histo_4mo-lung_" # used in output file names; REQUIRES CUSTOMIZING for each dataset/comparison
# End required parameters ###
#
## 0.3 Adapting to other tissue-timepoints ----
# This workflow was applied separately to the 4-month kidney in the Atlas.
# To run another tissue-timepoint, edit the following and re-run the script:
#   - section 0.2: counts_file, rpkms_file, out_file_prefix
#   - section 1.1: GEO SubSeries (4-month kidney: GSE347802; see GEO_SubSeries in meta data)
#   - section 1.3: the Tissue filter for the histopathology scores
#   - section 1.4: the Tissue and Age_Group filter, and any excluded samples
#   - section 3 / 4.2: covariates available at that timepoint (drop covariates that are missing or constant)
#   - plot titles ("4-month lung") and the GSEA pathway(s) to plot
#
# Alternatively, these values can be collected in one place, e.g. for 4-month kidney:
# tissue_timepoint <- list(
#   tissue = "Kidney",          # Tissue value in meta data and histopathology scores
#   timepoint = "4mo",          # Age_Group value in meta data
#   geo = "GSE347802",          # GEO SubSeries
#   label = "4-month kidney",   # used in plot titles
#   covariates = c("Sex", "RIN", "Age_in_days")
# )
# out_file_prefix <- paste0("integrated_mouse_Dp16_bulk-histo_", tissue_timepoint$timepoint, "-", str_to_lower(tissue_timepoint$tissue), "_")
# counts_file <- here("data", paste0(tissue_timepoint$geo, "_", tissue_timepoint$tissue, "_", tissue_timepoint$timepoint, "_counts_GEO.txt.gz"))
# rpkms_file <- here("data", paste0(tissue_timepoint$geo, "_", tissue_timepoint$tissue, "_", tissue_timepoint$timepoint, "_RPKMs_GEO.txt.gz"))
# histo <- histo_all |> filter(Tissue == tissue_timepoint$tissue, Score_name == "Composite_Score")
# meta_data <- meta_data |> filter(Tissue == tissue_timepoint$tissue, Age_Group == tissue_timepoint$timepoint)
# multivar_formula <- as.formula(paste("~", paste(tissue_timepoint$covariates, collapse = " + "), "+ Score_value"))
# and titles via paste0(tissue_timepoint$label, ": composite score")
#


# 1. Read in and inspect data ----
#

## 1.0 Load gene names and other annotation ----
# Derived from the Gencode vM24 basic annotation GTF; provided in the data/ directory of this repository
gene_anno <- gene_anno_file %>%
  read_tsv()
gene_anno
#
## 1.1 Download counts and RPKMs from GEO ----
# Requires the GEOquery package. Files are saved to data/ using their GEO file names;
# update counts_file and rpkms_file above to match.
GEOquery::getGEOSuppFiles(
  GEO = "GSE347804", # 4-month lung SubSeries of GSE347972 (see GEO_SubSeries in meta data for other tissue-timepoints)
  makeDirectory = FALSE,
  baseDir = here("data"),
  filter_regex = "counts|RPKM"
)
#
## 1.2 Read in counts and rpkms data ----
counts_data <- counts_file |>
  read_tsv()
counts_data
#
rpkms_data <- rpkms_file |>
  read_tsv()
rpkms_data
#

## 1.3 Read in and inspect histopathology scores ----
# provided in the data/ directory of this repository
# See T21_atlas_Dp16_histo_scores_v1.0_dictionary.tsv for column descriptions
histo_all <- histo_file |>
  read_tsv(na = c("", "NA")) # "NA" in Notes is read as a missing value
histo_all
histo_all |> count(Tissue, Score_name)
# Check for sample notes (e.g. physical outliers) - review before analysis and
# exclude affected samples in section 1.4 if required
histo_all |>
  filter(!is.na(Notes)) |>
  distinct(Sampleid, Tissue, Notes)
#
# Select composite scores for this tissue
histo <- histo_all |>
  filter(
    Tissue == "Lung",
    Score_name == "Composite_Score"
  )
histo
histo |> count(Score_value) # distribution of composite scores
#

## 1.4 Read in meta data ----
# provided in the data/ directory of this repository
# See T21_atlas_Dp16_metadata_v1.0_dictionary.tsv for column descriptions
# NOTE: samples are labelled by Sampleid INSIDE the counts and RPKM files, so
# Sampleid is the key for joining meta data with these data.
# Sampleid should not be reconstructed from other columns (e.g. Tissue), as the tissue
# label in Sampleid differs from Tissue for some samples (e.g. 4mo Whole_Blood, Lgint)
meta_data <- meta_data_file |>
  read_tsv() |>
  mutate(
    Harvest_Date = as.Date(Harvest_Date, format = "%m/%d/%y"),
    DOB = as.Date(DOB, format = "%m/%d/%y"),
    Age_in_days = as.numeric(Harvest_Date - DOB)
  ) |>
  # Select samples for this tissue-timepoint, run within a genotype
  filter(
    Genotype == "Dp16", 
    Tissue == "Lung",
    Age_Group == "4mo"
  ) |>
  # Sample(s) excluded from analysis
  filter(!Sampleid %in% c("T3_510091_Kidney_4mo")) |> # kidney physical outlier: histology section differed from the kidney used for RNA-seq (see Notes check in section 1.3)
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")) # convert to factor and set order (Female = reference)
  )
# Report samples in meta data that have no composite score (dropped by the join below)
meta_data |> anti_join(histo, by = join_by(Sampleid, Tissue))
# and join composite scores to meta data (keeps only samples with a score)
meta_data <- meta_data |>
  inner_join(histo, by = join_by(Sampleid, Tissue)) |>
  select(Sampleid, Genotype, everything())
# Report samples in meta data that are missing from the counts data (e.g. excluded due to QC)
meta_data |> filter(!Sampleid %in% colnames(counts_data))
# and remove them from meta data
meta_data <- meta_data |> filter(Sampleid %in% colnames(counts_data))
# inspect
meta_data
meta_data |> skimr::skim()
meta_data |> count(Genotype)
meta_data |> count(Genotype, Sex)
#

## 1.5 Convert unfiltered counts to tidy/long format and join with meta data ----
# (one row per observation, one column per variable)
#
# Subselect counts data and rpkms to match meta_data by Sampleid
# (ensures same samples and order as in meta_data)
counts_data <- counts_data %>% select(Geneid, any_of(meta_data %>% pull(Sampleid)))
rpkms_data <- rpkms_data %>% select(Geneid, any_of(meta_data %>% pull(Sampleid)))
#
# Gather to long format
counts_data_long <- counts_data %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "raw_count")
counts_data_long
#
rpkms_data_long <- rpkms_data %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM")
rpkms_data_long
#
# Join with metadata by Sampleid
counts_data_long <- counts_data_long %>%
  inner_join(meta_data, by = "Sampleid")
# Check number of samples before and after join
n_samples <- counts_data %>% select(-Geneid) %>% colnames() %>% length()
n_samples_join <- counts_data_long %>% distinct(Sampleid) %>% nrow()
cat("Number of samples in counts_data:", n_samples,
    "\nNumber of samples after join with meta data:", n_samples_join)
if(n_samples != n_samples_join) warning("Possible problem joining counts_data_long with meta data - different number of samples after join")
#
## 1.6 Check raw read count distributions across samples (not normalized) ----
counts_data_long %>%
  ggplot(aes(Sampleid, log2(raw_count + 0.1), color = Genotype)) +
  geom_sina(size = 0.01) +
  scale_color_manual(values = standard_colors) +
  labs(title = "Raw read count distributions across samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_counts_unfiltered.png")), width = 10, height = 5, units = "in")
#

# 2. Filter/remove genes with low expression ----
#
## 2.1 Summary of total reads per sample ----
counts_data %>%
  select(-Geneid) %>%
  colSums() %>%
  summary()
counts_data %>%
  select(-Geneid) %>%
  colSums() %>%
  enframe(name = "Sampleid", value = "Total_reads") |>
  arrange(Total_reads)
#

## 2.2 Filter by minimum counts per million ----
# Keep only rows (transcripts / genes) with greater than `min_cpm` cpm in `min_samples`
# NOTE: 10 counts=0.5 cpm for 20 million reads,  15 counts=0.5 cpm for 30 million reads...
# from Michael Love: https://support.bioconductor.org/p/95840/
# The independent filtering is designed only to filter out low count genes to
# the extent that they are not enriched with small p-values. Here the problem is
# not independent filtering, but that these two genes get a small p-value rather
# than being filtered or having an insignificant p-value. Datasets can be
# different in many ways, and for whatever reason, these two genes survive the
# filtering and get a counterintuitive small p-value. I'd recommend you just use
# a more strict filter in the very beginning, e.g. at least three samples with
# counts greater than 10
#
# Check min_samples and calculate if needed (ie if number is not supplied)
if (min_samples == "all") {
  min_samples=ncol(counts_data) - 1
} else if (min_samples == "auto") {
  min_samples=(ncol(counts_data) - 1) / 2
}
before <- counts_data %>%
  transmute(
    Geneid = Geneid,
    row_sum = rowSums(select(., -Geneid))
  ) %>%
  filter(row_sum > 0) %>%
  nrow()
cpm_data <- counts_data %>%
  column_to_rownames("Geneid") %>%
  cpm()
keep <- cpm_data |>
  as_tibble(rownames = "Geneid") |>
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "cpm") |>
  mutate(cpm > min_cpm) |> # check against min_cpm
  filter(`cpm > min_cpm` == TRUE) |> # and filter
  dplyr::count(Geneid) |> # count samples remaining per Geneid
  filter(n >= min_samples) # filter against min_samples
counts_filtered <- counts_data |>
  filter(Geneid %in% keep$Geneid)
#

## 2.3 Summarize rows before and after filtering ----
cat("Total number of rows:", counts_data %>% nrow(),
    "\nNumber of rows with non-zero total read counts before filtering:", before,
    "\nNumber of rows after filtering:", counts_filtered %>% nrow(), "\n")
# Gather to long format and join with meta data
counts_filtered_long <- counts_filtered %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "raw_count") |>
  inner_join(meta_data)
#
## 2.4 Check filtered read count distributions for each sample (not normalized) ----
counts_filtered_long %>%
  ggplot(aes(Sampleid, log2(raw_count + 0.1), color = Genotype)) +
  geom_sina(size = 0.01) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  labs(title = "CPM-filtered read count distributions across samples", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_counts_filtered.png")), width = 10, height = 5, units = "in")
#

# 3. Groups and/or Covariates setup ----
groups <- meta_data %>% 
  # CUSTOMIZE
  select(Sampleid, Genotype, Sex, Cage_ID, RIN, Age_in_days, Score_value)  %>%  
  mutate(
    Cage_ID = as.character(Cage_ID),
    RIN = (RIN - mean(RIN)) / sd(RIN), # convert to Z-scores to scale and center
    Age_in_days = (Age_in_days - mean(Age_in_days)) / sd(Age_in_days), # convert to Z-scores to scale and center
  )
#

# 4. Generate DESeqDataSet object(s) ----
## 4.1 simple model ----
simple_formula <- as.formula(paste0("~", "Score_value"))
#
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>%
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = simple_formula
)
#
# Check meta data read in to DESeqDataSet (simple version):
colData(dds)
#
## 4.2 multivariable model ----
# Put the variable of interest at the end of the formula. Thus the results
# function will by default pull the condition results unless contrast or name
# arguments are specified.
# see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#
multivar_formula <- as.formula(paste0("~ ", "Sex + RIN + Age_in_days +", "Score_value"))
#
dds_multi <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>%
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = multivar_formula
)
# Check meta data read in to DESeqDataSet (multivariable version):
colData(dds_multi)
#

## 4.3 Likelihood ratio test for multivariable model ----
# For how many genes does adding additional terms give `better` fit?
dds_multi_lrt <- DESeq(dds_multi, parallel = TRUE, test = "LRT", reduced = simple_formula)
#
dds_multi_lrt %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
dds_multi_lrt %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1)
#
# LRT and AIC (section 5.4) are based on the same per-gene deviance difference between
# nested models: the LRT asks in how many genes the added terms significantly improve fit
# (padj < 0.1), while AIC penalizes added terms and is used to choose a single model for
# all genes. Note that SVs are estimated from the same data, so both may overstate the
# improvement from SVA.

## 4.4 SVA model -------
# The goal of the sva is to remove all unwanted sources of variation while
# protecting the contrasts due to the primary variables included in mod. This
# leads to the identification of features that are consistently different
# between groups, removing all common sources of latent variation.
dds_sf <- estimateSizeFactors(dds) # size factors needed for normalized counts
dat  <- counts(dds_sf, normalized = TRUE)
# best to ensure genes with low counts have already been removed
mod  <- model.matrix(~ Score_value, colData(dds_sf)) # designate variable of interest to "protect"
mod0 <- model.matrix(~ 1, colData(dds_sf))
#
svseq <- svaseq(dat, mod, mod0, n.sv = NULL) # if NULL, number of factors will be estimated for you
#
groups_sva <- svseq %>%
  # extract + add additional coefficients to the estimated surrogate variables to allow comparison
  biobroom::tidy_sva(addVar = colData(dds_sf)) # includes sizeFactor
groups_sva
#
### Compare SVs with known categorical variables ----
groups_sva %>%
  pivot_longer(matches("^sv"), names_to = "SV", values_to = "estimate") %>%
  ggplot(aes(Sex, estimate)) +
  geom_sina() +
  facet_wrap(~SV)
#
### Compare SVs with known continuous variables ----
groups_sva %>%
  pivot_longer(matches("^sv"), names_to = "SV", values_to = "estimate") %>%
  ggplot(aes(sizeFactor, estimate)) +
  geom_point() +
  facet_wrap(~SV)
groups_sva %>%
  pivot_longer(matches("^sv"), names_to = "SV", values_to = "estimate") %>%
  ggplot(aes(RIN, estimate)) +
  geom_point() +
  facet_wrap(~SV)
groups_sva %>%
  pivot_longer(matches("^sv"), names_to = "SV", values_to = "estimate") %>%
  ggplot(aes(Age_in_days, estimate)) +
  geom_point() +
  facet_wrap(~SV)
groups_sva %>%
  pivot_longer(matches("^sv"), names_to = "SV", values_to = "estimate") %>%
  ggplot(aes(Score_value, estimate)) + # Score_value should be preserved ie not correlated to any SV
  geom_point() +
  facet_wrap(~SV)
#
### Build DESeq model with SVs -----
sva_formula <- as.formula(
  paste0("~", groups_sva %>% select(matches("^sv")) %>% colnames() %>% paste(collapse = " + "),
         "+ Score_value"
  ))
sva_formula
#
dds_sva <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>%
    select(Geneid, groups_sva %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups_sva,
  design = sva_formula
)
# Check metadata read in to DESeqDataSet (sva version):
colData(dds_sva)
#
### Likelihood ratio test for SVA model ----
# For how many genes does adding surrogate vectors give `better` fit?
dds_sva_lrt <- DESeq(dds_sva, parallel = TRUE, test = "LRT", reduced = simple_formula)
dds_sva_lrt %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
dds_sva_lrt %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1)
#
# 5. Run DESeq2 analysis ----
# Default analysis runs the following steps:
# 1. estimation of size factors
# 2. estimation of dispersion
# 3. Negative Binomial GLM fitting and Wald statistics
# note: since ~v1.16, shrinkage of log2foldChange is not run by default - this has moved to lfcShrink function
#

## 5.1 Run simple model ----
dds <- DESeq(dds, parallel=TRUE)
#

## 5.2 Run multivariable model ----
dds_multi <- DESeq(dds_multi, parallel=TRUE)
#

## 5.3 Run SVA model ----
dds_sva <- DESeq(dds_sva, parallel=TRUE)
#
## 5.4 Compare models by AIC ----
# Gene-wise AIC from model deviance (deviance = -2 * log-likelihood)
aic_dds <- get_aic(dds)
aic_dds_multi <- get_aic(dds_multi)
aic_dds_sva <- get_aic(dds_sva)
#
# Compare simple model with multivariable and SVA models
# Here, a positive AIC_diff for a gene = better fit for the more complex model (dashed lines at +/-10)
# A single parsimonious model is chosen and applied to all genes, rather than a
# per-gene model: the more complex model is preferred if it improves fit (lower AIC)
# for the majority of genes, accepting that fit may be worse for some genes by this metric
aic_simple_vs_multi <- plot_aic_diff(
  aic_dds, aic_dds_multi,
  subtitle = paste(deparse1(simple_formula), "vs.", deparse1(multivar_formula))
)
aic_simple_vs_sva <- plot_aic_diff(
  aic_dds, aic_dds_sva,
  subtitle = paste(deparse1(simple_formula), "vs.", deparse1(sva_formula))
)
#
aic_simple_vs_multi + aic_simple_vs_sva
ggsave(filename = here("plots", paste0(out_file_prefix, "AIC_comparison.png")), width = 13, height = 5, units = "in", bg = "white")
#
# 6. QC checks and overall sample grouping(s) ----
#
## 6.1 Check Size Factors used for normalization ----
get_size_fcts(dds) |>
  arrange(SizeFactor) |>
  summary()
#

## 6.2 Get normalized counts ----
nc <- dds %>% counts(normalized = TRUE) %>% as_tibble(rownames="Geneid")
nc
#

## 6.3 Get vst transformed values ----
# Variance-stablizing transformation and normalization ± covariate correction
vst_mat <- assay(vst(dds))
vst_mat |> as_tibble(rownames = "Geneid")
#
### Generate sv-adjusted vst values ----
vst_mat_sva_adj <-
  vst_mat %>%
  limma::removeBatchEffect(
    covariates = colData(dds_sva) %>% as_tibble() %>% select(matches("^sv")),
    design = colData(dds_sva) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
  )
vst_mat_sva_adj |> as_tibble(rownames = "Geneid")
#
### Distributions of vst values ----
vst_mat %>%
  as_tibble(rownames = "Geneid") %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "VST") %>%
  inner_join(meta_data) %>%
  ggplot(aes(Sampleid, VST, color = Genotype)) +
  geom_sina(size = 0.01) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  labs(title = "Unadjusted VST-normalized count distributions across samples", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
## 6.4 PCA plot(s) of normalized and variance-transformed count data ± covariate correction ----
#
### Overall groupings by PCA ----
plotPCA_custom2(vst_mat,
                PCA_by = "variance",
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title = "PCA plot",
                subtitle = "unadjusted",
                color_var = "Score_value",
                shapes = "Sex",
                labels = FALSE,
                x_lower_lim = -110, # CUSTOMIZE LIMITS
                x_upper_lim = 110,
                y_lower_lim = -110,
                y_upper_lim = 110)
#
plotPCA_custom2(x = vst_mat_sva_adj,
                PCA_by = "variance",
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title="PCA plot",
                subtitle = "SVA-adjusted",
                color_var = "Score_value",
                shapes = "Sex",
                labels = FALSE,
                x_lower_lim = -100, # CUSTOMIZE LIMITS
                x_upper_lim = 100,
                y_lower_lim = -100,
                y_upper_lim = 100)
#
unadjusted_pca_plot + `SVA-adjusted_pca_plot`
ggsave(filename = here("plots", paste0(out_file_prefix, "PCA_unadjusted_vs_SVA.png")), width = 13, height = 5, units = "in", bg = "white")
#
### Check meta data associations with PCs -------
# Unadjusted
c("Sex", "RIN", "Age_in_days", "Score_value") %>%
  paste("loadings ~", .) |>
  set_names() %>%
  map(., ~ pca_lm_function(unadjusted_PC_loadings, unadjusted_PC_percVar, .))
# SVA-adjusted
c("Sex", "RIN", "Age_in_days", "Score_value") %>%
  paste("loadings ~", .) |>
  set_names() %>%
  map(., ~ pca_lm_function(`SVA-adjusted_PC_loadings`, `SVA-adjusted_PC_percVar`, .))
#
### Verify sex groupings by PCA ----
plotPCA_custom2(x = vst_mat[rownames(vst_mat) %in% (gene_anno %>% filter(chr == "chrY" | chr == "chrX"))$Geneid,],
                PCA_by = "variance",
                save_PCs = FALSE,
                plot_title="PCA plot - chrY+X",
                subtitle = "unadjusted",
                color_var = "Sex",
                shapes = "Sex",
                labels = FALSE,
                x_lower_lim = -20,
                x_upper_lim = 20,
                y_lower_lim = -20,
                y_upper_lim = 20)
ggsave(filename = here("plots", paste0(out_file_prefix, "PCA_sex_genes.png")), width = 5, height = 5, units = "in", bg = "white")
#
# 7. Get DESeq2 results ----
# see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# To see available results: dds %>% resultsNames() list of results
# available does not show all combinations - only comparisons against reference
# level of last variable in the design formula are available by default
#
# The results function without any arguments will automatically perform a
# contrast of the last level of the last variable in the design formula over the
# first level; we specify the comparison of interest when calling results()
#
## NOTE cpm filtering does not make much difference with independent filtering
## on, but will get rid of some odd cases that get a small p-value despite low
## counts and may speed up analysis time
#

## 7.1 Define comparisons of interest ----
# NAMES VERSION (need to match elements in resultsNames())
# Score_value is a continuous predictor, so its effect is a single slope coefficient
# (log2 fold-change per one-unit increase in composite score) rather than a comparison
# between factor levels. Results are therefore extracted by coefficient name with
# get_results_tbl_names(); the contrasts version (get_results_tbl()) requires a factor
# with numerator and denominator levels to relevel and to build the coefficient name used
# for apeglm shrinkage, which only accepts a coefficient (coef) rather than a contrast.
comparisons <- list("Score_value")

## 7.2 Results summaries ----
#
for (i in comparisons) { # NAMES VERSION
  dds %>% get_results_sum_names(i, show_ind_filt_off=FALSE)
}
for (i in comparisons) { # NAMES VERSION
  dds_multi %>% get_results_sum_names(i, show_ind_filt_off=FALSE)
}
for (i in comparisons) { # NAMES VERSION
  dds_sva %>% get_results_sum_names(i, show_ind_filt_off=FALSE)
}
#
## 7.3 Assemble DESeq2 results table(s) ----
# NAMES VERSION (spaces or dashes in variables will cause problems here)
# Score_value is a continuous predictor, so its effect is a single slope coefficient
# (log2 fold-change per one-unit increase in composite score) rather than a comparison
# between factor levels. Results are therefore extracted by coefficient name with
# get_results_tbl_names(); the contrasts version (get_results_tbl()) requires a factor
# with numerator and denominator levels to relevel and to build the coefficient name used
# for apeglm shrinkage, which only accepts a coefficient (coef) rather than a contrast.
#
# Initialize empty vector to store names of results objects for later reference
comparisons_results <- character()
#
### 7.3.1 Simple model ------
for (comparison in comparisons) {
  name_res <- paste("res_simple", comparison, sep = "_") # CUSTOMIZE name of results table
  comparisons_results <- c(comparisons_results, name_res)
  res_temp <- dds %>% # SET to correct dds object
    get_results_tbl_names(
      name = comparison,
      shrink_type = "apeglm"
    )
  res_temp %>% assign(name_res, ., pos=1)
}
### 7.3.2 Multivariable model(s) ----
for (comparison in comparisons) {
  name_res <- paste("res_multi", comparison, sep = "_")  # CUSTOMIZE name of results table
  comparisons_results <- c(comparisons_results, name_res)
  res_temp <- dds_multi %>% # SET to correct dds object
    get_results_tbl_names(
      name = comparison,
      shrink_type = "apeglm"
    )
  res_temp %>% assign(name_res, ., pos=1)
}
### 7.3.3 SVA model(s) ----
for (comparison in comparisons) {
  name_res <- paste("res_sva", comparison, sep = "_") # CUSTOMIZE name of results table
  comparisons_results <- c(comparisons_results, name_res)
  res_temp <- dds_sva %>% # SET to correct dds object
    get_results_tbl_names(
      name = comparison,
      shrink_type = "apeglm"
    )
  res_temp %>% assign(name_res, ., pos=1)
}
#
# List of results tables:
comparisons_results
# Preview results tables:
res_sva_Score_value

## 7.4 Standard volcano plot(s) - no gene labels ----
res_sva_Score_value |>
  plotVolcano(
    sig = 0.1,
    title = "4-month lung: composite score",
    subtitle = paste("Model:", paste(sva_formula, collapse = ""))
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sva_Score_value", "_Volcano.png")), width = 5, height = 5, units = "in")
#
# 8. Export results ----
#
# Usually only export preferred model results, in this case:
comparisons_results |> str_subset("res_sva") |> export_res()
# To export all results tables:
# comparisons_results |> export_res()
#
# 9. Extra plots ----
#
## 9.1 Plot dispersion estimates ----
#
dds_sva |> ggplotDispEsts()
#

## 9.2 Plot p-value distributons ----
#
res_sva_Score_value |> plotPvals()
#
# 10. labelled volcano plots ----
#
## 10.1 label top significant genes ----
#
# get y limits
max(
  res_sva_Score_value %>% pull(padj) %>% -log10(.) %>% max(na.rm = TRUE)
)
# get x limits
max(
  res_sva_Score_value %>% pull(log2FoldChange_adj) %>% abs(.) %>% max(na.rm = TRUE) %>% ceiling()
)
#
v1 <- res_sva_Score_value %>%
  volcano_plot_lab(
    title = "4-month lung: Composite score (SVA)",
    subtitle = paste0(
      paste("Model:", paste(sva_formula, collapse = ""), "\n"),
      "[Down: ", (.) %>% filter(padj < 0.1 & FoldChange_adj <1) %>% nrow(), "; Up: ", (.) %>% filter(padj < 0.1 & FoldChange_adj >1) %>% nrow(), "]"
    ),
    labels = TRUE,
    n_labels = 3,
    raster = TRUE
  )
#
v1 <- v1 +
  coord_cartesian(xlim = c(-1, 1) * 2, ylim = c(0, 4)) # CUSTOMIZE limits using the values above
v1
ggsave(v1, filename = here("plots", paste0(out_file_prefix, "sva_Score_value_Volcano_labelled", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#


# 11. Individual gene plots ----
#
## Generate SV-adjusted RPKM values ----
# Can plot normalized counts or raw RPKMs but usually prefer RPKMs adjusted with same covariables as preferred model
rpkm_sva_adj <- rpkms_data %>%
  select(Geneid, (colData(dds_sva) %>% as_tibble() %>% pull(Sampleid))) %>%
  mutate_at(2:ncol(.), ~ log2(.)) %>%  # log2 transformation of RPKMs
  column_to_rownames(var = "Geneid") %>% # convert to dataframe to preserve Geneid
  limma::removeBatchEffect(
    covariates = colData(dds_sva) %>% as_tibble() %>% select(matches("^sv")),
    design = colData(dds_sva) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
  ) %>%
  as_tibble(rownames="Geneid") %>% # convert back to tibble
  mutate_at(2:ncol(.), ~ 2^(.)) # remove log2 transformation
# # Export adjusted RPKMs
# rpkms_data_long %>%
#   inner_join(
#     rpkm_sva_adj %>%
#       pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM_adj")
#   ) %>%
#   mutate(adjustment = "SVA") %>%
#   inner_join(gene_anno %>% select(Geneid, gene_name, chr, gene_type)) %>%
#   write_tsv(file = here("results", paste0(out_file_prefix, "RPKMs_sva_adj.txt.gz")))
# #
#
## Get genes of interest ----
# usually run only for preferred model
top_signif <- res_sva_Score_value %>%
  filter(padj < 0.1) %>%
  slice_min(pvalue, n = 10) %>%
  arrange(-FoldChange_adj) %>%
  slice_head(n = 10) %>% # ensure only 10 as slice_min() will sometimes have ties
  select(Geneid, Gene_name)
#
## Scatter plots ----
s1 <- rpkm_sva_adj %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM") %>% 
  inner_join(meta_data) %>%
  inner_join(top_signif) %>%
  mutate(Gene_name = fct_relevel(Gene_name, top_signif %>% pull(Gene_name))) %>%
  group_by(Geneid) %>%
  mutate(extreme = rstatix::is_extreme(log2(RPKM))) %>%
  filter(!extreme) %>%
  ungroup() %>%
  ggplot(aes(x = Score_value, y = log2(RPKM), color = Score_value)) +  # continuous gradient
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ Gene_name, scales = "free", nrow = 2) +
  scale_color_gradientn(
    colors = c("#bae4b3", "#2e7d32", "#004d1f"),
    # values = scales::rescale(c(0, 3, 4))
    ) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(
    title = "4mo Lung composite score: top significant genes",
    subtitle = "SVA-adjusted; Extreme outliers removed",
    x = "Score",
    y = "log2(RPKM)",
    color = "Score")
s1
ggsave(s1, filename = here("plots", paste0(out_file_prefix, "top_signif_RPKM_scatter", ".png")), width = 15, height = 5, units = "in")
# # May want to rasterize before saving as pdf
# s1r <- ggrastr::rasterize(s1, layers='Point', dpi = 600, dev = "ragg_png")
# ggsave(s1r, filename = here("plots", paste0(out_file_prefix, "top_signif_RPKM_scatter", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#
# 12. Heatmaps ----
#
## Generate gene-wise Z-scores ----
# usually need to transform in some way to compensate for wide range of expression levels
zscores_sva_adj <- rpkm_sva_adj %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM") |>
  group_by(Geneid) |>
  mutate(
    zscore = (log2(RPKM) - mean(log2(RPKM), na.rm = TRUE)) / sd(log2(RPKM), na.rm = TRUE)
  ) |>
  ungroup()
#
## Plot genes of interest as heatmap
# collect data for heatmap
hm_dat <- zscores_sva_adj |>
  inner_join(gene_anno) |>
  inner_join(meta_data) %>%
  inner_join(top_signif) %>%
  mutate(
    gene_name = fct_relevel(Gene_name, top_signif %>% pull(Gene_name)), # control row plotting order
    Sampleid = fct_reorder(Sampleid, Score_value) # order columns by composite score
  )
# Generate centered color scale
hm_lim <- hm_dat |>
  pull(zscore) |>
  abs() |>
  max() |>
  round(2)
breaks <- seq(-hm_lim, hm_lim, length.out = 11)
hm_palette <- circlize::colorRamp2(
  breaks,
  RColorBrewer::brewer.pal(11, "RdBu") |> rev()
)
# plot heatmap
hm_dat |>
  tidyHeatmap::heatmap(
    gene_name, # be careful of non-unique
    Sampleid,
    zscore,
    palette_value = hm_palette,
    heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_lim, hm_lim, length.out = 5)),
    cluster_rows = FALSE,
    cluster_columns = FALSE, # keep column order by composite score
    row_title = NULL,
    show_column_names = FALSE,
    column_title = NULL,
    border = TRUE
  ) |>
  annotation_tile(
    Score_value,
    palette = circlize::colorRamp2(
      seq(min(hm_dat$Score_value), max(hm_dat$Score_value), length.out = 3),
      c("#bae4b3", "#2e7d32", "#004d1f") # same colors as scatter plots
    )
  ) |>
  wrap_heatmap() +
  labs(
    title = "4-month lung: top significant genes",
    subtitle = "SVA-adjusted Z-scores",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "top_signif_zscore_heatmap", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")
#
# 13. GSEA Hallmarks analysis ----
#
# Mouse
hallmarks <- hallmarks_file %>%
  fgsea::gmtPathways(gmt.file = .)
#
## Generate ranks ----
ranks_sva_Score_value <- res_sva_Score_value %>%
  filter(!is.na(log2FoldChange_adj)) %>% # need to remove NA rows that will break plotEnrichment2()
  select(ID = Gene_name, t = log2FoldChange_adj) %>%
  arrange(-abs(t)) %>% # to keep strongest of any duplicates
  distinct(ID, .keep_all = TRUE) %>% # to avoid duplicates
  tibble::deframe() # convert to named numerical vector
#
## Run fgsea ----
# Unweighted enrichment statistic, as used for the Atlas
set.seed(1234) # fgseaMultilevel is stochastic; set seed for reproducible p-values
hallmarks_sva_Score_value <- run_fgsea2(geneset = hallmarks, ranks = ranks_sva_Score_value, weighted = FALSE)
hallmarks_sva_Score_value

## Export GSEA results ----
list(
  "Unweighted" = hallmarks_sva_Score_value %>% unnest(leadingEdge) %>% group_by(pathway, pval, padj, log2err, ES, NES, size) %>% summarize(leadingEdge = paste(leadingEdge, collapse = ",")) %>% arrange(padj, -abs(NES))
) |>
  export_excel(filename = "GSEA_Hallmarks")
#
## GSEA barplot(s) ----
hallmarks_sva_Score_value %>%
  filter(padj < 0.1) %>%
  slice_max(order_by = abs(NES), n = 20) %>%
  arrange(NES) %>%
  mutate(
    pathway = str_remove(pathway, "^HALLMARK_") %>% str_replace_all("_", " ") %>% str_to_title(),
    pathway = fct_inorder(pathway)
  ) %>%
  ggplot(aes(-log10(padj), pathway, fill = NES)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_col(color = "black") +
  scale_fill_gradient2(
    low = "#542788",
    mid = "#f7f7f7",
    high = "#b35806",
    midpoint = 0,
    guide = "colourbar",
  ) +
  labs(
    title = "4-month lung: Composite Score\nGSEA Hallmarks (Unweighted, Top 20 q < 0.1)",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "barplot_Hallmarks_sig_top20", ".pdf")), device = cairo_pdf, width = 6, height = 5, units = "in")
#
## GSEA enrichment plot(s) ----
plotEnrichment2(
  pathway = hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  stats = ranks_sva_Score_value,
  res = hallmarks_sva_Score_value,
  title = "4-month lung: Composite Score\nInterferon Gamma Response (Unweighted)"
)
ggsave(filename = here("plots", paste0(out_file_prefix, "enrichment_IFNg_response_unweighted", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#

# Save workspace ----
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# To reload:
# load(here("rdata", paste0(out_file_prefix, ".RData")))
#


# Session info ----
writeLines(capture.output(sessionInfo()), con = here("results", paste0(out_file_prefix, "session_info.txt")))
#