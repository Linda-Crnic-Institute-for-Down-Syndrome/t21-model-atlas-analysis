################################################
# Title: Dp16 mouse histopathology scores at 1 and 4 months: genotype comparisons
# Author(s):
#   - Srija Chillamcherla
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:
# Compare tissue-specific composite and component histopathology scores between
# Dp16 and wild-type (WT) mice at 1 and 4 months in the Trisomy 21 Model Atlas.
# Each tissue-timepoint is tested separately using one-sided, unpaired Wilcoxon
# rank-sum tests for higher scores in Dp16 mice, with rank-based effect sizes.
# See README.md for more details
#

### Data type(s):
#   A. Dp16 histopathology scores (provided in data/)
#      - T21_atlas_Dp16_histo_scores_1mo_v1.0.tsv
#      - T21_atlas_Dp16_histo_scores_4mo_v1.0.tsv
#   B. T21 Model Atlas mouse meta data (provided in data/)
#      - T21_atlas_Dp16_metadata_v1.0.tsv
#

### Workflow:
#   Step 1 - Read in and inspect histopathology scores and meta data
#   Step 2 - Composite score Wilcoxon tests and effect sizes
#   Step 3 - Composite score plots and effect size heatmaps
#   Step 4 - Component score Wilcoxon tests and effect sizes
#   Step 5 - Component score plots and effect size heatmaps
#   Step 6 - Selected component plots
#   Step 7 - Export results
#

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
library("tidyverse") # data import, transformation, and plotting
library("rstatix") # Wilcoxon tests and rank-based effect sizes
library("coin") # required by rstatix::wilcox_effsize()
library("tidyHeatmap") # effect size heatmaps
library("RColorBrewer") # heatmap color palettes
library("circlize") # continuous heatmap color mapping
library("openxlsx") # for exporting results as Excel workbooks
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::select,
  dplyr::count,
  dplyr::rename,
  base::paste
)
library("here")
source(here("helper_functions_histo.R")) # load helper functions
#

## 0.2 Set required parameters ----
# Input data files
histo_1mo_file <- here("data", "T21_atlas_Dp16_histo_scores_1mo_v1.0.tsv") # provided with repository
histo_4mo_file <- here("data", "T21_atlas_Dp16_histo_scores_4mo_v1.0.tsv") # provided with repository
meta_data_file <- here("data", "T21_atlas_Dp16_metadata_v1.0.tsv") # provided with repository
#
# Settings
standard_colors <- c("WT" = "grey30", "Dp16" = "#009b4e") # these should be named
#
out_file_prefix <- "histo_mouse_Dp16_tissue-timepoint_" # used in output file names
# End required parameters ###
#


# 1. Read in and inspect data ----
#

## 1.1 Read in and inspect histopathology scores ----
# provided in the data/ directory of this repository
# See T21_atlas_Dp16_histo_scores_v1.0_dictionary.tsv for column descriptions
# Scores are recorded per animal (Animal_No) within each age group
histo_scores <- list("1mo" = histo_1mo_file, "4mo" = histo_4mo_file) |>
  map(\(f) read_tsv(f, na = c("", "NA"))) |> # "NA" in Score_value and Notes is read as a missing value
  list_rbind(names_to = "Age")
histo_scores
histo_scores |> distinct(Age, Animal_No) |> count(Age)
histo_scores |> filter(Score_name == "Composite_Score") |> count(Age, Tissue)
# Check for sample notes (e.g. physical outliers, missing tissue) - review before
# analysis and exclude affected scores in section 1.2 if required
histo_scores |>
  filter(!is.na(Notes)) |>
  distinct(Age, Animal_No, Tissue, Notes) |>
  print(n = Inf)
#

## 1.2 Exclude scores ----
# Score(s) excluded from analysis
excluded_scores <- tribble(
  ~Age,  ~Animal_No, ~Tissue,
  "4mo", "T3_091",   "Kidney" # physical outlier: histology section differed from the kidney used for RNA-seq (see Notes)
)
histo_scores |> semi_join(excluded_scores, by = join_by(Age, Animal_No, Tissue)) # report
histo_scores <- histo_scores |> anti_join(excluded_scores, by = join_by(Age, Animal_No, Tissue))
#

## 1.3 Read in meta data ----
# provided in the data/ directory of this repository
# See T21_atlas_Dp16_metadata_v1.0_dictionary.tsv for column descriptions
# Meta data are per RNA-seq sample; reduce to one row per animal and build the
# histopathology key Animal_No = <Mouse_ID>_<last three digits of Cage_ID>
meta_data <- meta_data_file |>
  read_tsv() |>
  filter(Age_Group %in% c("1mo", "4mo")) |>
  mutate(Animal_No = paste0(Mouse_ID, "_", str_sub(as.character(Cage_ID), -3))) |>
  distinct(Age = Age_Group, Animal_No, Sex, Genotype) |>
  mutate(Genotype = fct_relevel(Genotype, "WT", "Dp16"))
# Animal_No should be unique within each age group
meta_data |> count(Age, Animal_No) |> filter(n > 1)
#

## 1.4 Join scores with meta data ----
# Report scored animals that are missing from the meta data (dropped by the join below)
histo_scores |> distinct(Age, Animal_No) |> anti_join(meta_data, by = join_by(Age, Animal_No))
# Report animals in meta data without histopathology scores
meta_data |> anti_join(histo_scores, by = join_by(Age, Animal_No)) |> count(Age)
# and join meta data to scores (keeps only scored animals present in meta data)
histo_meta_data <- histo_scores |>
  inner_join(meta_data, by = join_by(Age, Animal_No))
# inspect
histo_meta_data |> distinct(Age, Animal_No, Sex, Genotype) |> count(Age, Genotype)
histo_meta_data |> distinct(Age, Animal_No, Sex, Genotype) |> count(Age, Genotype, Sex)
#


# 2. Composite score Wilcoxon tests ----
# One-sided hypothesis: Dp16 histopathology scores exceed WT scores.
# A one-sided (alternative = "greater"), unpaired Wilcoxon rank-sum test is run
# separately for each tissue-timepoint; rstatix ref.group = "Dp16" makes Dp16 group1,
# so "greater" tests Dp16 > WT. A one-sided test requires a prespecified directional
# hypothesis; review before applying to other data.
# Correction within rstatix::wilcox_test() is turned off (p.adjust.method = "none");
# Benjamini-Hochberg adjusted p values (BH_padj) are calculated across tissues within
# each timepoint in section 2.4. Asterisks in plots use nominal p < 0.05, as in the
# original analysis.
#
composite_data <- histo_meta_data |>
  filter(Score_name == "Composite_Score") |>
  filter(!is.na(Score_value), !is.na(Genotype), !is.na(Tissue))
#

## 2.1 Prepare grouped data for testing ----
# One group per tissue-timepoint. Groups are removed if:
#   - either genotype is missing or has fewer than 2 animals
#   - scores are identical in all animals (e.g. all 0): no ranks can differ, so the
#     test and effect size are undefined (at 1 month: Small intestines, Colon, Spleen)
# Removed tissues remain in the individual score plots (section 3.1)
composite_test_data <- composite_data |>
  group_by(Age, Tissue, Score_name) |>
  add_count(Genotype) |>
  filter(n_distinct(Genotype) == 2, !any(n < 2), n_distinct(Score_value) > 1)
# Report removed tissue-timepoints
composite_data |>
  distinct(Age, Tissue) |>
  anti_join(composite_test_data |> ungroup() |> distinct(Age, Tissue), by = join_by(Age, Tissue))
#

## 2.2 Tests ----
wilcox_composite <- composite_test_data |>
  rstatix::wilcox_test(
    formula = Score_value ~ Genotype,
    ref.group = "Dp16", # Dp16 as group1 so that "greater" tests Dp16 > WT
    paired = FALSE,
    alternative = "greater",
    p.adjust.method = "none", # BH correction applied manually below
    detailed = TRUE
  )
wilcox_composite
#

## 2.3 Rank-based effect sizes ----
# r = |Z| / sqrt(N) (rstatix returns the absolute value, so r has no direction)
effsize_composite <- composite_test_data |>
  rstatix::wilcox_effsize(
    formula = Score_value ~ Genotype,
    ref.group = "Dp16",
    paired = FALSE,
    detailed = TRUE
  )
#

## 2.4 Combine tests and effect sizes ----
# signed_r adds direction to r: positive = higher scores in Dp16, negative = higher in WT.
# Direction comes from the rank-sum statistic W for Dp16 (group1) relative to its null
# expectation n1 * n2 / 2, i.e. whether Dp16 animals have higher or lower mean ranks.
composite_results <- wilcox_composite |>
  inner_join(
    effsize_composite |> select(Age, Tissue, Score_name, .y., group1, group2, effsize, magnitude),
    by = join_by(Age, Tissue, Score_name, .y., group1, group2)
  ) |>
  mutate(signed_r = effsize * sign(statistic - n1 * n2 / 2)) |>
  mutate(BH_padj = p.adjust(p, method = "BH"), .by = Age) |> # BH across tissues within each timepoint
  arrange(Age, p)
composite_results
#


# 3. Composite score plots and effect size heatmaps ----
# Asterisks indicate nominal one-sided p < 0.05
#
## 3.1 Individual score plots ----
composite_plot_data <- composite_data |>
  left_join(composite_results |> select(Age, Tissue, Score_name, p), by = join_by(Age, Tissue, Score_name)) |>
  mutate(Tissue_label = if_else(!is.na(p) & p < 0.05, paste0(Tissue, "*"), Tissue))
#
set.seed(1234) # jittered point positions
composite_1mo_plot <- composite_plot_data |>
  filter(Age == "1mo") |>
  individual_score_plot(
    facets = ~ Tissue_label,
    title = "Composite score by genotype across tissues - 1mo (one-sided Wilcoxon test, Dp16 > WT)",
    y_label = "Composite score",
    ncol = 11
  )
composite_1mo_plot
ggsave(filename = here("plots", paste0(out_file_prefix, "composite_1mo_genotype.pdf")),
       plot = composite_1mo_plot, device = cairo_pdf, width = 19, height = 3, units = "in")
#
set.seed(1234) # jittered point positions
composite_4mo_plot <- composite_plot_data |>
  filter(Age == "4mo") |>
  individual_score_plot(
    facets = ~ Tissue_label,
    title = "Composite score by genotype across tissues - 4mo (one-sided Wilcoxon test, Dp16 > WT)",
    y_label = "Composite score",
    ncol = 12
  )
composite_4mo_plot
ggsave(filename = here("plots", paste0(out_file_prefix, "composite_4mo_genotype.pdf")),
       plot = composite_4mo_plot, device = cairo_pdf, width = 21, height = 3, units = "in")
#

## 3.2 Effect size heatmaps ----
composite_heatmap_data <- composite_results |>
  filter(!is.na(effsize), !is.na(p)) |>
  mutate(
    Age_group = factor(Age, levels = c("1mo", "4mo"), labels = c("1 month", "4 month")),
    Tissue = fct_relevel(Tissue, composite_results |> filter(Age == "4mo") |> arrange(-effsize) |> pull(Tissue) |> unique())
  )
#
### Effect size (r, magnitude only) ----
# Color scale set to the data: 0 to the maximum r (rounded up to 2 decimals)
composite_max_r <- composite_heatmap_data |> pull(effsize) |> max() |> (\(x) ceiling(x * 100) / 100)()
composite_heatmap <- composite_heatmap_data |>
  as_tibble() |> # need to strip rstatix_test class to prevent tidyHeatmap error
  tidyHeatmap::heatmap(
    .row = Tissue,
    .column = Age_group,
    .value = effsize,
    scale = "none",
    palette_value = circlize::colorRamp2(
      seq(0, composite_max_r, length.out = 9),
      RColorBrewer::brewer.pal(9, "Reds")
    ),
    heatmap_legend_param = list(
      title = "effect size (r)",
      color_bar = "continuous",
      at = seq(0, composite_max_r, length.out = 5)
    ),
    cluster_columns = FALSE,
    cluster_rows = TRUE, # row order is set by clustering
    border = TRUE,
    column_names_rot = 45
  ) |>
  layer_asterisk(p < 0.05) |>
  tidyHeatmap::wrap_heatmap() +
  labs(title = "Wilcox effect size (unsigned): Dp16 vs. WT")
composite_heatmap
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_wilcox_results_combined.pdf")),
       plot = composite_heatmap, device = cairo_pdf, width = 5, height = 7, units = "in")
#
### Signed effect size (signed r) ----
# Diverging scale centered on 0 and symmetric about the largest |signed r|:
# red = higher scores in Dp16, blue = higher scores in WT
composite_max_signed_r <- composite_heatmap_data |> pull(signed_r) |> abs() |> max() |> (\(x) ceiling(x * 100) / 100)()
composite_heatmap_signed <- composite_heatmap_data |>
  as_tibble() |> # need to strip rstatix_test class to prevent tidyHeatmap error
  tidyHeatmap::heatmap(
    .row = Tissue,
    .column = Age_group,
    .value = signed_r,
    scale = "none",
    palette_value = circlize::colorRamp2(
      seq(-composite_max_signed_r, composite_max_signed_r, length.out = 11),
      RColorBrewer::brewer.pal(11, "RdBu") |> rev()
    ),
    heatmap_legend_param = list(
      title = "signed r",
      color_bar = "continuous",
      at = seq(-composite_max_signed_r, composite_max_signed_r, length.out = 5)
    ),
    cluster_columns = FALSE,
    cluster_rows = TRUE, # row order is set by clustering (may differ from the magnitude heatmap)
    border = TRUE,
    column_names_rot = 45
  ) |>
  layer_asterisk(p < 0.05) |>
  tidyHeatmap::wrap_heatmap() +
  labs(title = "Wilcox effect size (signed): Dp16 vs. WT")
composite_heatmap_signed
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_wilcox_results_combined_signed.pdf")),
       plot = composite_heatmap_signed, device = cairo_pdf, width = 5, height = 7, units = "in")
#


# 4. Component score Wilcoxon tests ----
# Same one-sided (Dp16 > WT), unpaired Wilcoxon rank-sum tests as for composite scores
# (section 2), run separately for each component within each tissue-timepoint.
# Correction within rstatix::wilcox_test() is turned off (p.adjust.method = "none");
# Benjamini-Hochberg adjusted p values (BH_padj) are calculated across components within
# each tissue-timepoint in section 4.4. Asterisks in plots use nominal p < 0.05, as in
# the original analysis.
#
component_data <- histo_meta_data |>
  filter(Score_name != "Composite_Score") |>
  filter(!is.na(Score_value), !is.na(Genotype), !is.na(Tissue))
#

## 4.1 Prepare grouped data for testing ----
# One group per component within each tissue-timepoint. Groups are removed if:
#   - either genotype is missing or has fewer than 2 animals
#   - scores are identical in all animals (e.g. all 0): no ranks can differ, so the
#     test and effect size are undefined
# Removed components are retained in the results with p = NA (section 4.4)
component_test_data <- component_data |>
  group_by(Age, Tissue, Score_name) |>
  add_count(Genotype) |>
  filter(n_distinct(Genotype) == 2, !any(n < 2), n_distinct(Score_value) > 1)
# Number of components tested and removed per timepoint
component_data |>
  distinct(Age, Tissue, Score_name) |>
  left_join(component_test_data |> ungroup() |> distinct(Age, Tissue, Score_name) |> mutate(tested = TRUE),
            by = join_by(Age, Tissue, Score_name)) |>
  count(Age, tested)
#

## 4.2 Tests ----
wilcox_component <- component_test_data |>
  rstatix::wilcox_test(
    formula = Score_value ~ Genotype,
    ref.group = "Dp16", # Dp16 as group1 so that "greater" tests Dp16 > WT
    paired = FALSE,
    alternative = "greater",
    p.adjust.method = "none", # BH correction applied manually below
    detailed = TRUE
  )
#

## 4.3 Rank-based effect sizes ----
# r = |Z| / sqrt(N) (rstatix returns the absolute value, so r has no direction)
effsize_component <- component_test_data |>
  rstatix::wilcox_effsize(
    formula = Score_value ~ Genotype,
    ref.group = "Dp16",
    paired = FALSE,
    detailed = TRUE
  )
#

## 4.4 Combine tests and effect sizes ----
# All components are kept; p, effsize and signed_r are NA for components removed in section 4.1
component_results <- component_data |>
  distinct(Age, Tissue, Score_name) |>
  left_join(wilcox_component, by = join_by(Age, Tissue, Score_name)) |>
  left_join(
    effsize_component |> select(Age, Tissue, Score_name, .y., group1, group2, effsize, magnitude),
    by = join_by(Age, Tissue, Score_name, .y., group1, group2)
  ) |>
  mutate(signed_r = effsize * sign(statistic - n1 * n2 / 2)) |>
  mutate(BH_padj = p.adjust(p, method = "BH"), .by = c(Age, Tissue)) |> # BH across components within each tissue-timepoint (NA p excluded)
  arrange(Age, Tissue, Score_name)
component_results
#


# 5. Component score plots and effect size heatmaps ----
#
## 5.1 Individual score plots for selected significant components ----
component_plot_data <- component_data |>
  inner_join(component_results |> select(Age, Tissue, Score_name, p), by = join_by(Age, Tissue, Score_name)) |>
  filter(p < 0.05) |>
  mutate(Score_label = paste0(Score_name, "*"))
#
set.seed(1234) # jittered point positions
component_1mo_plot <- component_plot_data |>
  filter(Age == "1mo", Tissue %in% c("Bone marrow", "Lung")) |> # CUSTOMIZE tissues to plot
  individual_score_plot(
    facets = Tissue ~ Score_label,
    title = "Component score by genotype across sig tissues - 1mo (Bone marrow and Lung)",
    y_label = "Component score",
    ncol = 3
  )
component_1mo_plot
ggsave(filename = here("plots", paste0(out_file_prefix, "component_1mo_genotype.pdf")),
       plot = component_1mo_plot, device = cairo_pdf, width = 10, height = 3, units = "in")
#
set.seed(1234) # jittered point positions
component_4mo_plot <- component_plot_data |>
  filter(Age == "4mo", Tissue %in% c("Lung", "Bone marrow", "Kidney")) |> # CUSTOMIZE tissues to plot
  individual_score_plot(
    facets = Tissue ~ Score_label,
    title = "Component score by genotype across sig tissues - 4mo (Bone marrow, Lung and Kidney)",
    y_label = "Component score",
    ncol = 6
  )
component_4mo_plot
ggsave(filename = here("plots", paste0(out_file_prefix, "component_4mo_genotype.pdf")),
       plot = component_4mo_plot, device = cairo_pdf, width = 12, height = 3, units = "in")
#

## 5.2 Effect size heatmaps ----
component_heatmap_data <- component_results |>
  filter(!is.na(effsize), !is.na(p)) |>
  mutate(
    Age_group = factor(Age, levels = c("1mo", "4mo"), labels = c("1 month", "4 month")),
    Tissue = factor(Tissue,
                    levels = c("Lung", "Bone marrow", "Kidney",
                               "Skeletal muscle", "Small intestines",
                               "Adipose", "Colon", "Bone", "Heart",
                               "Spleen", "Skin", "Brain")),
    Tissue_Score = paste(Tissue, Score_name, sep = "__")
  ) |>
  arrange(Tissue, Score_name) |>
  mutate(Tissue_Score = fct_inorder(Tissue_Score))
#
### Effect size (r, magnitude only) ----
# Color scale set to the data: 0 to the maximum r (rounded up to 2 decimals)
component_max_r <- component_heatmap_data |> pull(effsize) |> max() |> (\(x) ceiling(x * 100) / 100)()
component_heatmap <- component_heatmap_data |>
  as_tibble() |> # need to strip rstatix_test class to prevent tidyHeatmap error
  group_by(Tissue) |>
  tidyHeatmap::heatmap(
    .row = Tissue_Score,
    .column = Age_group,
    .value = effsize,
    scale = "none",
    palette_value = circlize::colorRamp2(
      seq(0, component_max_r, length.out = 9),
      RColorBrewer::brewer.pal(9, "Reds")
    ),
    heatmap_legend_param = list(
      title = "effect size (r)",
      color_bar = "continuous",
      at = seq(0, component_max_r, length.out = 5)
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    border = TRUE,
    row_labels = sub("^.*__", "", levels(component_heatmap_data$Tissue_Score)),
    column_names_rot = 45
  ) |>
  layer_asterisk(p < 0.05) |>
  tidyHeatmap::wrap_heatmap()
component_heatmap
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_wilcox_results_component.pdf")),
       plot = component_heatmap, device = cairo_pdf, width = 4.25, height = 7, units = "in")
#
### Signed effect size (signed r) ----
# Same row order as the magnitude heatmap; red = higher scores in Dp16, blue = higher scores in WT
component_max_signed_r <- component_heatmap_data |> pull(signed_r) |> abs() |> max() |> (\(x) ceiling(x * 100) / 100)()
component_heatmap_signed <- component_heatmap_data |>
  as_tibble() |> # need to strip rstatix_test class to prevent tidyHeatmap error
  group_by(Tissue) |>
  tidyHeatmap::heatmap(
    .row = Tissue_Score,
    .column = Age_group,
    .value = signed_r,
    scale = "none",
    palette_value = circlize::colorRamp2(
      seq(-component_max_signed_r, component_max_signed_r, length.out = 11),
      RColorBrewer::brewer.pal(11, "RdBu") |> rev()
    ),
    heatmap_legend_param = list(
      title = "signed r",
      color_bar = "continuous",
      at = seq(-component_max_signed_r, component_max_signed_r, length.out = 5)
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    border = TRUE,
    row_labels = sub("^.*__", "", levels(component_heatmap_data$Tissue_Score)),
    column_names_rot = 45
  ) |>
  layer_asterisk(p < 0.05) |>
  tidyHeatmap::wrap_heatmap()
component_heatmap_signed
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_wilcox_results_component_signed.pdf")),
       plot = component_heatmap_signed, device = cairo_pdf, width = 4.25, height = 7, units = "in")
#


# 6. Selected component plots ----
#
## 6.1 BALT hyperplasia ----
set.seed(1234) # jittered point positions
balt_hyperplasia_plot <- component_data |>
  filter(Score_name == "BALT_hyperplasia") |>
  left_join(component_results |> select(Age, Tissue, Score_name, p), by = join_by(Age, Tissue, Score_name)) |>
  mutate(Age_label = if_else(!is.na(p) & p < 0.05, paste0(Age, "*"), Age)) |>
  individual_score_plot(
    facets = ~ Age_label,
    title = "Lung - BALT hyperplasia",
    y_label = "Component score",
    ncol = 3
  )
balt_hyperplasia_plot
# Optional export
# ggsave(filename = here("plots", paste0(out_file_prefix, "selected_balt_hyperplasia.pdf")),
#        plot = balt_hyperplasia_plot, device = cairo_pdf, width = 8.5, height = 3, units = "in")
#


# 7. Export results ----
# Wilcoxon tests with effect sizes (effsize = r; signed_r > 0 = higher in Dp16)
write_tsv(composite_results, here("results", paste0(out_file_prefix, "composite_wilcox.txt")))
write_tsv(component_results, here("results", paste0(out_file_prefix, "component_wilcox.txt")))
list(
  "Composite" = composite_results,
  "Component" = component_results
) |>
  export_excel(filename = "wilcox_results")
#


# Save workspace ----
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace
# To reload:
# load(here("rdata", paste0(out_file_prefix, ".RData")))
#


# Session info ----
writeLines(capture.output(sessionInfo()), con = here("results", paste0(out_file_prefix, "session_info.txt")))
#
