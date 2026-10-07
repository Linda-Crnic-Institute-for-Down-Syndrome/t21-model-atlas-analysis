# Common helper functions for DESeq2 differential expression analysis

## Setting and modifying default theme for plots
theme_set(theme_gray(base_size = 12, base_family = "Arial") +
            theme(
              panel.border = element_rect(colour="black", fill = "transparent"), 
              plot.title = element_text(face="bold", hjust = 0),
              axis.text = element_text(color="black", size = 14), 
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              panel.background = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank(),
              strip.background = element_blank(), # facet label borders
              legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
            ))
#

## Density color function
getDenCols <- function(x, y, transform = TRUE) { # set to TRUE if using log2 transformation of data
  if(transform) {
    df <- data.frame(log2(x), log2(y))
  } else{
    df <- data.frame(x, y)
  }
  z <- grDevices::densCols(df, colramp = grDevices::colorRampPalette(c("black", "white")))
  df$dens <- grDevices::col2rgb(z)[1,] + 1L
  cols <-  grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  return(df$dens)
} # End of function

## Excel export function
export_excel <- function(named_list, filename = "") {
  wb <- openxlsx::createWorkbook()
  ## Loop through the list of split tables as well as their names
  ## and add each one as a sheet to the workbook
  Map(function(data, name){
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, data)
  }, named_list, names(named_list))
  ## Save workbook to working directory
  openxlsx::saveWorkbook(wb, file = here("results", paste0(out_file_prefix, filename, ".xlsx")), overwrite = TRUE)
  cat("Saved as:", here("results", paste0(out_file_prefix, filename, ".xlsx")))
} # end of function

# get size factors function ----
get_size_fcts <- function(x) {
  name <- deparse(substitute(x))
  x |>  
    sizeFactors() |>  
    enframe(name = "Sampleid", value = "SizeFactor")
}

# get standard ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# AIC functions ----
get_aic <- function(x) {
  n_params <- ncol(model.matrix(design(x), colData(x)))
  mcols(x) %>%
    as_tibble(rownames = "Geneid") %>%
    select(Geneid, deviance) %>%
    mutate(AIC = deviance + 2 * n_params)
}
#
plot_aic_diff <- function(aic1, aic2, subtitle = "") {
  aic1 %>% select(Geneid, AIC1 = AIC) %>%
    inner_join(aic2 %>% select(Geneid, AIC2 = AIC)) %>%
    mutate(AIC_diff = AIC1 - AIC2) %>%
    arrange(-AIC_diff) %>%
    mutate(Geneid = fct_inorder(Geneid)) %>%
    ggplot(aes(Geneid, AIC_diff)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-10, 10), linetype = 2) +
    geom_point() +
    theme(aspect.ratio = 0.8, axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = "AIC comparison", subtitle = subtitle, x = "Genes")
}
#

# Dendrogram and distance heatmap function -----
plotDendClust2 <- function(matrix, adjustment = "", color_var = "", clust_method = "complete", print_dend = TRUE) {
  if (dimnames(matrix)[[2]] %>% length() < 4) stop("Not enough samples - skipping dendogram and dist-heatmap")
  # Calc distance matrix for samples (not genes):
  dist <- matrix |> 
    t() |> 
    dist(method = "euclidean", diag = TRUE, upper = TRUE)
  dist_tidy <- dist |> 
    broom::tidy() |> 
    rename(Sampleid = item1, Sampleid2 = item2)
  # Hierarchical clustering and dendrogram
  dend <- dist |> 
    hclust(method = clust_method) |> # "complete", "ward.D" etc
    as.dendrogram()
  # Color labels of dendrogram object
  dendextend::labels_colors(dend) <- colnames(matrix) %>% 
    enframe(name = NULL, value = "Sampleid") %>% 
    inner_join(meta_data %>% select(Sampleid, color_var)) %>%
    inner_join(
      standard_colors %>% enframe(name = color_var, value = "color"), # REQUIRES CUSTOMIZATION
    ) %>% 
    pull(color) %>% 
    .[order.dendrogram(dend)]
  # Print dendrogram
  if (print_dend) dend %>% plot(main = paste("Sample Dendrogram -", adjustment))
  # Generate heatmap
  hm <- dist_tidy |> 
    inner_join(meta_data) |> 
    tidyHeatmap::heatmap(
      Sampleid,
      Sampleid2, 
      distance,
      palette_value = RColorBrewer::brewer.pal(3, "Blues") |> rev(),
      # palette_value = viridis::viridis(3),
      cluster_rows = dend,
      cluster_columns = dend,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_title = NULL,
      column_title = NULL
    ) |> 
    add_tile(
      !!color_var,
      palette = standard_colors
    ) |> 
    wrap_heatmap() +
    labs(
      title = paste("Sample-Sample Distances -", adjustment)
    )
  print(hm)
} # End of function

# PCA plot function:
plotPCA_custom2 <- function(x, PCA_by="variance", save_PCs = FALSE, plot_title="", subtitle="", color_var, labels=FALSE, shapes=NULL, x_lower_lim=-20, x_upper_lim=20, y_lower_lim=-20, y_upper_lim=20, top_n_variable="all") {
  components_to_plot <- c("PC1","PC2")
  # This section pulls the rows to be used for PCA (eg top 500 most-variable) from x (x should be variance-stabilized counts ± multivariable correction)
  if (PCA_by == "variance"){
    Pvars <- rowVars(x)
    if(top_n_variable=="all"){select <- order(Pvars, decreasing=TRUE)} else {
      select <- order(Pvars, decreasing=TRUE)[seq_len(min(top_n_variable, length(Pvars)))]
    }
    data_for_PCA <- x[select, ]
  }else if(PCA_by == "genelist"){
    genelist <- as.character(read.delim(genelist_file, header=FALSE)$V1)
    data_for_PCA <- subset(x, row.names(x) %in% genelist)
  }else if(PCA_by == "all"){
    data_for_PCA <- x
  }
  # calculate principal components
  prin_comp <- prcomp(t(data_for_PCA), scale=F) # set to FALSE as VST/rlog are already transformed
  # calculate percentage variance explained for each Component
  percentVar <- round(100*prin_comp$sdev^2/sum(prin_comp$sdev^2),2)
  prin_comp_percentage_explained <- data.frame(percentVar, row.names=colnames(prin_comp$x))
  # generate scree plot for fist 10 PCs
  top10_sum <- prin_comp_percentage_explained %>% 
    as_tibble(rownames = "Principal_component") %>% 
    slice_max(percentVar, n = 10) %>% 
    summarize(sum(percentVar)) %>% 
    pull() %>% 
    round()
  scree <- prin_comp_percentage_explained %>% 
    as_tibble(rownames = "Principal_component") %>% 
    slice_max(percentVar, n = 10) %>% 
    arrange(-percentVar) %>% 
    mutate(Principal_component = fct_inorder(Principal_component)) %>% 
    ggplot(aes(Principal_component, percentVar)) +
    geom_col() +
    labs(
      x = NULL, y = "% variance",
      title = "Top 10 principal components",
      subtitle = paste0(subtitle, "; variance acccounted for: ", top10_sum, "%")
    )
  # make a tbl with the principal components and append design variables
  prin_comp_x <- prin_comp$x %>% 
    as_tibble(rownames="Sampleid") %>%
    inner_join(meta_data)
  if(save_PCs) {
    assign(paste0(subtitle, "_PC_loadings"), prin_comp_x, pos = 1)
    assign(paste0(subtitle, "_PC_percVar"), prin_comp_percentage_explained %>% as_tibble(rownames = "PC"), pos = 1)
  }
  # SHAPES - needs to be a factor
  if (!is_empty(shapes)) {
    prin_comp_x <- prin_comp_x %>% 
      mutate(
        !!shapes := as_factor(.data[[shapes]])
      )
    if((prin_comp_x %>% pull(shapes) %>% levels() %>% length()) > 5) warning("Too many levels for shapes")
  }
  #
  xlabel <- paste0(components_to_plot[1], " (", prin_comp_percentage_explained[components_to_plot[1],],"%)")
  ylabel <- paste0(components_to_plot[2], " (", prin_comp_percentage_explained[components_to_plot[2],],"%)")
  # generate plot
  pca_plot <- ggplot(
    data = prin_comp_x,
    aes_string(x=components_to_plot[1], y=components_to_plot[2], color = color_var, shape = shapes)) +
    # scale_color_manual(values=standard_colors) +
    # scale_color_brewer(palette = "Paired") +
    geom_hline(yintercept=0, colour="black", linetype=1, size=0.3, alpha =1) + 
    geom_vline(xintercept=0, colour="black", linetype=1, size=0.3, alpha =1) +
    geom_point(size=4, alpha=0.8) +
    # scale_y_continuous(limits=c(y_lower_lim, y_upper_lim), breaks=seq(-120,120,by=10), expand=c(0, 0)) +
    # scale_x_continuous(limits=c(x_lower_lim, x_upper_lim), breaks=seq(-120,120,by=10), expand=c(0, 0)) +
    scale_y_continuous(limits=c(y_lower_lim, y_upper_lim)) +
    scale_x_continuous(limits=c(x_lower_lim, x_upper_lim)) +
    # scale_y_continuous() +
    # scale_x_continuous() +
    labs(title=plot_title,
         subtitle=subtitle,
         y=ylabel, 
         x=xlabel) +
    theme(aspect.ratio = 1)
  # add sample labels (off by default)
  if (labels) {
    pca_plot <- pca_plot + 
      geom_text_repel(aes(label=Sampleid),
                      color="black",
                      force=10,
                      size=3,
                      min.segment.length = 0.1
      )
  }
  assign(paste0(subtitle, "_pca_plot"), pca_plot, pos = 1)
  # Draw the plots
  scree %>% print()
  pca_plot %>% print() # print required to work with for loops 
} # End of function

# Function to run lm for variables against PCA loadings ----
pca_lm_function <- function(loadings_data, percVar, f) {
  loadings_data %>% 
    select(Sampleid, matches("PC")) %>% # or subselect PCs if needed
    inner_join(groups) %>% 
    inner_join(colData(dds) %>% as_tibble()) %>% 
    pivot_longer(cols = matches("PC"), names_to = "PC", values_to = "loadings") %>% 
    nest(data = -PC) %>% 
    inner_join(percVar) %>% 
    mutate(
      formula = f, # not required, but stores model formula
      fit = map(data, ~ lm(f, data = .)),
      tidied = map(fit, broom::tidy)
    ) %>% 
    unnest(tidied) %>% 
    filter(!str_detect(term, "Intercept")) %>% 
    arrange(p.value)
}


# Functions to print results summaries ± independent filtering -------
# CONTRASTS VERSION
get_results_sum <- function(x, contrast = c(predictor, "denom", "num"), show_ind_filt_off=TRUE) {
  message("------------\nResults summary for ", contrast[2], " vs. ", contrast[3])
  message("Model formula: ", design(x), "\n------")
  if (show_ind_filt_off) {
    cat("independentFiltering OFF, Cooks cutoff ON")
    x %>%
      results(
        contrast,
        # name = name,
        cooksCutoff=TRUE,
        independentFiltering=FALSE
      ) %>%
      summary()
    message("----\n")
  } else {
    cat("independentFiltering ON, Cooks cutoff ON")
    x %>%
      results(
        contrast,
        # name = name,
        cooksCutoff=TRUE,
        independentFiltering=TRUE
      ) %>%
      summary()
    message("----\n")
  }
}
#

# Functions to generate final results tbls in our format (use with for loops to get all comparisons) ------
# CONTRASTS VERSION
get_results_tbl <- function(x, contrast = "", cooks = TRUE, ind_filt = TRUE, shrink_type = "apeglm") {
  #
  message("------\nGenerating DEseq2 results for ", paste0(contrast[1], ": ", contrast[2], " vs. ", contrast[3]))
  message("Model formula: ", design(x), "\n------")
  ## relevel variable of interest (and re-run nbinomWaldTest()) to ensure that
  # coefficient of interest is available (check with resultsNames()) for
  # lfcShrink() - otherwise need to use contrast argument which does not allow
  # use of apeglm shrinkage
  message(paste0("1. Re-leveling '", contrast[1], "' with '", contrast[3], "' as reference level"))
  # no longer hardcoded - see ?`SummarizedExperiment-class` for info on accessors for colData
  x[[contrast[1],]] <- x[[contrast[1],]] |> relevel(contrast[3])
  # x[[contrast[1],]] <- x[[contrast[1],]] |> fct_relevel(c(contrast[3], contrast[2])) # Safer with more levels?
  #
  message("2. Running negative binomial Wald test")
  x <- x %>% nbinomWaldTest(quiet=TRUE)
  ## Get results without LFC shrinkage (default)
  res <- x %>%
    results(
      contrast,
      cooksCutoff = cooks,
      independentFiltering = ind_filt
    )
  res_tbl <- res %>% # Could also use biobroom::tidy() on results or dds but has less useful colnames
    as.data.frame() %>%
    as_tibble(rownames = "Geneid")
  ## Apply LFC shrinkage
  message("3. Calculating log2 fold-change shrinkage")
  # NOTES ON LOG FOLD CHANGE SHRINKAGE
  # earlier DESeq2 versions (<1.16) carried out LFC shrinkage by
  # default(betaPrior=TRUE); more recent versions set betaPrior=FALSE and use
  # lfcShrink() in a separate step.
  # https://support.bioconductor.org/p/95695/
  # However: 
  # with betaPrior=TRUE, p-values are calculated for the shrunken LFC, while
  # betaPrior=FALSE + subsequent lfcShrink() calculates p-values based on
  # un-shrunken LFC and only shrinks them afterwards.
  # Also:
  # Difference between using lfcShrink() with the coef argument and using lfcShrink() with the contrast argument
  # lfcShrink(dds=dds, coef=2, res=res) OR lfcShrink(dds=dds, contrast=c("condition","B","A"), res=res)
  # https://support.bioconductor.org/p/98833/#98837 From Michael Love:
  # They are not identical. Using contrast is similar to what DESeq2 used to do:
  # it forms an expanded model matrix, treating all factor levels equally, and
  # averages over all distances between all pairs of factor levels to estimate the
  # prior. Using coef, it just looks at that column of the model matrix (so
  # usually that would be one level against the reference level) and estimates the
  # prior for that coefficient from the distribution of those MLE of coefficients.
  # I implemented both for lfcShrink, because 'contrast' provides backward support
  # (letting people get the same coefficient they obtained with previous
  # versions), while future types of shrinkage estimators will use the 'coef'
  # approach, which is much simpler.
  # my current recommendation would be to use the p-values from un-shrunken LFC
  # and then use the shrunken LFC for visualization or ranking of genes. This is
  # the table you get with default DESeq => results => lfcShrink. If you want to
  # be future-proof, I'd go with 'coef' with lfcShrink. All the methods I'm
  # planning on adding (ours and others) really just want to shrink one
  # coefficient at a time, not do the expanded model matrix thing.
  res_shrink <- x %>%
    lfcShrink(
      coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep="_"), # REPLACE dashes or this fails
      type = shrink_type, # "apeglm" or "normal" or "ashr"; apeglim and ashr are better at preserving large LFCs
      parallel = TRUE,
      res = res
    )
  res_shrink_tbl <- res_shrink %>%
    as.data.frame() %>%
    as_tibble(rownames="Geneid")
  # Combine and mutate to get final results tbl
  message("4. Assembling final results table\n------")
  res_tbl %>%
    inner_join(res_shrink_tbl %>% select(Geneid, log2FoldChange), by = "Geneid") %>%
    dplyr::rename(
      log2FoldChange = log2FoldChange.x,
      log2FoldChange_adj = log2FoldChange.y
    ) %>%
    mutate(
      FoldChange = 2^log2FoldChange,
      FoldChange_adj = 2^log2FoldChange_adj,
      Model = design(x) |> paste(collapse = "")
    ) %>%
    dplyr::select(
      Geneid:baseMean,
      Model,
      FoldChange,
      log2FoldChange,
      FoldChange_adj,
      log2FoldChange_adj,
      pvalue,
      padj
    ) %>%
    arrange(padj) %>%
    inner_join(gene_anno, by="Geneid") %>%
    select(Gene_name = gene_name, chr, everything())
}


# MA plot functions -----
plotDEgg <- function(res, sig, title, subtitle){
  num_sig_up=length(which(!is.na(res$padj[res$padj < sig & res$log2FoldChange > 0])))
  num_sig_down=length(which(!is.na(res$padj[res$padj < sig & res$log2FoldChange < 0])))
  # get max for y-axis
  y_lim <- res %>% summarize(max = max(log2FoldChange), min = min(log2FoldChange)) %>% abs() %>% max() %>% ceiling()
  # generate plot
  res %>% 
    mutate(Density = getDenCols(baseMean, FoldChange, transform = TRUE)) %>% # set to TRUE if using log2 transformation of data
    arrange(Density) |> 
    ggplot(aes(log2(baseMean), log2(FoldChange))) +
    geom_hline(yintercept=0, color="black", linetype=2) +
    # geom_hex(bins=120) +
    # scale_fill_gradient(low="black", high="white") +
    geom_point(aes(color = Density), size = 1) +
    scale_color_viridis_c() +
    geom_point(data=(res %>% filter(padj < 0.1)), aes(log2(baseMean), log2(FoldChange)), color="red", size = 1) +
    # ylim(c(-6, 6)) +
    ylim(c(-y_lim, y_lim)) + # set symmetrical y-axis
    labs(title=title,
         subtitle=paste0(subtitle, "; [Up: ", num_sig_up, ", Down: ", num_sig_down, "]"),
         fill="density"
    ) + 
    theme(aspect.ratio=1)
} # plots unadjusted FoldChange

plotDEgg_adj <- function(res, sig, title, subtitle){
  num_sig_up=length(which(!is.na(res$padj[res$padj < sig & res$log2FoldChange > 0])))
  num_sig_down=length(which(!is.na(res$padj[res$padj < sig & res$log2FoldChange < 0])))
  # get max for y-axis
  y_lim <- res %>% summarize(max = max(log2FoldChange_adj), min = min(log2FoldChange_adj)) %>% abs() %>% max() %>% ceiling()
  # generate plot
  res %>% 
    mutate(Density = getDenCols(baseMean, FoldChange_adj, transform = TRUE)) %>% # set to TRUE if using log2 transformation of data
    ggplot(aes(log2(baseMean), log2(FoldChange_adj))) +
    geom_hline(yintercept=0, color="black", linetype=2) +
    # geom_hex(bins=120) +
    # scale_fill_gradient(low="black", high="white") +
    geom_point(aes(color = Density), size = 1) +
    scale_color_viridis_c() +
    geom_point(data=(res %>% filter(padj < 0.1)), aes(log2(baseMean), log2(FoldChange_adj)), color = "red", size = 1) +
    # ylim(c(-6, 6)) +
    ylim(c(-y_lim, y_lim)) + # set symmetrical y-axis
    labs(title=title,
         subtitle=paste0(subtitle, "; [Up: ", num_sig_up, ", Down: ", num_sig_down, "]"),
         fill="density"
    ) + 
    theme(aspect.ratio=1)
} # plots adjusted FoldChange


# signif highlight volcano plot function -----
plotVolcano <- function(res, sig, title, subtitle) {
  num_sig_up=length(which(!is.na(res$padj[res$padj < sig & res$log2FoldChange > 0])))
  num_sig_down=length(which(!is.na(res$padj[res$padj < sig & res$log2FoldChange < 0])))
  # get max for x-axis
  x_lim <- res %>% summarize(max = max(log2FoldChange), min = min(log2FoldChange)) %>% abs() %>% max() %>% ceiling()
  # generate plot
  res %>%
    ggplot(aes(log2(FoldChange_adj), -log10(padj))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = -log10(0.1), linetype = 2) +
    geom_point(data=(res %>% filter(padj >= 0.1)), aes(log2(FoldChange), -log10(padj)), color="black", size=0.9) +
    geom_point(data=(res %>% filter(padj < 0.1)), aes(log2(FoldChange), -log10(padj)), color="red", size=0.9) +
    xlim(c(-x_lim, x_lim)) + # set symmetrical y-axis
    # ylim(c(-0, -log10(res$padj) %>% max(na.rm=TRUE) %>% plyr::round_any(10, f=ceiling))) + # not needed - ggplot automaticall adds ~5%
    # ylim(c(-0, -log10(res$padj) %>% max(na.rm=TRUE))) +
    labs(title=title,
         subtitle=paste0(subtitle, "; [Up: ", num_sig_up, ", Down: ", num_sig_down, "]")
    ) + 
    theme(aspect.ratio=1.2)
}

# labelled Volcano plot function ----
volcano_plot_lab <- function(res, 
                             labels = TRUE,
                             n_labels = 2,
                             title = "", 
                             subtitle = "",
                             y_lim = c(0, NA),
                             raster = FALSE
){res <- res %>% 
  filter(!is.na(padj)) %>% # genes without padj (independent filtering / outliers) cannot be plotted
    mutate(
      color = if_else(padj < 0.1, "padj < 0.1", "All")
    )
  # Get max finite -log10(pval) and replace 0s if needed
  max_finite <- -log10(min(res$padj[res$padj > 0], na.rm = TRUE))
  res <- res %>% mutate(
    shape = if_else(padj == 0, "infinite", "finite"),
    padj = if_else(padj == 0, 10^-(max_finite * 1.05), padj)
  )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange_adj), na.rm = TRUE), min = min(log2(FoldChange_adj), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    ceiling()
  p <- res %>% 
    ggplot(aes(log2(FoldChange_adj), -log10(padj), color = color, shape = shape)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
    scale_shape_manual(values = c("infinite" = 2, "finite" = 16)) +
    guides(shape = "none") +
    # xlim(-x_lim, x_lim) +
    # ylim(y_lim) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
  if(labels == TRUE) {
    p <- p +
      geom_text_repel(data = res %>% slice_min(order_by = padj, n = n_labels), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE) + # , nudge_y = -max_finite / 5
      geom_text_repel(data = res %>% slice_min(order_by = FoldChange_adj, n = n_labels) %>% filter(padj < 0.1), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = -x_lim / 5, nudge_y = max_finite / 20, ylim = c(max_finite / 10, NA)) +
      geom_text_repel(data = res %>% slice_max(order_by = FoldChange_adj, n = n_labels) %>% filter(padj < 0.1), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = x_lim / 5, nudge_y =  max_finite / 20, ylim = c(max_finite / 10, NA))
  }
  #
  if(raster) {
    # to rasterize all points:
    # rasterize(p, layers='Point', dpi = 600)
    p <- rasterize(p, layers='Point', dpi = 600, dev = "ragg_png")
    # otherwize can rasterize individual layers using rasterize(geom_sina())
  }
  #
  return(p)
} # end of function


# Highlight chr21 genes Volcano plot function -----
volcano_plot_chr21 <- function(res, title = "", 
                               subtitle = "down in T21                                                  up in T21",
                               y_lim = c(0, NA),
                               raster = FALSE){
  res <- res %>% 
    filter(!is.na(padj)) %>% # genes without padj (independent filtering / outliers) cannot be plotted
    mutate(
      color = if_else(chr == "chr21", "chr21", "All")
    )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange_adj), na.rm = TRUE), min = min(log2(FoldChange_adj), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    ceiling()
  p <- res %>% 
    ggplot(aes(log2(FoldChange_adj), -log10(padj), color = color)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(data = . %>% filter(chr != "chr21")) +
    geom_point(data = . %>% filter(chr == "chr21")) +
    # scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
    scale_color_manual(values = c("chr21" = "#009b4e", "All" = "black")) + 
    # xlim(-x_lim, x_lim) +
    ylim(y_lim) +
    # geom_text_repel(data = res %>% filter(!is.na(padj)) %>% slice_max(order_by = FoldChange_adj, n = 4), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = 0.5, nudge_y = 10) +
    # geom_text_repel(data = res %>% filter(!is.na(padj)) %>% slice_min(order_by = FoldChange_adj, n = 4), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = -0.5, nudge_y = 10) +
    geom_text_repel(data = res %>% filter(!is.na(padj)) %>% slice_min(order_by = padj, n = 4), aes(label = Gene_name), min.segment.length = 0, show.legend = FALSE, nudge_x = 0.5, nudge_y = 0.1) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
  #
  if(raster) {
    # to rasterize all points:
    # rasterize(p, layers='Point', dpi = 600)
    p <- rasterize(p, layers='Point', dpi = 600, dev = "ragg_png")
    # otherwize can rasterize individual layers using rasterize(geom_sina())
  }
  #
  return(p)
} # end of function

# Highlight chr16 triplicated genes Volcano plot function -----
# Dp16-SPECIFIC: triplicated genes defined by chromosome AND coordinates of the
# Dp(16)1Yey duplication (GRCm38); change for other mouse models or human data
volcano_plot_chr16trip <- function(res, title = "", 
                                   subtitle = "",
                                   y_lim = c(0, NA),
                                   raster = FALSE){
  res <- res %>% 
    filter(!is.na(padj)) %>% # genes without padj (independent filtering / outliers) cannot be plotted
    mutate(
      color = case_when(
        chr == "chr16" & start >= 75540514 & end <= 97962622 ~ "Triplicated", # triplicated region of chr16
        padj < 0.1 ~ "padj < 0.1",
        TRUE ~ "Other"
      )
    )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange_adj), na.rm = TRUE), min = min(log2(FoldChange_adj), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    ceiling()
  p <- res %>% 
    ggplot(aes(log2(FoldChange_adj), -log10(padj), color = color)) + 
    geom_hline(yintercept = -log10(0.1), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(data = . %>% filter(color == "Other")) +
    geom_point(data = . %>% filter(color == "padj < 0.1")) +
    geom_point(data = . %>% filter(color == "Triplicated")) +
    scale_color_manual(values = c("Triplicated" = "#0000CC", "padj < 0.1" = "grey50", "Other" = "black")) +
    coord_cartesian(xlim = c(-x_lim, x_lim)) + # symmetrical x-axis; override by adding another coord_cartesian()
    ylim(y_lim) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
  #
  if(raster) {
    # to rasterize all points:
    # rasterize(p, layers='Point', dpi = 600)
    p <- rasterize(p, layers='Point', dpi = 600, dev = "ragg_png")
    # otherwize can rasterize individual layers using rasterize(geom_sina())
  }
  #
  return(p)
} # end of function

# Export DEseq results function -----
export_res <- function(results_vector) { 
  for (result in results_vector) {
    if(!str_detect(result, "multi")) {
      # generate filename
      filename <- paste0(
        out_file_prefix,
        result %>% str_replace("res_", "results_"),
        "_",
        min_cpm,
        "cpm"
      )
      path <- paste0(here("results/"), filename)
      res <- get(result)
      # write out tab-delim text file
      message("Exporting simple results to tab-delimited text file:\n", paste0(path, ".txt"))
      res %>% 
        write_tsv(file = paste0(path, ".txt"))
      # write out xlsx
      message("Exporting simple results to xlsx:\n", paste0(path, ".xlsx"))
      sheetname <- "DESeq2_results" # max length 31 char
      wb=createWorkbook()
      addWorksheet(wb, sheetname)
      textStyle=createStyle(numFmt="@")
      HS=createStyle(textDecoration="Bold")
      writeData(wb, sheetname, res, headerStyle=HS, keepNA=TRUE) # need to include NAs
      addStyle(wb, sheetname, style=textStyle, rows=1:nrow(res)+1, cols=1, gridExpand=TRUE)
      freezePane(wb=wb, sheet=sheetname, firstActiveRow=2)
      saveWorkbook(wb, file = paste0(path, ".xlsx"), overwrite=TRUE)
    } else if (str_detect(result, "multi")) {
      # generate filename
      filename <- paste0(
        out_file_prefix,
        result %>% str_replace("res_", "results_"),
        "_",
        min_cpm,
        "cpm"
      )
      path <- paste0(here("results/"), filename)
      res <- get(result) 
      # write out tab-delim text file
      message("Exporting multivariable results to tab-delimited text file:\n", paste0(path, ".txt"))
      res %>% 
        write_tsv(file = paste0(path, ".txt"))
      # write out xlsx
      message("Exporting multivariable results to xlsx:\n", paste0(path, ".xlsx"))
      sheetname <- "DESeq2_results" # max length 31 char
      wb=createWorkbook()
      addWorksheet(wb, sheetname)
      textStyle=createStyle(numFmt="@")
      HS=createStyle(textDecoration="Bold")
      writeData(wb, sheetname, res, headerStyle=HS, keepNA=TRUE) # need to include NAs
      addStyle(wb, sheetname, style=textStyle, rows=1:nrow(res)+1, cols=1, gridExpand=TRUE)
      freezePane(wb = wb, sheet=sheetname, firstActiveRow=2)
      saveWorkbook(wb, file=paste0(path, ".xlsx"), overwrite=TRUE)
    }
  }
}


# Plot dispersion estimates function -----
ggplotDispEsts <- function(x) {
  title <- paste("Dispersions for", deparse(substitute(x)))
  x_row_calcs <- mcols(x, use.names=TRUE) %>% as_tibble(rownames="Geneid")
  x_row_calcs %>% 
    ggplot(aes(baseMean, dispGeneEst)) +
    geom_point(aes(color="geneEst"), size=0.8) +
    geom_point(aes(baseMean, dispersion, color="final"), shape=1, size=0.8) +
    geom_point(aes(baseMean, dispFit, color="fitted"), size=0.8) +
    scale_x_continuous(trans="log10") +
    ylab("dispersion") +
    scale_y_continuous(trans="log10") +
    scale_colour_manual(name="", values=c(geneEst="black", final="steelblue", fitted="red")) +
    labs(title=title)
}


# Plot p-value distributons function
plotPvals <- function(res, title="") {
  title <- paste("Distribution of adjusted p-values:\n", title)
  res %>% 
    ggplot(aes(padj)) +
    geom_histogram(bins=100) +
    labs(title=title)
}

# GSEA FUNCTIONS -----
# function to get combined pos and neg GSEA results -----
run_fgsea2 <- function(geneset, ranks, weighted = FALSE) {
  library("fgsea")
  # with gseaParam = 0, results are VERY similar to original GSEA # this seems to not be operating as expected as N^0 = 1, so all ranking stats would be 1
  weight = 0
  if(weighted) weight = 1
  # Run positive enrichment
  fgseaRes_POSITIVE <- fgseaMultilevel(
    geneset, 
    ranks, 
    minSize=15, 
    maxSize=500,
    gseaParam = weight,
    # nperm = 1000,
    eps = 0.0, # fgsea has a default lower bound eps=1e-10 for estimating P-values. If you need to estimate P-value more accurately, you can set the eps argument to zero
    scoreType = "pos"
  )
  # Run negative enrichment
  fgseaRes_NEGATIVE <- fgseaMultilevel(
    geneset,
    ranks,
    minSize=15,
    maxSize=500,
    gseaParam = weight,
    # nperm = 1000,
    eps = 0.0, # fgsea has a default lower bound eps=1e-10 for estimating P-values. If you need to estimate P-value more accurately, you can set the eps argument to zero
    scoreType = "neg"
  )
  # Combine positive and negative results + re-adjust pvals
  fgseaRes_POS_NEG <- inner_join(
    fgseaRes_POSITIVE %>% 
      as_tibble(),
    fgseaRes_NEGATIVE %>% 
      as_tibble(),
    by = c("pathway"),
    suffix = c("_POS", "_NEG")
  )
  fgseaRes_COMBINED <- bind_rows(
    fgseaRes_POS_NEG %>% filter(ES_POS > abs(ES_NEG)) %>% select(pathway) %>% inner_join(fgseaRes_POSITIVE),
    fgseaRes_POS_NEG %>% filter(ES_POS < abs(ES_NEG)) %>% select(pathway) %>% inner_join(fgseaRes_NEGATIVE)
  ) %>% 
    mutate(padj = p.adjust(pval, method = "BH"))%>% 
    arrange(padj, -abs(NES))
  return(fgseaRes_COMBINED)
} # end of function

# function for moustache plot ----
moustache <- function(x, title = "", subtitle = "GSEA Hallmarks", top_n = 5, direction = "pos") {
  xlim <- x %>% select(NES) %>% range() %>% abs() %>% max() * 1.05 # set and expand xlim
  p <- x %>% 
    ggplot(aes(NES, padj)) + 
    geom_hline(yintercept = 0.1, linetype = 2) +
    geom_point() +
    theme(aspect.ratio = 1) +
    # scale_x_continuous(limits = c(-5,5)) +
    scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-20, 20, 2)) + # make x-axis symmetrical about 0
    labs(title = title,
         # subtitle = "GSEA / Hallmarks",
         subtitle = subtitle,
         y = "FDR q-value")
  if (direction == "pos" | direction == "both") {
    p <- p + geom_text_repel(data = (x %>%
                                       filter(padj < 0.1 & NES > 0) %>% 
                                       mutate(pathway = str_replace(pathway, "^HALLMARK_", ""),
                                              pathway = str_replace_all(pathway, "_", " "),
                                              pathway = str_to_title(pathway)) %>% 
                                       arrange(-NES) %>% 
                                       head(n = top_n)),
                             aes(NES, padj, label = pathway),
                             nudge_y = 0.2,
                             nudge_x = -0.5,
                             direction = "both",
                             angle = 90,
                             ylim = c(0.1, NA),
                             # hjust = 0,
                             # vjust = 1,
                             # force = 10
    )
  }
  if (direction == "neg" | direction == "both") {
    p <- p + geom_text_repel(data = (x %>%
                                       filter(padj < 0.1 & NES < 0) %>% 
                                       mutate(pathway = str_replace(pathway, "^HALLMARK_", ""),
                                              pathway = str_replace_all(pathway, "_", " "),
                                              pathway = str_to_title(pathway)) %>% 
                                       arrange(-NES) %>% 
                                       tail(n = top_n)),
                             aes(NES, padj, label = pathway),
                             nudge_y = 0.2,
                             nudge_x = -2,
                             direction = "both",
                             angle = 90,
                             ylim = c(0.1, NA),
                             # hjust = 0,
                             # vjust = 1,
                             # force = 10
    )
  }
  print(p)
  ggsave(filename = here("plots", paste0(out_file_prefix, deparse(substitute(x)), ".png")),
         height = 6,
         width = 5,
         units = "in")
  ggsave(filename = here("plots", paste0(out_file_prefix, deparse(substitute(x)), ".pdf")),
         device = cairo_pdf,
         height = 6,
         width = 5,
         units = "in")
  cat("Plot saved as: ", here("plots", paste0(out_file_prefix, deparse(substitute(x)), ".*")))
} # End of function

# Customized version of plotEnrichment -----
plotEnrichment2 <- function (pathway, stats, res, title = "") 
{
  # pathway = hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE # for testing
  # stats = ranks # for testing
  # title = "Hallmark Interferon Gamma Response"
  gseaParam = 0
  ticksSize = 0.4
  pathname <- deparse(substitute(pathway)) %>% str_remove("\\w+\\$")
  label = paste0(
    "NES = ", (res %>% filter(pathway == pathname))$NES %>% round(2),
    "\n",
    "Q = ", (res %>% filter(pathway == pathname))$padj %>% rstatix::p_format()
  )
  x_label <- length(stats)*0.99
  y_label <- ((res %>% filter(pathway == pathname))$ES)*0.95
  # Setting and modifying default theme for plots
  theme_set(theme_gray(base_size=12, base_family="Arial") +
              theme(panel.border=element_rect(colour="black", fill="transparent"), 
                    plot.title=element_text(face="bold", hjust=0),
                    axis.text=element_text(color="black", size=12), 
                    axis.text.x=element_text(angle=0, hjust=0.5),
                    # axis.text.x=element_text(angle=90, hjust=0.5),
                    # axis.text.x=element_text(angle=45, hjust=1),
                    panel.background=element_blank(),
                    panel.grid=element_blank(),
                    plot.background=element_blank()
              ) +
              # theme(strip.background=element_rect(colour="black", fill="light grey", size=1))
              theme(strip.background = element_blank()) # adjusts facet label borders
  )
  #
  rnk <- rank(-stats) # rank highest values first
  ord <- order(rnk) # get correct order
  statsAdj <- stats[ord] # ensure ranked list is ordered correctly
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam) # gets sign and multiplies by absolute value ^ gsea param
  statsAdj <- statsAdj/max(abs(statsAdj))
  zero_cross <- statsAdj[statsAdj > 0] %>% length() # New; get Zero crossing point
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(stats = statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + 
    # geom_point(color = "green", size = 0.1) + # why bother?
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0, colour = "black") + 
    geom_vline(xintercept = zero_cross, linetype = 2, color = "grey50") +
    geom_line(color = "green") + # this is the main running ES score
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + # gene set ticks
    theme(aspect.ratio = 0.7) +
    annotate(geom = "text", x = zero_cross + (n / 50), y = 0.025, label = paste0("Zero cross at ", zero_cross), hjust = 0, size=4, family="Arial") +
    annotate(geom = "text", x_label, y_label, label = label, hjust = 1, family="Arial") +
    labs(title = title, x = "Rank", y = "Enrichment score")
  g
} # End of function


# Add additional symbols for tidyHeatmap -----
layer_symbol2 <- function(.data,
                          ...,
                          symbol = "point"){
  
  .data_drame = .data@data
  
  
  symbol_dictionary = 
    list(
      point = 21,
      square = 22,
      diamond = 23,
      arrow_up = 24,
      arrow_down = 25,
      star = 8,
      asterisk = 42
    )
  
  if(!symbol %in% names(symbol_dictionary) | length(symbol) != 1) 
    stop(sprintf("tidyHeatmap says: the symbol argument must be a character string - one of %s", names(symbol_dictionary) |> paste(collapse = ", ")))
  
  # Comply with CRAN NOTES
  . = NULL
  column = NULL
  row = NULL
  
  # Make col names
  # Column names
  .horizontal = .data@arguments$.horizontal
  .vertical = .data@arguments$.vertical
  .abundance = .data@arguments$.abundance
  
  # Append which cells have to be signed
  .data@layer_symbol= 
    .data@layer_symbol |>
    bind_rows(
      .data_drame |>
        droplevels() |>
        mutate(
          column = !!.horizontal %>%  as.factor()  %>%  as.integer(),
          row = !!.vertical  %>%  as.factor() %>% as.integer()
        ) |>
        filter(...) |>
        select(column, row) |>
        mutate(shape = symbol_dictionary[[symbol]])
    )
  
  .data
  
  
}



