suppressPackageStartupMessages({
  require(pheatmap);
  require(RColorBrewer);
  require(viridis)
})
source('src/R/util.load.cts.R')

plot.pheatmap = function(cts_file, sample.meta.file, variables, gene_fn, outFig, logTrans=T) {
  # this file plots heatmap for interested gene expression in all of the samples
  # cts_file: the normalized count file, first column is gene name, column names are samples
  # sample.meta.file: the meta data file, first column needs to be sample names corresponds to cts_file
  # variables: columns in sample.meta.file that want to merge as groups to plot
  # gene_fn: genes that want to show in the heatmap, each line needs to be one gene.
  # outFig: output figure
  # logTrans: if do log2 transform of the count data
  cts.mtx = read.combined_cts(cts_file)
  sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')
  genes = read.table(gene_fn)$V1
  genes_df = cts.mtx[rownames(cts.mtx) %in% genes,]
  col_groups = sample.meta$mergeCond
  mat_col <- data.frame(group = col_groups)
  rownames(mat_col) <- colnames(genes_df)
  
  genes_df = genes_df[,rownames(mat_col)[order(mat_col$group)]] # sort by column
  mat_col = data.frame(group=mat_col[rownames(mat_col)[order(mat_col$group)],])
  rownames(mat_col) = colnames(genes_df)
  
  mat_colors <- list(group = brewer.pal(length(unique(col_groups)), "Set1"))
  names(mat_colors$group) <- unique(col_groups)
  
  if (logTrans) {genes_df = log2(genes_df + 0.01)}
  png(outFig)
  pheatmap(
    mat               = genes_df,
    color             = inferno(9),
    border_color      = NA,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "",
    cluster_cols=FALSE
  )
  dev.off()
}




