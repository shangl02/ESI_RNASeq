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
  cts.mtx = read.table(cts_file, sep='\t', header=T,row.names=1)
  sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')
  # annotate ensembl id to symbol
  anno <- addAnnoByVector(truncateEnsemblID(rownames(cts.mtx)), species)
  symbol = make.unique(anno$symbol)
  keep = complete.cases(symbol)
  cts.mtx = cts.mtx[keep,]
  rownames(cts.mtx) = symbol[keep]
  genes = read.table(gene_fn)$V1
  overlap_genes = c()
  for (g in genes) {
    if (g %in% rownames(cts.mtx)) {
      overlap_genes = c(overlap_genes, g)
    }
  }
  genes_df = cts.mtx[c(overlap_genes),]
  col_groups = sample.meta[, c('mergeCond','plot_order')]
  mat_col = sample.meta[,c('Sample','mergeCond','plot_order')]
  rownames(mat_col) <- mat_col$Sample
  
  genes_df = genes_df[,rownames(mat_col)[order(mat_col$plot_order)]] # sort by column
  mat_col = mat_col[order(mat_col$plot_order),]
  mat_col = data.frame(group=mat_col$mergeCond,row.names=rownames(mat_col))
  
  num_colors = length(unique(col_groups$mergeCond))
  mat_colors <- list( group= colorRampPalette(brewer.pal(9, "Dark2"))(num_colors))
  names(mat_colors$group) <- unique(col_groups$mergeCond[order(col_groups$plot_order)])
  
  if (logTrans) {genes_df = log2(genes_df + 0.01)}
  png(outFig)
  pheatmap(
    mat               = genes_df,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color      = NA,
    show_colnames     = FALSE,
    show_rownames     = FALSE,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "",
    cluster_cols=FALSE,
    cluster_rows=T
  )
  dev.off()
}




