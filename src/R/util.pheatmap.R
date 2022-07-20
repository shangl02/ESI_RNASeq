suppressPackageStartupMessages({
  require(pheatmap);
  require(RColorBrewer);
  require(viridis)
})
source('src/R/util.load.cts.R')

plot.pheatmap = function(norm_cts_fn, sample_meta_fn, variables, species, gene_fn, outFig, logTrans=T,show_rownames=T, z_score='none') {
  ## notice: In meta data file, a column nameed 
  # this file plots heatmap for interested gene expression in all of the samples
  # cts_file: the normalized count file, first column is gene name, column names are samples
  # sample.meta.file: the meta data file, first column needs to be sample names corresponds to cts_file
  # variables: columns in sample.meta.file that want to merge as groups to plot
  # gene_fn: genes that want to show in the heatmap, each line needs to be one gene.
  # outFig: output figure
  # logTrans: if do log2 transform of the count data
  cts.mtx = read.table(norm_cts_fn, sep='\t', header=T,row.names=1)
  sample.meta <- read.csv(sample_meta_fn, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')
  # annotate ensembl id to symbol
  anno <- addAnnoByVector(truncateEnsemblID(rownames(cts.mtx)), species)
  symbol = toupper(make.unique(anno$symbol))
  keep = complete.cases(symbol)
  cts.mtx = cts.mtx[keep,]
  rownames(cts.mtx) = symbol[keep]
  genes = toupper(read.table(gene_fn)$V1)
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
  mat_colors <- list( group= colorRampPalette(brewer.pal(8, "Dark2"))(num_colors))
  names(mat_colors$group) <- unique(col_groups$mergeCond[order(col_groups$plot_order)])
  
  if (logTrans) {genes_df = log2(genes_df + 0.01)}
  if (z_score) {
    z_scale = 'row'
  } else {
    z_scale='none'
  } 
  
  png(outFig)
  pheatmap(
    mat               = genes_df,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color      = NA,
    show_colnames     = FALSE,
    show_rownames     = show_rownames,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "",
    cluster_cols=F,
    cluster_rows=T,
    scale = z_scale
  )
  dev.off()
}



plot_pheatmap_group_mean = function(norm_cts_fn, sample_meta_fn, variables, species, gene_fn, deseq2_fn, outFig, logTrans=T,show_rownames=T, z_score='none') {
  # this function plot heatmap for averaged expression for each group
  # genes are ranked based on the log2foldchange of the input deseq2 result file
  
  # get gene expression data
  cts.mtx = read.table(norm_cts_fn, sep='\t', header=T,row.names=1)
  sample.meta <- read.csv(sample_meta_fn, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')
  # annotate ensembl id to symbol
  anno <- addAnnoByVector(truncateEnsemblID(rownames(cts.mtx)), species)
  symbol = toupper(make.unique(anno$symbol))
  keep = complete.cases(symbol)
  cts.mtx = cts.mtx[keep,]
  rownames(cts.mtx) = symbol[keep]
  genes = toupper(read.table(gene_fn)$V1)
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
  if (logTrans) {genes_df = log2(genes_df + 0.01)}
  if (z_score) {
    z_scale = 'row'
  } else {
    z_scale='none'
  }
  
  # read in DESeq2 results and sort by log2foldchange
  df = read.csv(fn)
  df = df[order(df$log2FoldChange, decreasing=T),]
  df$symbol = toupper(df$symbol)
  df = df[df$symbol %in% genes, ]
  genes_lfc = df$symbol
  genes_df = genes_df[genes_lfc,]
  by_group <- sample.meta %>% group_by(plot_order,Condition) %>%
    summarize(samples=list(Sample))
  
  # get average expression for each group
  ave_df = matrix(nrow = nrow(genes_df), ncol = length(by_group))
  colnames(ave_df) = by_group$Condition
  rownames(ave_df) = rownames(genes_df)
  for (i in 1:nrow(by_group)){
    name = by_group$Condition[[i]]
    sps = by_group$samples[[i]]
    ave_df[,name] = rowMeans(genes_df[,sps])
  }
  
  # plot
  mat_col = data.frame(group=colnames(ave_df))
  rownames(mat_col) <- mat_col$group
  mat_colors <- list(group= colorRampPalette(brewer.pal(8, "Dark2"))(ncol(ave_df)))
  names(mat_colors$group) <- colnames(ave_df)
  
  png(outFig)
  pheatmap(
    mat               = ave_df,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color      = NA,
    show_colnames     = FALSE,
    show_rownames     = show_rownames,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "",
    cluster_cols=F,
    cluster_rows=F,
    scale = z_score
  )
  dev.off()
}

