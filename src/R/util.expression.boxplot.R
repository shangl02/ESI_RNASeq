source('src/R/util.DE.R')
suppressPackageStartupMessages({
  require(ggplot2);
  require(forcats)
})


expression_boxplot = function(norm_ctx_fn, gene_fn, sample_meta_fn, variables, species, figure, logTrans){
  # norm_ctx_fn: normalized expression matrix
  # gene_fn: file with gene list, each line is a gene
  # sample_meta_fn: condition samples, first column name needs to be Sample
  # variables: columns in sample_meta_fn to merge the conditions
  # species: human or mouse.
  # figure: output figure name
  # logTrans: log transform of the expression matrix or not.
  genes = read.table(gene_fn)$V1
  sample.meta <- read.csv(sample_meta_fn, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')
  samples = sample.meta$Sample

  norm_df = read.table(norm_ctx_fn,row.names = 1,header=T)
  norm_df = norm_df[,samples]
  if (logTrans) {norm_df = log2(norm_df + 0.1)}
  norm_df = addAnno(norm_df, species)[[1]]
  gene_df = as.data.frame(t(norm_df[norm_df$symbol %in% genes,]))
  colnames(gene_df) = gene_df['symbol',]
  gene_df$Sample = rownames(gene_df)
  meta_col = c('Sample','mergeCond','plot_order')
  gene_df = merge(gene_df, sample.meta[,c(meta_col)], by='Sample')
  # in case there's duplicate gene names
  colnames(gene_df) = make.unique(colnames(gene_df))
  genes = setdiff(colnames(gene_df), meta_col)

  count_df = data.frame(matrix(ncol=4, nrow=0))
  for (g in genes) {
    tryCatch(
      {
        cts = cbind(g,gene_df$mergeCond, gene_df[,g],gene_df$plot_order)
        count_df = rbind(count_df, cts)
      },
      error = function(err) {
        print(paste("Error: ", err))
      }
    )
  }
  colnames(count_df) = c('gene', 'condition', 'tpm', 'order')
  count_df$tpm = as.numeric(count_df$tpm)
  count_df = count_df[order(count_df$order),]
  levels = unique(count_df$condition)
  count_df$condition = factor(count_df$condition, levels=levels)
  
  png(figure, width = 800, height = 480)
  p = ggplot(count_df, aes(x=gene, y=tpm, fill=condition)) + geom_boxplot() +
        facet_wrap(~gene, scale="free")
  if (logTrans) {p = p + labs(y = 'log2 of norm count')}
  print(p)
  dev.off()
}
