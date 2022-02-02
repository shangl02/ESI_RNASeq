source('src/R/util.DE.R')
suppressPackageStartupMessages({
  require(ggplot2);
})


expression_boxplot = function(tpm_fn, gene_fn, sample_meta_fn, species, figure){
  genes = read.table(gene_fn)$V1
  sample.meta <- read.csv(sample_meta_fn, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')
  samples = sample.meta$Sample
  
  norm_df = read.csv(norm_ctx_fn,row.names = 1)
  norm_df = norm_df[,samples]
  norm_df = addAnno(norm_df, species)[[1]]
  gene_df = as.data.frame(t(norm_df[norm_df$symbol %in% genes,]))
  colnames(gene_df) = gene_df['symbol',]
  gene_df$Sample = rownames(gene_df)
  gene_df = merge(gene_df, sample.meta[,c('Sample','mergeCond')], by='Sample')
  
  count_df = data.frame(matrix(ncol=3, nrow=0))
  for (g in genes) {
    tryCatch(
      {
        cts = cbind(g,gene_df$mergeCond, gene_df[,g])
        count_df = rbind(count_df, cts)
      },
      error = function(err) {
        print(paste("Error: ", err))
      }
    )
  }
  colnames(count_df) = c('gene', 'condition', 'tpm')
  count_df$tpm = as.numeric(count_df$tpm)
  
  png(figure, width = 600, height = 480)
  ggplot(count_df, aes(x=gene, y=tpm,fill=condition)) + geom_boxplot() + 
    facet_wrap(~gene, scale="free")
  dev.off()
}
