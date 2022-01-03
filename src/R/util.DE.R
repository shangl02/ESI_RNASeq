suppressPackageStartupMessages({
  require(dplyr);
  require(ggplot2);
  require(DESeq2);
  require(limma);
  require(edgeR);
})

DE.DESeq = function(dds) {
  dds <- DESeq(dds)
  res <- results(dds)
  return(list(dds, res))
}

DE.Voom = function(dge, model_matrix) {
  v <- voom(counts=dge, design=model_matrix) 
  fit <- lmFit(v, model_matrix)
  cont.matrix <- makeContrasts(GroupGroup1-GroupControl, levels=model_matrix)
  fit1 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit1)
  p <- plotSA(fit2, main="Final model: Mean-variance trend")
  return(list(fit1, fit2, p))
}