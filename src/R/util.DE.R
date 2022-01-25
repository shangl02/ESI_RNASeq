suppressPackageStartupMessages({
  require(dplyr);
  require(ggplot2);
  require(gplots);
  require(DESeq2);
  require(limma);
  require(edgeR);
  require(stringr);
  require(apeglm);
  require(gage);
  require(gageData);
  require(EnhancedVolcano)
})

DE.DESeq = function(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, species) {
  
  col <- colnames(compare_df)
  row <- nrow(compare_df)
  
  for (i in 1:row) {
    tryCatch({
      # 1st part, list the control sample and test sample
      ctrl = as.vector(compare_df[i,'control'])
      test = as.vector(compare_df[i,'test'])
      sub_cond_df = sample.meta[sample.meta$mergeCond %in% c(ctrl,test),]
      row.names(sub_cond_df) = sub_cond_df$Sample
      # sub expression df
      sub_expr_df = cts.mat[,row.names(sub_cond_df)]
      print(paste0('Process conparison ', ctrl, ' vs ', test))
      
      
      # 2nd part, do the DE for one compare
      dds <-DESeqDataSetFromMatrix(countData=round(sub_expr_df), colData=sub_cond_df,design=~mergeCond)
      keep <- rowSums(counts(dds)) >= min_total_count
      dds <- dds[keep,]
      # define control and test
      dds$mergeCond <- relevel(dds$mergeCond, ref = ctrl)
      dds <- DESeq(dds)
      res <- results(dds)
      resOrdered <- res[order(res$padj),]
      result <- resOrdered[complete.cases(resOrdered),]
      
      # create sub-folder and save results
      comparison = str_replace_all(paste(test, ctrl, sep="_VS_"), '[ :,>]','_')
      path = file.path(outdir, comparison)
      dir.create(path)
      
      l = addAnnoByVector(truncateEnsemblID(rownames(result)), species)
      result$symbol <- l$symbol
      write.csv(result, file.path(path, paste0(comparison, '.result.csv')),quote=F)
      
      # output cluster figure
      #pdf(paste(strsplit(outputFile,'\\,')[[1]][1],'.pdf',sep=""))
      
      pdf(file.path(path, paste0(comparison,'.plots.pdf')))
      rld <- rlog(dds)
      # Heat map
      library("RColorBrewer")
      library("gplots")
      #heat map of distance matrix
      distsRL <- dist(t(assay(rld)))
      mat <- as.matrix(distsRL)
      rownames(mat) <- colnames(mat) <- with(colData(dds),sub_cond_df$mergeCond)
      colours = colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
      heatmap.2(mat, trace="none", col = colours,margin=c(17,17))
      
      ## MA plot
      c <- resultsNames(dds)[2]
      res <- lfcShrink(dds, coef =  c, type="apeglm")
      plotMA(res, ylim=c(-5,5))
      
      ## topN gene cluster heatmap
      l<-addAnnoByVector(truncateEnsemblID(rownames(rld)), species)
      rownames(rld) <- l$symbol
      plot.topNVar.heatmap(rld, 20);
      
      ## Volcano plot
      p<-plot.enhancedVolcano(res, species, lfc_cutoff, alpha, paste0(test, " vs ", ctrl));
      print(p)
      #plot.volcano(result, lfc_cutoff, alpha, TRUE)
      dev.off()
      
      ## Significant Up and Down
      sig_up <- result[result$padj < alpha & result$log2FoldChange >= log2(lfc_cutoff),]
      sig_dn <- result[result$padj < alpha & result$log2FoldChange <= -log2(lfc_cutoff),]
      
      print(paste("upreg",ctrl,test,nrow(sig_up),sep=" "))
      print(paste("dnreg",ctrl,test,nrow(sig_dn),sep=" "))
      
      ## pathway analysis
      path2 <- file.path(path,'Pathway')
      dir.create(path2)
      analysis.pathway(sig_up, species, path2, 'Up')
      analysis.pathway(sig_dn, species, path2, 'Down')
    },
    error = function(err) {
      print(paste("Error: ", err))
    }, 
    finally = function(f) {
      dev.off()
    })
  }
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

plot.volcano = function(resultsObject, species, lfc_cutoff, alpha, showText=false) {
  ## Add annotation information
  addAnno(resultsObject, species)
  
  topT <- as.data.frame(resultsObject)
  
  par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
  with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col="grey", cex=0.5, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
  with(subset(topT, padj< alpha & abs(log2FoldChange)>lfc_cutoff), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.6))
  
  #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-1*lfc_cutoff, col="black", lty=4, lwd=2.0)
  abline(v=lfc_cutoff, col="black", lty=4, lwd=2.0)
  abline(h=-log10(max(topT$pvalue[topT$padj<alpha], na.rm=TRUE)), col="brown", lty=4, lwd=2.0)
  
  gn.selected <- abs(resultsObject$log2FoldChange) > lfc_cutoff & resultsObject$padj < alpha 
  if (showText){
    text(resultsObject$log2FoldChange[gn.selected],
         -log10(resultsObject$padj)[gn.selected],
         lab=rownames(resultsObject)[gn.selected ], cex=0.4)
  }
}

plot.enhancedVolcano = function(res, species, lfc_cutoff, alpha, title) {
  ## Add annotation information
  l <- addAnnoByVector(truncateEnsemblID(rownames(res)), species)
  
  EnhancedVolcano(res,
                  lab = l$symbol,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = title,
                  pCutoff = alpha,
                  FCcutoff = lfc_cutoff,
                  pointSize = 1.0,
                  labSize = 3.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75)
}

plot.topNVar.heatmap = function(vsd, topN) {
  ## Gene cluster
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), topN)
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vsd)[, c("mergeCond")])
  row.names(anno) <- colnames(vsd)
  pheatmap(mat, annotation_col = anno)
}


isEntrezID = function(v) {
  v2<-sapply(v, function(s) startsWith(s, 'ENS'))
  return(length(v2[v2==TRUE]) / length(v2) > 0.8)
}

truncateEnsemblID <- function(v) {
  r <- sapply(v, function(s) {t=strsplit(s,'\\.'); return(t[[1]][1])})
  return(r)
}

addAnno = function(res, species) {
  dbID = loadOrg.pathway(species)
  
  #ens.str <- substr(str_split(rownames(res), '.')[0], 1, 15)  
  ens.str<-truncateEnsemblID(rownames(res))
  
  l<-addAnnoByVector(ens.str, species)
  res$symbol <- l$symbol
  res$entrez <- l$entrez
  res$name <- l$name
  res$ensembl <- l$ensembl
  
  # if (isEntrezID(ens.str)) {
  #   res$symbol <- mapIds(get(dbID), keys = ens.str, column = "SYMBOL", keytype = "ENSEMBL", multiVals="first")
  #   res$entrez <- mapIds(get(dbID), keys = ens.str, column = "ENTREZID", keytype = "ENSEMBL", multiVals="first")
  #   res$name <- mapIds(get(dbID), keys = ens.str, column= "GENENAME", keytype = "ENSEMBL", multiVals = "first")
  #   res$ensembl <- mapIds(get(dbID), keys = ens.str, column= "ENSEMBL", keytype = "ENSEMBL", multiVals = "first")
  # } else {
  #   res$entrez <- mapIds(get(dbID), keys = ens.str,column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  #   res$symbol <- mapIds(get(dbID), keys = ens.str, column= "SYMBOL", keytype = "SYMBOL", multiVals = "first")
  #   res$name <- mapIds(get(dbID), keys = ens.str, column= "GENENAME", keytype = "SYMBOL", multiVals = "first")
  #   res$ensembl <- mapIds(get(dbID), keys = ens.str, column= "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  # }
  return(list(res,dbID))
}

addAnnoByVector = function(v, species) {
  dbID = loadOrg.pathway(species)
  r = list()
  if (isEntrezID(v)) {
    r$symbol <- mapIds(get(dbID), keys = v, column = "SYMBOL", keytype = "ENSEMBL", multiVals="first")
    r$entrez <- mapIds(get(dbID), keys = v, column = "ENTREZID", keytype = "ENSEMBL", multiVals="first")
    r$name <- mapIds(get(dbID), keys = v, column= "GENENAME", keytype = "ENSEMBL", multiVals = "first")
    r$ensembl <- mapIds(get(dbID), keys = v, column= "ENSEMBL", keytype = "ENSEMBL", multiVals = "first")
  } else {
    r$entrez <- mapIds(get(dbID), keys = v, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    r$symbol <- mapIds(get(dbID), keys = v, column= "SYMBOL", keytype = "SYMBOL", multiVals = "first")
    r$name <- mapIds(get(dbID), keys = v, column= "GENENAME", keytype = "SYMBOL", multiVals = "first")
    r$ensembl <- mapIds(get(dbID), keys = v, column= "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  }
  return(r)
}

create.lfcGeneList = function(res) {
  fc <- res$log2FoldChange
  names(fc) <- res$entrez
  fc<-na.omit(fc)
  fc = sort(fc, decreasing = TRUE)
  # head(fc)
}

loadOrg.pathway = function(species) {
  if (tolower(species) == 'human'){
    dbID = 'org.Hs.eg.db'
    library(org.Hs.eg.db)
  } 
  else if (tolower(species) == 'mouse') {
    dbID = 'org.Mm.eg.db'
    library(org.Mm.eg.db)
  }
  else {
    stop("unsupported species")
  }
  return(dbID)
}

loadKEGG.pathway = function(species){
  if (tolower(species) == 'human'){
    data(kegg.sets.hs)
    data(sigmet.idx.hs)
    return(kegg.sets.hs[sigmet.idx.hs])
  }
  else if (tolower(species) == 'mouse') {
    data(kegg.sets.mm)
    data(sigmet.idx.mm)
    return(kegg.sets.mm[sigmet.idx.mm])
  }
  else {
    stop("unsupported species")
  }
}

analysis.pathway = function(res, species, path, prefix) {
  # add annotation columns in res  
  l = addAnno(res, species)
  res=l[[1]]
  dbID=l[[2]]
  
  # generate named vector
  gene_list = create.lfcGeneList(res)
  
  # pathway analysis
  plot_pathway("KEGG", gene_list, species, dbID, path, prefix)
  plot_pathway("GO", gene_list, species, dbID, path, prefix)
}
