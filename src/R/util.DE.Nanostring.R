suppressPackageStartupMessages({
  require(dplyr);
  require(ggplot2);
  require(gplots);
  require(DESeq2);
  require(stringr);
  require(RUVSeq);
  require(limma);
  require(matrixStats);
  require(MASS);
  require(RColorBrewer);
  require(PCAtools)
})


DE.DESeq2.NanoString = function(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, species) {
  # cts.mat: count matrix loaded using function read.combined_nano_cts
  # sample.meta: meta data, first column name needs to be "Sample".
  # compare_df: two columns, first is control, second is test
  # min_total_count: threshold to filter genes.
  # outdir: output directory
  # lfc_cutoff, alpha: log2foldchange and FDR threshold.
  # species: 'human' or 'mouse', used to annotate genes
  col <- colnames(compare_df)
  row <- nrow(compare_df)
  rownames(sample.meta) = sample.meta$Sample
  # remove low quality house keeping genes
  genes = rownames(cts.mat)[grep("^Endoge",rownames(cts.mat))]
  hk = rownames(cts.mat)[grep("^House",rownames(cts.mat))]
  condition = factor(sample.meta$mergeCond)
  samples = factor(sample.meta$Sample)
  set = newSeqExpressionSet(as.matrix(cts.mat), phenoData = sample.meta)
  hk_raw = cts.mat[grep("^House",rownames(cts.mat)),]
  pval = vector(length = nrow(hk_raw))
  for (i in 1:nrow(hk_raw)){
    reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(sample.meta$mergeCond)) # negative binomial 
    pval[i] = coef(summary(reg))[2,4]
  }
  sum(pval <= .05)
  # summary(reg)
  idx <- pval <= .05
  exc <- row.names(hk_raw[idx, ])
  hk = hk[!(hk %in% exc)] # remove low quality house keeping genes
  # RLE QC plot using RUVSeq to decide how many hidden factors to adjust.
  pdf(paste0(outdir,'/','RUVSeq_QC.pdf'))
  boxplot(log2(cts.mat+1))
  colors <- brewer.pal(6, "Set2")
  plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[sample.meta$group])
  plotPCA(set, col=colors[sample.meta$group], cex=1.2)
  # normalize the data using RUVSeq
  set1 = betweenLaneNormalization(set, which="upper")
  set1 = RUVg(set1, hk, k=1)
  dds = DESeqDataSetFromMatrix(counts(set1),colData=pData(set1),design=~1)
  dds = estimateSizeFactors(dds)
  dds = estimateDispersionsGeneEst(dds)
  cts <- counts(dds, normalized=TRUE)
  disp <- pmax((rowVars(cts) - rowMeans(cts)),0)/rowMeans(cts)^2
  mcols(dds)$dispGeneEst <- disp
  dds <- estimateDispersionsFit(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  mat <- assay(vsd)
  covars <- as.matrix(colData(dds)[,grep("W",colnames(colData(dds))),drop=FALSE])
  mat = removeBatchEffect(mat, covariates = covars)
  assay(vsd) = mat
  normalizedcount <- set1@assayData[["normalizedCounts"]]
  write.table(normalizedcount, paste0(outdir,'/RUVSeq_norm.tsv'), sep='\t', quote=F)
  
  dev.off()
  
  # DESeq analysis
  for (i in 1:row) {
    tryCatch({
      # 1st part, list the control sample and test sample
      ctrl = as.vector(compare_df[i,'control'])
      test = as.vector(compare_df[i,'test'])
      print(paste0('Process conparison ', ctrl, ' vs ', test))
      # 2. Do the DE analysis
      dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                                    colData = pData(set1),
                                    design = ~ W1 + condition)
      dds$condition <- relevel(dds$condition, ref = ctrl)
      dds <- DESeq(dds)
      contrast = c('condition',test,ctrl)
      result = as.data.frame(results(dds, contrast= contrast))
      
      # create sub-folder and save results
      comparison = str_replace_all(paste(test, ctrl, sep="_VS_"), '[ :,>]','_')
      path = file.path(outdir, comparison)
      dir.create(path)
      write.table(result, file.path(path, paste0(comparison, '.result.tsv')),sep='\t',quote=F)
      
      
      # output cluster figure
      pdf(file.path(path, paste0(comparison,'plots.pdf')))
      rld <- rlog(dds)
      rld.sub <- rld[ ,rld$condition %in% c(ctrl, test) ]
      #heat map of distance matrix
      distsRL <- dist(t(assay(rld.sub)))
      mat <- as.matrix(distsRL)
      sub_dds <- dds[, dds$condition %in% c(ctrl, test)]
      rownames(mat) <- colnames(mat) <- with(colData(sub_dds), sub_dds$condition)
      colours = colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
      heatmap.2(mat, trace="none", col = colours,margin=c(17,17))
      # print(p1)
      ## MA plot
      c <- resultsNames(sub_dds)[2]
      res <- lfcShrink(sub_dds, coef =  c, type="apeglm")
      DESeq2::plotMA(results(dds,contrast =contrast))
      # print(p2)
      ## Volcano plot
      result$symbol = str_split_fixed(row.names(result),"_", 3)[,2]
      p = plot.enhancedVolcano(result, species, lfc_cutoff, alpha, paste0(test, " vs ", ctrl));
      print(p)
      dev.off()
      
      ## Significant Up and Down
      sig_up <- result[result$padj < alpha & result$log2FoldChange >= log2(lfc_cutoff) & +
                         result$baseMean >= 20,]
      sig_dn <- result[result$padj < alpha & result$log2FoldChange <= -log2(lfc_cutoff) & +
                         result$baseMean >= 20,]
      print(paste("upreg",ctrl,test,nrow(sig_up),sep=" "))
      print(paste("dnreg",ctrl,test,nrow(sig_dn),sep=" "))
      
      ## pathway analysis
      path2 <- file.path(path,'Pathway')
      dir.create(path2)
      analysis.pathway.nano(sig_up, species, path2, 'Up')
      analysis.pathway.nano(sig_dn, species, path2, 'Down')
    },error = function(err) {
      print(paste("Error: ", err))
    }, 
    finally = function(f) {
      dev.off()
    }
    )
  }
}

plot.enhancedVolcano = function(res, species, lfc_cutoff, alpha, title) {
  ## Add annotation information
  l <- res
  
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

addAnnoNano = function(res, species) {
  dbID = loadOrg.pathway(species)
  
  #ens.str <- substr(str_split(rownames(res), '.')[0], 1, 15)  
  symbol = res$symbol
  
  l<-addAnnoByVectorNano(symbol, species)
  res$entrez <- l$entrez
  res$name <- l$name
  res$ensembl <- l$ensembl
  
  return(list(res,dbID))
}


addAnnoByVectorNano = function(v, species) {
  dbID = loadOrg.pathway(species)
  r = list()
  r$entrez <- mapIds(get(dbID), keys = v, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  r$name <- mapIds(get(dbID), keys = v, column= "GENENAME", keytype = "SYMBOL", multiVals = "first")
  r$ensembl <- mapIds(get(dbID), keys = v, column= "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  return(r)
}


analysis.pathway.nano = function(res, species, path, prefix) {
  # add annotation columns in res  
  l = addAnnoNano(res, species)
  res=l[[1]]
  dbID=l[[2]]
  
  # generate named vector
  gene_list = create.lfcGeneList(res)
  
  # pathway analysis
  plot_pathway("KEGG", gene_list, species, dbID, path, prefix)
  # plot_pathway("GO", gene_list, species, dbID, path, prefix)
}
