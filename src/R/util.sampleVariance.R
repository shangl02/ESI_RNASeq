suppressPackageStartupMessages({
  require(dplyr);
  require(ggplot2);
  require(DESeq2);
  require(glmpca);
  require(sva);
  require(pheatmap)
  require(RColorBrewer)
})

pasteMultiCol = function(vsd, variable, s) {
  delim=''
  r=vector("character")
  for (v in variables) {
    r = paste(r, vsd[[v]], sep=delim)
    delim=s
  }
  return(r)
}

plot.vst = function(dds, vsd, rld) {
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
  colnames(df)[1:2] <- c("x", "y")  
  lvls <- c("log2(x + 1)", "vst", "rlog")
  df$transformation <- factor(df$transformation, levels=lvls)
  ggplot(df, aes(x = x, y = y)) + 
    geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) +ggtitle("Comparison of un-transformed vs vst vs rlog")
}

calc.sampleDist.vst <- function(vsd, variables, delim) {
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- pasteMultiCol(vsd, variables, delim)
  colnames(sampleDistMatrix) <- NULL
  return(list(dist=sampleDists, dist.mat=sampleDistMatrix))
}

## Sample distance using VST
plot.sampleDist.vst<- function(sampleDists, sampleDistMatrix){
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors, fontsize=8)
}

plot.pca.vst = function(vsd, v, labelCol='Sample') {
  ## PCA plot
  pcaData <- plotPCA(vsd, intgroup = v, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  meta<-data.frame(vsd@colData@listData)
  meta<-meta[,-which(names(meta) %in% v)]
  pcaData <- merge(pcaData, meta, by=0, all.x=TRUE)
  ggplot(pcaData, aes(x = PC1, y = PC2, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])),label=get(labelCol))) +
    geom_point(size =3) +
    geom_text_repel(hjust=0,vjust=0) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle("PCA with VST data") + 
    stat_ellipse()+stat_ellipse(level=0.8)
}

plot.mds.vst = function(vsd, v, sampleDistMatrix, labelCol='Sample') {
  mds <- as.data.frame(colData(vsd)) %>%
    cbind(cmdscale(sampleDistMatrix))
  ggplot(mds, aes(x = `1`, y = `2`, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])), label=get(labelCol))) +
    geom_point(size = 3) + geom_text_repel(hjust=0,vjust=0) +
    coord_fixed() + 
    ggtitle("MDS with VST data")
}


plot.glmpca = function(dds, v, labelCol='Sample') {
  ## PCA plot using generalized PCA
  gpca <- glmpca(counts(dds), L=2)
  gpca.dat <- gpca$factors
  # for (variable in v) {
  #   gpca.dat[[variable]] <- dds[[variable]]
  # }
  gpca.dat <- merge(gpca.dat, data.frame(dds@colData@listData), by='row.names',all=TRUE)
  
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])), label=get(labelCol))) +
    geom_point() + geom_text_repel() +
    coord_fixed() +
    ggtitle("glmpca - Generalized PCA")
}

plot.mds = function(dge, v, labelCol='Sample') {
  dge$samples$tmpGroup = pasteMultiCol(dge$samples, v, '_')
  mds <- plotMDS(dge, pch="*", col=as.numeric(dge$samples$tmpGroup), plot=F)
  names(mds)
  data.frame(x=mds$x, y=mds$y, dge$samples) %>%
    ggplot(aes(x=x, y=y, col=tmpGroup, label=get(labelCol)))+
    geom_point() + geom_text_repel()+
    ggtitle("MDS plot")
}


plot.svaseq = function(dge, label_col='Sample'){
  cpm_new <- cpm(dge)
  mod  <- model.matrix(~ mergeCond, dge$samples)
  mod0 <- model.matrix(~   1, dge$samples)
  svseq <- svaseq(cpm_new, mod, mod0)
  df <- data.frame(svseq$sv, dge.filter$samples)
  if (showLabel) {
    ggplot(aes(x=X1, y=X2, col=Group, label=get(label_col)))+geom_label() + ggtitle("SVA plot")
  } else {
    ggplot(aes(x=X1, y=X2, col=Group, label=get(label_col)))+geom_label() + ggtitle("SVA plot")
    
  }
=======
  data.frame(svseq$sv, dge.filter$samples) %>%
    ggplot(aes(x=X1, y=X2, col=mergeCond, label=get(label_col)))+geom_label() + ggtitle("SVA plot")
}


process.sampleVariance.all = function(cts.mat, sample.meta, variables,min_count=5, min_total_count=30) {
  design.formula = build.formula(variables)
  
  ## build dge and dds
  dge <- DGEList(cts.mat, samples=sample.meta)
  design <- model.matrix(design.formula, data=dge$samples)
  keep <- filterByExpr(dge, design, min.count=min_count, min.total.count=min_total_count)
  dge.filter <- dge[keep, , keep.lib.sizes=F]
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(dge.filter$counts)),
                                colData = dge.filter$samples,
                                design = design.formula)
  
  ## Sample Variance Section
  vsd <- vst(dds, blind = FALSE)
  rld <- rlog(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds) #required for log2 approach which need size factors to account for sequencing depth, and specify normalized=TRUE
  l <- calc.sampleDist.vst(vsd, variables, '_') ## calculate sample distance matrix
  norm_count = counts(dds,normalized=TRUE)
  norm_count = data.frame(cbind(geneid=rownames(norm_count),norm_count))
  anno <- addAnnoByVector(truncateEnsemblID(rownames(norm_count)), species)
  norm_count$symbol = anno$symbol
  write.table(norm_count, 'norm_count.tsv',sep='\t',quote=F,row.names=F)
  ## plotting
  tryCatch(
    {
      print("Generating pdf")
      pdf(file='SampleVariance.all.pdf')
      p1<-plot.mds(dge.filter, variables)
      p2<-plot.glmpca(dds, variables)
      p3<-plot.vst(dds,vsd, rld)
      test = l$dist.mat
      colnames(test) = sample.meta$Sample
      plot.sampleDist.vst(l$dist, test) # doesn't need to print it
      p5<-plot.mds.vst(vsd, variables, l$dist.mat)
      p6<-plot.pca.vst(vsd, variables)
      p7<-plot.svaseq(dge.filter, 'Sample')
      # plot.topNVar.vst(vsd, 20)
      print(p1)
      print(p2)
      print(p3)
      print(p5)
      print(p6)
      print(p7)
      dev.off()
    }, 
    error = function(err) {
      print(paste("Error: ", err))
      dev.off()
    }, 
    finally =function(msg){
      print("in finally")
      dev.off()
    }
  )
}
