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

plot.countTransformation = function(dds, vsd, rld) {
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

calc.sampleDist <- function(obj, variables, delim) {
  sampleDists <- dist(t(assay(obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- pasteMultiCol(obj, variables, delim)
  colnames(sampleDistMatrix) <- NULL
  return(list(dist=sampleDists, dist.mat=sampleDistMatrix))
}

## Sample distance using VST
plot.sampleDist<- function(sampleDists, sampleDistMatrix, title="Sample distance heatmap"){
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors, fontsize=8,
           main=title)
}

plot.pca.transCount = function(obj, v, labelCol='Sample', title='PCA plot with transformed counts', showEcllipse=TRUE) {
  ## PCA plot
  pcaData <- plotPCA(obj, intgroup = v, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  meta<-data.frame(obj@colData@listData)
  meta<-meta[,-which(names(meta) %in% v)]
  pcaData <- merge(pcaData, meta, by=0, all.x=TRUE)
  p<-ggplot(pcaData, aes(x = PC1, y = PC2, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])),label=get(labelCol))) +
    geom_point(size =3) +
    geom_text_repel(hjust=0,vjust=0) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle(title) 
  
  if (showEcllipse) {
    p <- p + stat_ellipse()+stat_ellipse(level=0.8)
  }
  return(p)
}

plot.mds.transCount = function(obj, v, sampleDistMatrix, labelCol='Sample', title="MDS with transformed counts", showEcllipse=TRUE) {
  mds <- as.data.frame(colData(obj)) %>%
    cbind(cmdscale(sampleDistMatrix))
  ggplot(mds, aes(x = `1`, y = `2`, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])), label=get(labelCol))) +
    geom_point(size = 3) + geom_text_repel(hjust=0,vjust=0) +
    coord_fixed() + 
    ggtitle(title)
  
  if (showEcllipse) {
    p <- p + stat_ellipse()+stat_ellipse(level=0.8)
  }
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
    ggtitle("Generalized PCA of normalized counts")
}

# TS1: Lei@5/12
# This function is deprecated
# plot.mds = function(dge, v, labelCol='Sample') {
#   dge$samples$tmpGroup = pasteMultiCol(dge$samples, v, '_')
#   mds <- plotMDS(dge, pch="*", col=as.numeric(dge$samples$tmpGroup), plot=F)
#   names(mds)
#   data.frame(x=mds$x, y=mds$y, dge$samples) %>%
#     ggplot(aes(x=x, y=y, col=tmpGroup, label=get(labelCol)))+
#     geom_point() + geom_text_repel()+
#     ggtitle("MDS plot")
# }

plot.mds = function(dds, varaibles, delim="_", labelCol='Sample') {
  dds@colData$mergeCond<-pasteMultiCol(dds@colData, variables, delim)
  mds <- plotMDS(dds, pch="*", col=as.numeric(as.factor(dds@colData$mergeCond)), plot=F)
  names(mds)
  data.frame(x=mds$x, y=mds$y, dds@colData) %>%
    ggplot(aes(x=x, y=y, col=mergeCond, label=get(labelCol)))+
    geom_point() + geom_text_repel()+
    ggtitle("MDS plot of normalized counts")
}

plot.svaseq = function(dge, labelCol){
  cpm_new <- cpm(dge)
  mod  <- model.matrix(~ mergeCond, dge$samples)
  mod0 <- model.matrix(~   1, dge$samples)
  svseq <- svaseq(cpm_new, mod, mod0)
  df <- data.frame(svseq$sv,dge.filter$samples)
  ggplot(aes(x=X1, y=X2, col=Group, label=get(label_col)))+geom_label() + ggtitle("SVA plot")

  data.frame(svseq$sv, dge.filter$samples) %>%
    ggplot(aes(x=X1, y=X2, col=mergeCond, label=get(labelCol)))+geom_label() + ggtitle("SVA plot")
}


process.sampleVariance.all = function(cts.mat, sample.meta, variables, min_count=5, min_total_count=30, labelCol='Sample') {
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
  
  ## calculate sample distance matrix
  l_vst <- calc.sampleDist(vsd, variables, '_') 
  l_rld <- calc.sampleDist(rld, variables, '_')
  
  norm_count = counts(dds,normalized=TRUE)
  norm_count = data.frame(cbind(geneid=rownames(norm_count),norm_count))
  anno <- addAnnoByVector(truncateEnsemblID(rownames(norm_count)), species)
  norm_count$symbol = anno$symbol
  write.table(norm_count, 'norm_count.tsv',sep='\t',quote=F,row.names=F)
  ## plotting
  tryCatch(
    {
      print("Generating pdf")
      pdf(file='SampleVariance.all.pdf', width=20)
      p1<-plot.mds(dds, variables, '_', labelCol)
      p2<-plot.glmpca(dds, variables, labelCol)
      p3<-plot.countTransformation(dds, vsd, rld)
      
      test = l_vst$dist.mat
      colnames(test) = sample.meta$Sample
      plot.sampleDist(l_vst$dist, test, "Sample distance with VST transformed counts") # doesn't need to print it
      
      p5<-plot.mds.transCount(vsd, variables, l_vst$dist.mat, labelCol, "MDS plot with VST transformed counts", FALSE)
      p6<-plot.mds.transCount(rld, variables, l_rld$dist.mat, labelCol, "MDS plot with Rlog transformed counts", FALSE)
      
      p7<-plot.pca.transCount(vsd, variables, labelCol, "PCA plot with VST transformed counts", FALSE)
      p8<-plot.pca.transCount(rld, variables, labelCol, "PCA plot with Rlog transformed counts", FALSE)
      #p9<-plot.svaseq(dge.filter, labelCol)
      # plot.topNVar.vst(vsd, 20)
      print(p1+p2)
      print(p3)
      print(p5+p6)
      print(p7+p8)
      #print(p9)
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
