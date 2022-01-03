suppressPackageStartupMessages({
  require(dplyr);
  require(ggplot2);
  require(DESeq2);
  require(glmpca);
  require(sva);
})

pasteMultiCol = function(vsd, variables, s) {
  delim = ''
  r=vector("character")
  for (v in variables) {
    r = paste(r, vsd[[v]], sep=delim)
    delim=s
  }
  return(r)
}

plot.vst = function(vsd, rld) {
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  
  colnames(df)[1:2] <- c("x", "y")  
  lvls <- c("log2(x + 1)", "vst", "rlog")
  df$transformation <- factor(df$transformation, levels=lvls)
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) +ggtitle("Comparison of un-transformed vs vst vs rlog")
}

## Sample distance using VST
plot.sampleDist.vst<- function(vsd, variables, s){
  sampleDists <- dist(t(assay(vsd)))
  library("pheatmap")
  library("RColorBrewer")
  sampleDistsMatrix <- as.matrix(sampleDists)
  rownames(sampleDistsMatrix) <- pasteMultiCol(vsd, variables, s)
  colnames(sampleDistsMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistsMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
}

plot.pca.vst = function(vsd, v) {
  ## PCA plot
  pcaData <- plotPCA(vsd, intgroup = v, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(x = PC1, y = PC2, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])))) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    ggtitle("PCA with VST data")
}

plot.mds.vst = function(vsd, v) {
  mds <- as.data.frame(colData(vsd)) %>%
      cbind(cmdscale(sampleDistsMatrix))
    ggplot(mds, aes(x = `1`, y = `2`, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1])))) +
      geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
}

plot.glmpca = function(dds, v) {
## PCA plot using generalized PCA
  gpca <- glmpca(counts(dds), L=2)
  gpca.dat <- gpca$factors
  for (variable in v) {
    gpca.dat[[variable]] <- dds[[variable]]
  }
  ggplot(gpca.dat, aes(x = dim1, y = dim2, color = get(v[1]), shape = get(ifelse(length(v)>1, v[2], v[1]))))+
    geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
}

plot.mds = function(dge, v) {
  dge$samples$tmpGroup = pasteMultiCol(dge$samples, v, '_')
  mds <- plotMDS(dge, pch="*", col=as.numeric(dge$samples$tmpGroup), plot=F)
  names(mds)
  data.frame(x=mds$x, y=mds$y, dge$samples)%>%
    ggplot(aes(x=x, y=y, col=tmpGroup))+geom_point() + ggtitle("MDS plot")
}


plot.svaseq = function(dge.filter, label_col){
  cpm_new <- cpm(dge.filter)
  mod  <- model.matrix(~ Group, dge.filter$samples)
  mod0 <- model.matrix(~   1, dge.filter$samples)
  svseq <- svaseq(cpm_new, mod, mod0)
  data.frame(svseq$sv, dge.filter$samples) %>%
    ggplot(aes(x=X1, y=X2, col=Group, label=get(label_col)))+geom_label() + ggtitle("SVA plot")
}

