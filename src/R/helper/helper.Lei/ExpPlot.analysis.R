library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(stringr)
library(GOfuncR)
library('fs')

source('src/R/util.load.cts.R')
source('src/R/util.verify.R')
source('src/R/util.sampleVariance.R')
source('src/R/util.DE.R')
source('src/R/util.pathway.R')

## Parameter setting
tpm.file='X:/projects/p036_LIRB1_2/Output/PreProcess/Gene_FPKM_37Samples.txt'
wd='X:/projects/p036_LIRB1_2/Output/CustomizedAnalysis/SelectedGeneExpressionPlot'
sample.meta.file = "X:/projects/p036_LIRB1_2/Output/PreProcess/SampleMeta_37Samples.txt"
variables<-c('Treatment')  ## Variables in design
species="human"

## Function to perform beeswarm and boxplot of expression for genes under interest
expBarplot<-function(TPM.mat, genes, species, sample.meta, prefix) {
  ## Formalize rownames with ensemble ID
  if (isEntrezID(rownames(TPM.mat))) {
    eid<-truncateEnsemblID(rownames(TPM.mat))
    dbID = loadOrg.pathway(species)
    gid<-mapIds(get(dbID),
                keys=eid, 
                column="SYMBOL",
                keytype="ENSEMBL",
                multiVals="first")
    rownames(TPM.mat)=gid
  }
  
  ## extract sub matrix, melt and merge sample metadata
  a<-TPM.mat[rownames(TPM.mat) %in% genes,]
  write.csv(a, file.path(wd, paste0(prefix,".TPM.csv")), quote=TRUE, row.names=TRUE, col.names=TRUE)
  b=melt(a, measure.vars=colnames(a))
  colnames(b)<-c("Gene", "Sample", "TPM")
  c=merge(x=b, y=sample.meta, by="Sample", all.x=TRUE)
  l=levels(c$Gene)
  l=l[order(l)]
  c=c %>% mutate(across(Gene, factor, levels=l)) ## order alphabetically
  
  ## plotting
  pdf(file.path(wd, paste0(prefix, ".TPM.beeswarm.pdf")), width=15, height=3*length(genes)/4)
  g1=ggplot(c, aes(x=mergeCond, y=TPM)) + 
    geom_beeswarm(cex = 3) + theme(axis.text.x=element_text(angle=45,vjust=0.1,hjust=NULL)) +
    facet_wrap(. ~ Gene, ncol=5, scale='free')
  print(g1)
  dev.off()
  
  pdf(file.path(wd, paste0(prefix, ".TPM.barplot.pdf")), width=15, height=3*length(genes)/4)
  g2=ggplot(c, aes(x=mergeCond, y=TPM)) + 
    geom_boxplot(aes(fill=Treatment), position=position_dodge(0.9)) + 
    scale_fill_viridis_d() + 
    theme(legend.position="top") + 
    facet_wrap(.~Gene, ncol=5,  scale='free')
  print(g2)
  dev.off()
}


## Load TPM matrix
tpm.mat = read.combined_cts(tpm.file)

## read sample metadata table
sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')

#################################################################
## Option #1 Read a input gene list and plot expression (TPM)
selected.gene.file='X:/projects/p036_LIRB1_2/Output/CustomizedAnalysis/Up.selectedGene.csv'
geneList<-read.csv(selected.gene.file, row.names=1)
expBarplot(tpm.mat, geneList$x, species, sample.meta, path_file(selected.gene.file))

#################################################################
## Option #2 Read a list of pathways under interest and plot expression
go_df=read.csv('X:/projects/p041_PMS1/Metadata/GOterms_for_profiling.txt', sep='\t', header=TRUE)
goIDs=go_df[[1]]

for (go in goIDs) {
  print(go)
  tryCatch ({
    genes=getGenesFromGO(c(go), species)$gene
    prefix=str_replace(go, ":", "_")
    expBarplot(tpm.mat, genes, species, sample.meta, prefix)
  },
  error = function(err) {
    print(paste("Error: ", err))
  }, 
  finally = function(f) {
  })
}


