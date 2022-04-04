library(reshape2)
library(ggplot2)
library(ggbeeswarm)
library(stringr)
library(GOfuncR)
library('fs')

setwd('C:/Users/shangl02/source/repos/ESI_RNASeq/')
source('src/R/util.load.cts.R')
source('src/R/util.verify.R')
source('src/R/util.sampleVariance.R')
source('src/R/util.DE.R')
source('src/R/util.pathway.R')

param_file='X:\\projects\\p046_KDM5\\RNASeq.02012022\\analysis.REdiscoverTE\\param.R'
source(param_file)
tpm.file=cts_file
wd=outdir


## Parameter setting
# tpm.file='X:/projects/p036_LIRB1_2/Output/PreProcess/Gene_FPKM_37Samples.txt'
# wd='X:/projects/p036_LIRB1_2/Output/CustomizedAnalysis/SelectedGeneExpressionPlot'
# sample.meta.file = "X:/projects/p036_LIRB1_2/Output/PreProcess/SampleMeta_37Samples.txt"
# variables<-c('Treatment')  ## Variables in design
# species="human"

## Function to perform beeswarm and boxplot of expression for genes under interest
expBarplot<-function(tpm.mat, genes, species, sample.meta, outdir, prefix) {
  ## Formalize rownames with ensemble ID
  if (isEntrezID(rownames(tpm.mat))) {
    eid<-truncateEnsemblID(rownames(tpm.mat))
    dbID = loadOrg.pathway(species)
    gid<-mapIds(get(dbID),
                keys=eid, 
                column="SYMBOL",
                keytype="ENSEMBL",
                multiVals="first")
    rownames(tpm.mat)=gid
  }
  
  ## add a human housekeeping geen if only one gene in the list
  if (length(genes)==1){
    genes=c(genes,'ACTB')  
  }
  
  ## extract sub matrix, melt and merge sample metadata
  a<-tpm.mat[rownames(tpm.mat) %in% genes, ]
  
  write.csv(a, file.path(outdir, paste0(prefix,".selectedGene.Expression.csv")), quote=TRUE, row.names=TRUE)
  b=melt(a, measure.vars=colnames(a))
  colnames(b)<-c("Gene", "Sample", "TPM")
  c=merge(x=b, y=sample.meta, by="Sample", all.x=TRUE)
  l=levels(c$Gene)
  l=l[order(l)]
  c=c %>% mutate(across(Gene, factor, levels=l)) ## order alphabetically
  
  wUnit=8
  wNum=3
  hUnit=4
  plotWidth=ifelse(length(genes)<wNum, wUnit*length(genes), wUnit*wNum)
  plotHeight=hUnit*(round(length(genes)/wNum)+1)
  ## plotting
  pdf(file.path(outdir, paste0(prefix, ".selectedGene.Expression.beeswarm.pdf")), width=plotWidth, height=plotHeight)
  g1=ggplot(c, aes(x=mergeCond, y=TPM)) + 
    geom_beeswarm(cex = 1) + theme(axis.text.x=element_text(angle=90,vjust=0.3,hjust=0.3)) +
    facet_wrap(. ~ Gene, ncol=wNum, scale='free') + 
    ggtitle(paste0(prefix,' Expression'))
  print(g1)
  dev.off()
  
  pdf(file.path(outdir, paste0(prefix, ".selectedGene.Expression.barplot.pdf")), width=plotWidth, height=plotHeight)
  g2=ggplot(c, aes(x=mergeCond, y=TPM)) + 
    geom_boxplot(aes(fill=mergeCond), position=position_dodge(0.9)) + 
    scale_fill_viridis_d() + 
    theme(legend.position="top", axis.text.x=element_text(angle=90,vjust=0.3,hjust=0.3)) + 
    facet_wrap(.~Gene, ncol=wNum,  scale='free')+
    ggtitle(paste0(prefix,' Expression'))
  
  print(g2)
  dev.off()
  
  pdf(file.path(outdir, paste0(prefix, ".selectedGene.Expression.violin.pdf")), width=plotWidth, height=plotHeight)
  p <- ggplot(c, aes(x=mergeCond, y=TPM, fill=mergeCond)) + 
    geom_violin(trim=FALSE) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
    #geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme(legend.position="top", axis.text.x=element_text(angle=90,vjust=0.3,hjust=0.3)) + 
    facet_wrap(.~Gene, ncol=wNum,  scale='free') + 
    ggtitle(paste0(prefix,' Expression'))
  print(p)
  dev.off()
}



## Load TPM matrix
tpm.mat = read.combined_cts(tpm.file)

## read sample metadata table
sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')

#################################################################
## Option #1 Read a input gene list and plot expression (TPM)
selected.gene.file='X:/projects/p046_KDM5/RNASeq.02012022/metadata/GeneList.Guler2017.txt'
#geneList<-read.csv(selected.gene.file, row.names=1)  ## Two column file
#geneList<-geneList[[1]]

geneList<-scan(selected.gene.file, character(), quote = "")

expBarplot(tpm.mat, geneList, species, sample.meta, outdir, prefix)

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


