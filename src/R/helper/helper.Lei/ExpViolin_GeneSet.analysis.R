##
## Function: generating expression plots of geneset(s) in comparison with the rest of the genes in a comparison
## 


library(readxl)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyr)

inputDir='X:/projects/p046_KDM5/RNASeq.02012022/analysis.REdiscoverTE/Step2.2.Downstream_from_Limma/LFC.Gene'
outputDir='X:/projects/p046_KDM5/RNASeq.02012022/analysis.REdiscoverTE/Step2.2.Downstream_from_Limma/ExprViolin.Geneset'

dirs <- list.dirs(path=inputDir, full.names=TRUE, recursive=FALSE)

## load geneset excel file (each column is a gene list)
geneset_file='X:/projects/p046_KDM5/RNASeq.02012022/metadata/2022_IFN_signature.refined.xlsx'
df_geneset<-read_excel(geneset_file)


lapply(dirs, function(x) {
  ## get the comparison stat result file and load
  dirname<-basename(x)
  file=file.path(x, paste0(dirname,".result.csv"))
  print(file)
  
  df<-read.csv(file, stringsAsFactors=F, na.strings="unknown")
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  df$padj<-as.numeric(df$padj)
  df<-na.omit(df)

  for (i in 1:ncol(df_geneset)){
    geneset_name<-colnames(df_geneset)[i]
    geneset<- na.omit(df_geneset[[i]])
    mapped<-length(df[df$symbol %in% geneset, "symbol"])
    group_name<-paste0(geneset_name, " (", mapped,"/", length(geneset), " genes mapped)")
    df$Group <- ifelse(df$symbol %in% geneset, group_name, 'Others')

    img<-file.path(outputDir, paste(gsub('[\\. ]', '_', geneset_name), gsub('[\\. ]', '_', dirname), "exprViolin.png", sep="."))
    png(file=img)
    p<-ggplot(df, aes(x=Group, y=log2FoldChange, fill=Group)) + 
      geom_violin(trim=FALSE, alpha=0.3) + 
      geom_boxplot(width=0.2) + 
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      #geom_jitter(color="black", size=0.1, alpha=0.1) +
      theme_ipsum()+
      theme(legend.position = "top", axis.text.x = element_blank()) + 
      ggtitle(paste0(dirname, " geneset expression"))+
      ylim(-1,5)
    print(p)
    dev.off()
  }
})
