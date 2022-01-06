# Created on 2021/11/03 by Shangzhong.Li@pfizer.com

#!/usr/bin/env Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

# library(gtools)
if (!require("DESeq2")) {
  install.packages("DESeq2", dependencies = TRUE)
  library(DESeq2)
}
#====== 1. set working directory and list files ========
# filePath <- '/data/shangzhong/Contamination_Result'
args <- commandArgs(TRUE)
expr_mtx_fn <- args[1]
condition_fn <- args[2]
compare_fn <- args[3]
filePath <- args[4]
setwd(filePath)


# expr_mtx_fn <- '/hpc/grid/wip_drm_targetsciences/users/shangzhong/CTI/ILT3/merge_sub.tsv'
# condition_fn <- '/hpc/grid/wip_drm_targetsciences/users/shangzhong/CTI/ILT3/condition.txt'
# compare_fn <- '/home/lis262/Code/Scripts/HPC_Scripts/RNAseq/p01_DEseq2_pairs.tsv'
# filePath <- '/hpc/grid/wip_drm_targetsciences/users/shangzhong/CTI/ILT3/DESeq2'
setwd(filePath)

expr_mtx_df = read.table(expr_mtx_fn,row.names=1,header=T)
condition_df = read.table(condition_fn,header=T)
compare_df = read.table(compare_fn,header=T)

remove_samples = as.vector(condition_df[condition_df$condition == 'remove',]$SAMPLE_ID)
expr_mtx_df = expr_mtx_df[,!(names(expr_mtx_df) %in% remove_samples)]
condition_df = condition_df[!(condition_df$SAMPLE_ID %in% remove_samples),]
#====== 2. Do the differential expression ========
col <- colnames(compare_df)
row <- nrow(compare_df)

for (i in 1:row) {
    # 1st part, list the control sample and test sample
    ctrl = as.vector(compare_df[i,'control'])
    test = as.vector(compare_df[i,'test'])
    sub_cond_df = condition_df[condition_df$condition %in% c(ctrl,test),]
    row.names(sub_cond_df) = sub_cond_df$SAMPLE_ID
    # sub expression df
    sub_expr_df = expr_mtx_df[,row.names(sub_cond_df)]
    
    # 2nd part, do the DE for one compare
    dds <-DESeqDataSetFromMatrix(countData=sub_expr_df, colData=sub_cond_df,design=~condition)
    keep <- rowSums(counts(dds)) >=10
    dds <- dds[keep,]
    # define control and test
    dds$condition <- relevel(dds$condition, ref = ctrl)
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    result <- resOrdered[complete.cases(resOrdered),]
    sig_up <- result[result$padj < 0.05 & result$log2FoldChange >= log2(2),]
    sig_dn <- result[result$padj < 0.05 & result$log2FoldChange <= -log2(2),]
    print(paste("upreg",ctrl,test,nrow(sig_up),sep=" "))
    print(paste("dnreg",ctrl,test,nrow(sig_dn),sep=" "))
    outputFile <- paste(test,'_VS_',ctrl,'.csv',sep="")
    write.csv(result,outputFile,quote=F)
    # output cluster figure
    pdf(paste(strsplit(outputFile,'\\,')[[1]][1],'.pdf',sep=""))
    rld <- rlog(dds)
    # Heat map
    library("RColorBrewer")
    library("gplots")
    #heat map of distance matrix
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- with(colData(dds),sub_cond_df$condition)
    colours = colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
    heatmap.2(mat, trace="none", col = colours,margin=c(17,17))
    dev.off()
}
