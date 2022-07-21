# this file defines softwares that need to be installed
packages = c('dplyr','ggplot2','DESeq2','glmpca','sva','pheatmap',
             'RColorBrewer','EnhancedVolcano','org.Mm.eg.db',
             'org.Hs.eg.db','org.Rn.eg.db','PCAtools','RUVSeq','gplots','stringr',
             'limma','matrixStats','MASS','pheatmap','viridis','GOfuncR',
             'clusterProfiler','apeglm','gage','gageData','forcats',
             'ggnewscale','ggridges','msigdbr','ggupset')

BiocManager::install(packages)