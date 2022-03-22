# this file defines softwares that need to be installed
packages = c('dplyr','ggplot2','DESeq2','glmpca','sva','pheatmap',
              'RColorBrewer','EnhancedVolcano','org.Mm.eg.db',
             'org.Hs.eg.db','PCAtools','RUVSeq','gplots','stringr',
             'limma','matrixStats','MASS','pheatmap','viridis','GOfuncR')
BiocManager::install(packages)


