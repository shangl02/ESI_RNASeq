library("ggVennDiagram")

## IO setting
rootDir='X:/projects/p036_LIRB1_2/Output/Round2_top2000'
outdir='X:/projects/p036_LIRB1_2/Output/CustomizedAnalysis'

## Species
species="human"

## Define direction and wanted treatment/control
direction="Up"
treat=c('1880_RN888_VS_antiHA', '1882_RN888_VS_antiHA')
ctrl=c('1E1_RN888_VS_antiHA', '741_RN888_VS_antiHA', 'antiHA_RN888_VS_antiHA')

## Cutoff of significance
lfc_cutoff=2
alpha=0.05
#topN=500


## function to load treatment_vs_ctrl DE result csv
read.result<-function(rootDir, list) {
  l=list()
  for (t in list) {
    f=file.path(rootDir,t,paste0(t,".result.csv"))
    df=read.csv(f, row.names=1, header=TRUE)
    l[[length(l)+1]] = df
  }
  return(l)
}

## function to filter DE result table
table.filter<-function(df, direction, lfc_cutoff, alpha) {
  if (tolower(direction) == "up") {
    df_sub = df[df$padj < alpha & df$log2FoldChange >= lfc_cutoff,]
  } 
  else if (tolower(direction) == "down") {
    df_sub = df[df$padj < alpha & df$log2FoldChange <= -lfc_cutoff, ]
  }
  else {
    stop(paste0("Direction not recognized: ", direction))  
  }
  return(df_sub)
}

## function to intersect/union a list of vectors 
merge.gene=function(list, operation) {
  v=c()
  i=0
  for (l in list) {
    v2=l$symbol
    if (i==0){
      v=v2
    }
    else {
      if (tolower(operation) == "intersect") {
        v=intersect(v, v2)
      } else if (tolower(operation) == "union") {
        v=union(v, v2)
      } else {
        stop(paste0("Unrecognizable operation ", operation))
      }
    }
    i=i+1
  }
  return(v)
}

## function to merge multiple DE result tables
merge.table=function(df_list, name_list) {
  stopifnot(length(df_list)==length(name_list))
  
  df_all=data.frame()
  for (i in 1:length(name_list)){
    l=df_list[[i]]
    nc=lapply(colnames(l), function(x) paste(name_list[i], x, sep="."))
    colnames(l)=nc
    
    row.names(l)=l$bycol
    if (i==1){
      df_all=l
    }
    else{
      
      df_all=merge(df_all, l, by='row.names', all=TRUE)
    }
  }
  return(df_all)
}

## load all result files
l_treat.all=read.result(rootDir, treat)
l_ctrl.all=read.result(rootDir, ctrl)

l_treat=lapply(l_treat.all, function(x) table.filter(x, direction, lfc_cutoff, alpha))
l_ctrl=lapply(l_ctrl.all, function(x) table.filter(x, direction, lfc_cutoff, alpha))
dim(l_treat[[1]])
dim(l_ctrl[[1]])

## CHANGE ACCORDINGLY
## intersect the candidates in treatment
g_treat=merge.gene(l_treat, "intersect")

## CHANGE ACCORDINGLY
## union the candidates in control 
g_ctrl=merge.gene(l_ctrl, "union")

## Diff genes between treatment and ctrl groups
diff=setdiff(g_treat, g_ctrl)
length(diff)
write.csv(diff, file.path(outdir, paste0(direction,'.selectedGene.csv')))



suppressPackageStartupMessages({
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2);
  require(GOfuncR);
})
dbID = loadOrg.pathway(species)
entrez <- mapIds(get(dbID), keys = diff, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
kegg_enrich <- enrichKEGG(gene = entrez,
                       organism = tolower(species),
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.10)
write.table(kegg_enrich@result, file.path(outdir, paste0(direction,'.KEGG_Enrichment.txt')), sep="\t", quote=FALSE, row.names = FALSE)


go_enrich <- enrichGO(gene = entrez, OrgDb = dbID, 
                     readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
write.table(go_enrich@result, file.path(outdir, paste0(direction,'.GO_Enrichment.txt')), sep="\t", quote=FALSE, row.names = FALSE)




## merging all DE results
df_all=merge.table(c(l_treat.all, l_ctrl.all), 
                   c(treat, ctrl))

df_diff=df_all[df_all$`1880_RN888_VS_antiHA.symbol` %in% diff,]
write.csv(df_diff, file.path(outdir, paste0(direction,'selectedGene.DEResult.csv')))
# 
# 
# df_all=data.frame()
# for (i in 1:length(treat)){
#   ##
#   l=l_treat[[i]]
#   nc=lapply(colnames(l), function(x) paste(treat[i], x, sep="."))
#   colnames(l)=nc
#   
#   if (i==1){
#     df_all=l
#   }
#   else{
#     df_all=merge(df_all, l, by='row.names', all=TRUE)
#   }
# }
# 
# 
# 
# for (l in l_treat)
# {
#   
#   l_all[[length(l_all)+1]] = rownames(l)
# }
# for (l in l_ctrl)
# {
#   l_all[[length(l_all)+1]] = rownames(l)
# }
# #names(l_all) = c(treat, ctrl)
# 
# ggVennDiagram(l_all) +
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
# 
# 
# # 
# # 
# # 
# # 
# # t<-l_treat[[1]] df1$symbol %in% diff,]
# # t<-t[order(t$padj),]
# # n<-min(c(length(t[[1]]), topN))
# # t_final<-t[1:n,]
# # analysis.pathway(t_final, species, path2, direction)
# # 
# 
# 
# 
# 
# 
# 
# ggVennDiagram(list(HO=df1$symbol,WT=df_ctrl$symbol)) +
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
# 
# 
# 
# ## tag of this analysis
# tag=paste0(treat, "_except_", ctrl)
# 
# ## load significant genes in treat and ctrl
# file1=file.path(rootDir, treat, "Cerebellum_M_VS_Striatum_M", paste0(direction, "Gene_lfc1_PValue0.05.txt"))
# file2=file.path(rootDir, ctrl, "Cerebellum_M_VS_Striatum_M", paste0(direction, "Gene_lfc1_PValue0.05.txt"))
# 
# df1=read.csv(file1, sep=",", row.names = 1, header=TRUE)
# df_ctrl=read.csv(file2,sep=",", row.names = 1, header=TRUE)
# 
# ## venn diagram
# pdf(file.path(outdir, paste0(direction,"Gene.", tag, ".venn.pdf")))
# ggVennDiagram(list(HO=df1$symbol,WT=df_ctrl$symbol)) +
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
# dev.off()
# 
# diff=setdiff(df1$symbol, df_ctrl$symbol)
# 
# t<-df1[df1$symbol %in% diff,]
# t<-t[order(t$padj),]
# n<-min(c(length(t[[1]]), topN))
# t_final<-t[1:n,]
# write.csv(t_final, file.path(outdir, paste0(direction,"Gene.", tag,".top",topN,".csv")), row.names=TRUE, col.names=TRUE, quote=FALSE)
# 
# path2 <- file.path(outdir, paste0('Pathway.', tag))
# dir.create(path2)
# analysis.pathway(t_final, species, path2, direction)
