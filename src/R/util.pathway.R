
suppressPackageStartupMessages({
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2);
})

plot_pathway = function(pathDB, gene_list, species, dbID, dir, prefix) {
  
  if (tolower(pathDB) == "kegg") {
    p_enrich <- enrichKEGG(gene = names(gene_list),
                           organism = tolower(species),
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.10)
    
    p_gse <- gseKEGG(geneList = gene_list, organism = species,
                     nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.1, verbose = FALSE)
  }
  else if (tolower(pathDB) == "go") {
    p_enrich <- enrichGO(gene = names(gene_list), OrgDb = dbID, 
                         readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
    
    p_gse <- gseGO(geneList=gene_list, ont ="ALL", keyType = "ENTREZID", OrgDb = get(dbID),
                   nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
  }
  
  # write GSE result
  write.table(p_gse@result, file.path(dir, paste0(prefix, ".", pathDB, '.GSE.result.txt')), sep="\t", quote=FALSE, row.names = FALSE)
  
  ## Save plots in pdf
  tryCatch(
    {
      pdf(file.path(dir, paste0(prefix, ".", pathDB, ".plots.pdf")))
      p1<-dotplot(p_enrich)
      
      x <- pairwise_termsim(p_enrich)
      p2<-emapplot(x, showCategory = 10)
      
      p_enrich <- setReadable(p_enrich, dbID, keyType="ENTREZID")
      p3 <- cnetplot(p_enrich, categorySize="pvalue", foldChange=gene_list, showCategory = 10, shadowtext='category',fixed=FALSE)
      
      p4<-ridgeplot(p_gse) + labs(x = "enrichment distribution")
      
      print(p1)
      print(p2)
      print(p3)
      print(p4)
      dev.off()
    }, 
    error = function(err) {
      print(paste("Error: ", err))
    },
    finally = function(){
      dev.off()
    }
  )
}

