
suppressPackageStartupMessages({
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2);
  require(GOfuncR);
})

plot_pathway = function(pathDB, gene_list, species, dbID, dir, prefix) {
  
  if (tolower(pathDB) == "kegg") {
    ## enrichment analysis
    print("Pathway analysis with KEGG")
    p_enrich <- enrichKEGG(gene = names(gene_list),
                           organism = tolower(species),
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.10)
    
    ## GSEA
    p_gse <- gseKEGG(geneList = gene_list, organism = species,
                     nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.1, verbose = FALSE)
  }
  else if (tolower(pathDB) == "go") {
    print("Pathway analysis with GO")
    p_enrich <- enrichGO(gene = names(gene_list), OrgDb = dbID, 
                         readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
    
    p_gse <- gseGO(geneList=gene_list, ont ="ALL", keyType = "ENTREZID", OrgDb = get(dbID),
                   nPerm = 10000, minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
  }
  
  # write Enrichment result
  write.table(p_enrich@result, file.path(dir, paste0(prefix, ".", pathDB, '.Enrichment.result.txt')), sep="\t", quote=FALSE, row.names = FALSE)
  
  # write GSE result
  write.table(p_gse@result, file.path(dir, paste0(prefix, ".", pathDB, '.GSE.result.txt')), sep="\t", quote=FALSE, row.names = FALSE)
  
  p_enrich@result <- jamba::renameColumn(p_enrich@result,
                                             from="description",
                                             to="Description")
  
  p_gse@result <- jamba::renameColumn(p_gse@result, 
                                      from="description",
                                             to="Description")
  
  
  ## Save plots in pdf
  file=file.path(dir, paste0(prefix, ".", pathDB, ".enrichment.plots.PVal0.05_QVal0.1.pdf"))
  print("Plot pathway with enrichment")
  plot.enrich(file, p_enrich, dbID, gene_list)
  
  file=file.path(dir, paste0(prefix, ".", pathDB, ".GSE.plots.PVal0.05_QVal0.1.pdf"))
  print("plot pathway with GSE")
  plot.enrich(file, p_gse, dbID, gene_list)
  # 
  # pdf(file.path(dir, paste0(prefix, ".", pathDB, ".enrichment.plots.PVal0.05_QVal0.1.pdf")))
  # tryCatch(
  #   {
  #     p1<-dotplot(p_enrich)
  #     print(p1)
  #   }, 
  #   error = function(err) {
  #     print(paste("Error: ", err))
  #   })
  # tryCatch(
  #   {
  #     x <- pairwise_termsim(p_enrich)
  #     p2<-emapplot(x, showCategory = 10)
  #     
  #     p_enrich <- setReadable(p_enrich, dbID, keyType="ENTREZID")
  #     print(p2)
  #     
  #   },
  #   error = function(err) {
  #     print(paste("Error: ", err))
  #   })
  # tryCatch(
  #   {
  #     p3 <- cnetplot(p_enrich, categorySize="pvalue", foldChange=gene_list, showCategory = 10, shadowtext='category',fixed=FALSE)
  #     print(p3)
  #   },  
  #   error = function(err) {
  #     print(paste("Error: ", err))
  #   })
  # tryCatch(
  #   {
  #     p4<-ridgeplot(p_gse) + labs(x = "Gene set Enrichment Distribution")
  #     print(p4)
  #   }, 
  #   error = function(err) {
  #     print(paste("Error: ", err))
  #   }
  # )
  # dev.off()
}

plot.enrich<-function(file, pathway, dbID, gene_list) {
  ## Save plots in pdf
  pdf(file)
  tryCatch(
    {
      p1<-dotplot(pathway)
      print(p1)
    }, 
    error = function(err) {
      print(paste("Error in dotplot: ", err))
    })
  tryCatch(
    {
      x <- pairwise_termsim(pathway)
      p2<-emapplot(x, showCategory = 10)
      
      pathway <- setReadable(pathway, dbID, keyType="ENTREZID")
      print(p2)
      
    },
    error = function(err) {
      print(paste("Error in emapplot: ", err))
    })
  tryCatch(
    {
      p3 <- cnetplot(pathway, categorySize="pvalue", foldChange=gene_list, showCategory = 10, shadowtext='category',fixed=FALSE)
      print(p3)
    },  
    error = function(err) {
      print(paste("Error in cnetplot: ", err))
    })
  tryCatch(
    {
      p4<-ridgeplot(pathway) # + labs(x = "Gene set Enrichment Distribution")
      print(p4)
    }, 
    error = function(err) {
      print(paste("Error in ridgeplot: ", err))
    }
  )
  dev.off()
}


getGenesFromGO = function(goIDs, species) {
    if (tolower(species) == 'human'){
      dbID = 'Human.sapiens'
    } 
    else if (tolower(species) == 'mouse') {
      dbID = 'Mus.musculus'
    }
    else {
      stop("unsupported species")
    }
  
  get_anno_genes(goIDs, database=dbID)
}

