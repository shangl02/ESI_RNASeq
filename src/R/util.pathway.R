
suppressPackageStartupMessages({
  require(clusterProfiler);
  require(enrichplot);
  require(ggplot2);
  require(GOfuncR);
  require(org.Mm.eg.db);
  require(org.Hs.eg.db)
})

plot_pathway = function(pathDB, gene_list, species, dbID, dir, prefix) {
  gene_list = gene_list[unique(names(gene_list))] # get unique genes
  gene_list = sort(gene_list, decreasing=T)
  if (tolower(pathDB) == "kegg") {
    ## enrichment analysis
    print("Pathway analysis with KEGG")
    p_enrich <- enrichKEGG(gene = names(gene_list),
                           organism = tolower(species),
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.10)
    
    ## GSEA
    p_gse <- gseKEGG(geneList = gene_list, organism = species,
                      pvalueCutoff = 0.1, verbose = FALSE)
    # set readable
    if (species == 'mouse'){
      p_enrich = setReadable(p_enrich, org.Mm.eg.db, keyType="ENTREZID")
    }
    if (species == 'human'){
      p_enrich = setReadable(p_enrich, org.Hs.eg.db, keyType="ENTREZID")
    }
  }
  else if (tolower(pathDB) == "go") {
    print("Pathway analysis with GO")
    p_enrich <- enrichGO(gene = names(gene_list), OrgDb = dbID, 
                         readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
    p_gse <- gseGO(geneList=gene_list, ont ="ALL", keyType = "ENTREZID", OrgDb = get(dbID),
                   pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
  }
  
  # write Enrichment result
  write.table(p_enrich@result, file.path(dir, paste0(prefix, ".", pathDB, '.Enrichment.result.txt')), sep="\t", quote=FALSE, row.names = FALSE)
  
  # write GSE result
  write.table(p_gse@result, file.path(dir, paste0(prefix, ".", pathDB, '.GSE.result.txt')), sep="\t", quote=FALSE, row.names = FALSE)

  
  ## Save plots in pdf
  file=file.path(dir, paste0(prefix, ".", pathDB, ".enrichment.plots.PVal0.05_QVal0.1.pdf"))
  print("Plot pathway with enrichment")
  plot.enrich(file, p_enrich, dbID, gene_list)
  
  file=file.path(dir, paste0(prefix, ".", pathDB, ".GSE.plots.PVal0.05_QVal0.1.pdf"))
  print("plot pathway with GSE")
  plot.enrich(file, p_gse, dbID, gene_list)
  
}

plot.enrich<-function(file, pathway, dbID, gene_list) {
  write.table(p_enrich@result, file.path(dir, paste0(prefix, ".", pathDB, '.enrich.result.txt')), sep="\t", quote=FALSE, row.names = FALSE)
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

