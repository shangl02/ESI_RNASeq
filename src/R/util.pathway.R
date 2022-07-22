suppressPackageStartupMessages({
  require(clusterProfiler);
  require(enrichplot);
  require(ggplot2);
  require(GOfuncR);
  require(org.Mm.eg.db);
  require(org.Hs.eg.db);
  require(org.Rn.eg.db);
  require(msigdbr)
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
    p_enrich = setReadable(p_enrich, dbID, keyType="ENTREZID")
  }
  else if (tolower(pathDB) == "go") {
    print("Pathway analysis with GO")
    p_enrich <- enrichGO(gene = names(gene_list), OrgDb = dbID, 
                         readable = T, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
    p_gse <- gseGO(geneList=gene_list, ont ="ALL", keyType = "ENTREZID", OrgDb = get(dbID),
                   pvalueCutoff = 0.05, verbose = TRUE, pAdjustMethod = "none")
  }
  
  print('Start writing')
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

plot.enrich <- function(file, pathway, dbID, gene_list) {
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



gsea = function(deseq2_fn, species, out_fn, plot_fn, category='C2') {
  ## this function runs gsea for deseq2 results
  # get species database
  dbID = loadOrg.pathway(species)
  # sort genes based on pi score
  df = read.csv(deseq2_fn, row.names=1)
  df$pi = -log10(df$pvalue) * sign(df$log2FoldChange)
  df = df[order(df$pi,decreasing = T),]
  anno <- addAnnoByVector(truncateEnsemblID(rownames(df)), species)
  df$entrez = anno$entrez
  
  gene_list = df$pi
  names(gene_list) = df$entrez
  # run GSEA
  m_t2g <- msigdbr(species=species, category=category) %>%
      dplyr::select(gs_name, entrez_gene)
  em <- GSEA(gene_list, TERM2GENE = m_t2g)
  em = setReadable(em, dbID, keyType="ENTREZID")  
  write.table(em@result, out_fn, sep="\t", quote=FALSE, row.names = FALSE)
  saveRDS(em, file=paste0(tools::file_path_sans_ext(out_fn),'.rds'))
  # plot results
  plot.enrich(plot_fn, em, dbID, gene_list)
}


gsea_plot = function(rds_fn, ids, out_pdf){
  # this function generates gsea plots for specific pathway
  pdf(out_pdf)
  em = readRDS(rds_fn)
  for (id in ids){
    tryCatch(
      {
        p <- gseaplot2(em, geneSetID = em$ID[id], title = em$Description[id])
        print(p)
      }, 
      error = function(err) {
        print(paste("Error in GSEA plot: ", err))
      }
    )
  }
  dev.off()
}


GSEA_wrapper = function(compare_df, sample.meta, deseq2_path, species){
  # this function run GSEA for all DESeq2 results automatically
  col <- colnames(compare_df)
  row <- nrow(compare_df)
  for (i in 1:row) {
    tryCatch({
      # 1st part, extract control and test sample metadata
      ctrl = as.vector(compare_df[i,'control'])
      test = as.vector(compare_df[i,'test'])
      print(paste0('Process comparison ', ctrl, ' vs ', test))
      sub_cond_df = sample.meta[sample.meta$mergeCond %in% c(ctrl,test),]
      row.names(sub_cond_df) = sub_cond_df$Sample
      
      # 2nd part, create sub-folder for GSEA
      comparison = str_replace_all(paste(test, ctrl, sep="_VS_"), '[ :,>]','_')
      path = glue('{deseq2_path}/{comparison}/GSEA')
      dir.create(path)
      deseq2_fn = glue('{deseq2_path}/{comparison}/{comparison}.result.csv')
      # 3rd part, run GSEA
      out_fn = glue('{path}/gsea.tsv')
      plot_fn = glue('{path}/gsea_plots.pdf')
      gsea(deseq2_fn, species, out_fn, plot_fn, category='C2')
    },
    error = function(err) {
      print(paste("Error:GSEA", err))
    },
    finally = function(f) {
    })    
  }
}
