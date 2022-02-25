source('src/R/util.load.cts.R')
source('src/R/util.load.gtf.R')
source('src/R/util.verify.R')
source('src/R/util.design.R')
source('src/R/util.sampleVariance.R')
source('src/R/util.DE.R')
source('src/R/util.pathway.R')
source('src/R/util.pheatmap.R')
source('src/R/util.expression.boxplot.R')
## Load parameter◘
source('src/R/Helper/local_param.R')

process.sampleVariance.all = function(cts.mat, sample.meta, variables,min_count=5, min_total_count=30) {
  design.formula = build.formula(variables)
  
  ## build dge and dds
  dge <- DGEList(cts.mat, samples=sample.meta)
  design <- model.matrix(design.formula, data=dge$samples)
  keep <- filterByExpr(dge, design, min.count=min_count, min.total.count=min_total_count)
  dge.filter <- dge[keep, , keep.lib.sizes=F]
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(dge.filter$counts)),
                                colData = dge.filter$samples,
                                design = design.formula)
  
  ## Sample Variance Section
  vsd <- vst(dds, blind = FALSE)
  # rld <- rlog(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds) #required for log2 approach which need size factors to account for sequencing depth, and specify normalized=TRUE
  l <- calc.sampleDist.vst(vsd, variables, '_') ## calculate sample distance matrix
  norm_count = counts(dds,normalized=TRUE)
  norm_count = data.frame(cbind(geneid=rownames(norm_count),norm_count))
  anno <- addAnnoByVector(truncateEnsemblID(rownames(norm_count)), species)
  norm_count$symbol = anno$symbol
  write.table(norm_count, 'norm_count.tsv',sep='\t',quote=F,row.names=F)
  ## plotting
  tryCatch(
    {
      print("Generating pdf")
      pdf(file='SampleVariance.all.pdf')
      # p1<-plot.mds(dge.filter, variables)
      # p2<-plot.glmpca(dds, variables)
      # p3<-plot.vst(dds,vsd, rld)
      # test = l$dist.mat
      # colnames(test) = sample.meta$Sample
      # p4<-plot.sampleDist.vst(l$dist, test)
      # p5<-plot.mds.vst(vsd, variables, l$dist.mat)
      p6<-plot.pca.vst(vsd, variables)
      # p7<-plot.svaseq(dge.filter, 'Sample')
      # # plot.topNVar.vst(vsd, 20)
      # print(p1)
      # print(p2)
      # print(p3)
      # print(p4)
      # print(p5)
      print(p6)
      # print(p7)
      dev.off()
    }, 
    error = function(err) {
      print(paste("Error: ", err))
      dev.off()
    }, 
    finally =function(msg){
      print("in finally")
      dev.off()
    }
  )
}

## Load counts
## Option 1: load Sample counts per file
cts.mat = read_merge.cts(pattern, metric.vec, patName)
# set sample names

## Option 2: load count matrix
cts.mat = read.combined_cts(cts_file)
dim(cts.mat)

## Option 3: load nanostring count matrix
cts.mat = read.csv(cts_file, sep='\t')


## Output directory
dir.create(outdir)
setwd(outdir)

## Load sample metadata
sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')

## Verify
if (!verify.meta(cts.mat, sample.meta)) {
  print("Ordering columns of count matrix")
  cts.mat<-cts.mat[,unlist(sample.meta[1])]
  stopifnot(verify.meta(cts.mat.ordered, sample.meta))
}

# Generate sample variance report
process.sampleVariance.all(cts.mat, sample.meta, variables)

## Read comparison file
compare_df = read.table(comparison.file,header=T, sep="\t")

## DE analysis
DE.DESeq(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, species)


#=============================================
#        pheatmap for subset of genes
#=============================================
gene_fn = '/media/EphB2/genes.txt'
outFig = '/media/EphB2/pheatmap.png'
norm_cts_fn = '/media/EphB2//DESeq2/norm_count_sub.tsv'
sample_meta_fn = '/media/EphB2/condition_sub.tsv'
variables = c('merge')
logTrans = T
plot.pheatmap(norm_cts_fn, sample.meta.file, variables, gene_fn, outFig, logTrans)

#=============================================
#        pathway analysis for genes
#=============================================
gene_fn = '/media/EphB2/genes.txt'
genes = read.table(gene_fn)$V1
res_fn = '/media/EphB2/DESeq2/2A_Cre_EphB2_sg_VS_2B_Cre_SC_sg/2A_Cre_EphB2_sg_VS_2B_Cre_SC_sg.result.csv'
res = read.csv(res_fn, row.names=1)
res = res[rownames(res) %in% genes,]
path = '/media/EphB2/'
prefix = 'pathway'
species = 'mouse'
analysis.pathway(res, species, path, prefix)


#=============================================
#            merge all DESeq results
#=============================================
fn_paths = '/media/EphB2/DESeq2/'
fns  = Sys.glob(file.path(fn_paths,'*','*.result.csv'))
dfs = lapply(fns, read.csv)
res =  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "X", all.x = TRUE),
              dfs)
write.table(res,'/media/CD38/DESeq2/merge_res.tsv',sep='\t',quote=F)
df = read.table('/media/CD38/DESeq2/merge_lfc_qval.tsv',header=T)
sub_df = df[df$symbol %in% genes,]
sub_df = sub_df[order(sub_df$symbol),]
write.table(sub_df, '/media/CD38/DESeq2/sub_lfc_qval.tsv',sep='\t', 
            quote=F,row.names=F)

#================================================
#           gene expression boxplot
#================================================
norm_ctx_fn = '/media/CD38/tpm.tsv'
gene_fn = '/media/CD38/genes.txt'
sample_meta_fn = '/media/CD38/condition.tsv'
figure = '/media/CD38//genes_boxplot.png'
species = 'mouse'
logTrans = F
variables = c('Strain','Diet')
expression_boxplot(norm_ctx_fn, gene_fn, sample_meta_fn, variables, species, figure, logTrans)

#================================================
#           Enhanced volcano plot
#================================================
fn = '/media/EphB2/DESeq2/2A_Cre_EphB2_sg_VS_2B_Cre_SC_sg/2A_Cre_EphB2_sg_VS_2B_Cre_SC_sg.result.csv'
res = read.csv(fn,row.names=1)
png('/media/EphB2/volcano.png')
plot.enhancedVolcano(res, 'mouse', 1, 0.05, '2A_Cre_EphB2_sg_VS_2B_Cre_SC_sg')
dev.off()

plot.enhancedVolcano = function(res, species, lfc_cutoff, alpha, title) {
  ## define y axis
  y_max  = -log10(min(res$padj,na.rm=T)) + 1
  
  ## Add annotation information
  l <- addAnnoByVector(truncateEnsemblID(rownames(res)), species)
  
  EnhancedVolcano(res,
                  lab = l$symbol,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = title,
                  pCutoff = alpha,
                  FCcutoff = lfc_cutoff,
                  pointSize = 1.0,
                  labSize = 3.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75,
                  ylim = c(0, 6))
}
#================================================
#          Nanostring DE analysis
#================================================
cts.mat = read.combined_nano_cts(cts_file)
dim(cts.mat)

## Load sample metadata
sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')

## Verify
if (!verify.meta(cts.mat, sample.meta)) {
  print("Ordering columns of count matrix")
  cts.mat<-cts.mat[,unlist(sample.meta[1])]
  stopifnot(verify.meta(cts.mat.ordered, sample.meta))
}

## Read comparison file
compare_df = read.table(comparison.file,header=T, sep="\t")

DE.DESeq2.NanoString(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, species)

