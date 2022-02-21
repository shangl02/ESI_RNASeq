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
source('src/R/local_param.R')

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
  rld <- rlog(dds, blind = FALSE)
  dds <- estimateSizeFactors(dds) #required for log2 approach which need size factors to account for sequencing depth, and specify normalized=TRUE
  l <- calc.sampleDist.vst(vsd, variables, '_') ## calculate sample distance matrix
  
  ## plotting
  tryCatch(
    {
      print("Generating pdf")
      pdf(file='SampleVariance.all.pdf')
      # p1<-plot.mds(dge.filter, variables)
      # p2<-plot.glmpca(dds, variables)
      # p3<-plot.vst(dds,vsd, rld)
      # p4<-plot.sampleDist.vst(l$dist, l$dist.mat)
      # p5<-plot.mds.vst(vsd, variables, l$dist.mat)
      p6<-plot.pca.vst(vsd, variables)
      # p7<-plot.svaseq(dge.filter, 'Sample')
      ## plot.topNVar.vst(vsd, 20)
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


#====== pheatmap for subset of genes========
gene_fn = '/media/C3aR/pheatmap//genes.txt'
outFig = '/media/C3aR/pheatmap//heatmap.png'
logTrans = F
plot.pheatmap(cts_file, sample.meta.file, variables, gene_fn, outFig, logTrans=F)

#==== merge all DESeq results ================
fn1 = '/media/CD38/DESeq2/WT_CDA-HFD_VS_WT_Regular/WT_CDA-HFD_VS_WT_Regular.result.csv'
fn2 = '/media/CD38/DESeq2/CD38KO_CDA-HFD_VS_CD38KO_Regular/CD38KO_CDA-HFD_VS_CD38KO_Regular.result.csv'
fn3 = '/media/CD38/DESeq2/CD38KO_CDA-HFD_VS_WT_CDA-HFD/CD38KO_CDA-HFD_VS_WT_CDA-HFD.result.csv'
df1 = read.csv(fn1)
df2 = read.csv(fn2)
df3 = read.csv(fn3)
res =  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "X", all.x = TRUE),
              list(df1,df2,df3))
write.table(res,'/media/CD38/DESeq2/merge_res.tsv',sep='\t',quote=F)
df = read.table('/media/CD38/DESeq2/merge_lfc_qval.tsv',header=T)
sub_df = df[df$symbol %in% genes,]
sub_df = sub_df[order(sub_df$symbol),]
write.table(sub_df, '/media/CD38/DESeq2/sub_lfc_qval.tsv',sep='\t', 
            quote=F,row.names=F)

#====== gene expression boxplot =================
norm_ctx_fn = '/media/IL26/tpm.tsv'
gene_fn = '/media/IL26/genes.txt'
sample_meta_fn = '/media/IL26/condition.tsv'
figure = '/media/IL26/genes.png'
species = 'human'
expression_boxplot(norm_ctx_fn, gene_fn, sample_meta_fn, variables, species, figure)

#====== Nanostring DE analysis ==================
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
