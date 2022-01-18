source('src/R/util.load.cts.R')
source('src/R/util.load.gtf.R')
source('src/R/util.verify.R')
source('src/R/util.design.R')
source('src/R/util.sampleVariance.R')
source('src/R/util.DE.R')
source('src/R/util.pathway.R')

## Parameter
cts_file = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\AlignmentNCounting\\ORD_Pipeline\\RNASeq_140.mm10.12_07_2021_10_46_07\\RSEM\\RSEM.RAW.table"
sample.meta.file = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\SampleMeta.txt"
comparison.file = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\Comparisons.txt"
min_count=5
min_total_count=30
variables<-c('TimePoint', 'CellLine')  ## Variables in design
outdir <- 'X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\output2'
lfc_cutoff=2
alpha=0.05
species='mouse'

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
      p1<-plot.mds(dge.filter, variables)
      p2<-plot.glmpca(dds, variables)
      p3<-plot.vst(dds,vsd, rld)
      p4<-plot.sampleDist.vst(l$dist, l$dist.mat)
      p5<-plot.mds.vst(vsd, variables, l$dist.mat)
      p6<-plot.pca.vst(vsd, variables)
      # p7<-plot.svaseq(dge.filter, 'Sample')
      ## plot.topNVar.vst(vsd, 20)
      print(p1)
      print(p2)
      print(p3)
      print(p4)
      print(p5)
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
pattern = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\RNASeq_140.mm10.12_07_2021_10_46_07\\RSEM\\*.rsem.genes.results"
metric.vec = c(1,6) # gene + metric
cts.mat = read_merge.cts(file, metric.vec)
# set sample names
dim(cts.mat)

## Option 2: load count matrix
cts.mat = read.combined_cts(cts_file)


## Load genome annotation file
gtf.file <- 'X:\\projects\\Ensembl\\release-98\\gtf\\Mus_musculus.GRCm38.98.chr.gtf'
gene_anno = load.gtf(gtf.file)

## Output directory
dir.create(outdir)
setwd(outdir)

## Load sample metadata
sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, ':')

# sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE) %>%
#    select(Sample, TimePoint, CellLine, Group)%>%
#    arrange(Sample)
#table(sample.meta$CellLine, sample.meta$TimePoint)

## Verify
stopifnot(verify.meta(cts.mat, sample.meta))

# Generate sample variance report
process.sampleVariance.all(cts.mat, sample.meta, variables)

## DE analysis
compare_df = read.table(comparison.file,header=T, sep="\t")
DE.DESeq(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, species)
