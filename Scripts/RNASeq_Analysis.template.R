scriptPath = 'X:\\projects\\p038_RNASeq_ESI_Framework\\Scripts'
source(file.path(scriptPath, 'util.load.cts.R'))
source(file.path(scriptPath, 'util.load.gtf.R'))
source(file.path(scriptPath, 'util.verify.R'))
source(file.path(scriptPath, 'util.design.R'))
source(file.path(scriptPath, 'util.PCA.R'))

## Load counts
## Option 1: load Sample counts per file
pattern = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\RNASeq_140.mm10.12_07_2021_10_46_07\\RSEM\\*.rsem.genes.results"
metric.vec = c(1,6) # gene + metric
cts.mat = read_merge.cts(file, metric.vec)
# set sample names
dim(cts.mat)

## Option 2: load count matrix
cts_file = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\RNASeq_140.mm10.12_07_2021_10_46_07\\RSEM\\RSEM.RAW.table"
cts.mat = read.combined_cts(cts_file)


## Load genome annotation file
gtf.file <- 'X:\\projects\\Ensembl\\release-98\\gtf\\Mus_musculus.GRCm38.98.chr.gtf'
gene_anno = load.gtf(gtf.file)

## Output directory
out_dir <- 'X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\output'
setwd(out_dir)

## Load sample metadata
sample.meta.file = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\SampleMeta.txt"
sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE) %>%
  select(Sample, TimePoint, CellLine, Group)%>%
  arrange(Sample)
table(sample.meta$CellLine, sample.meta$TimePoint)

## Verify
stopifnot(verify.meta(cts.mat, sample.meta))

## Variables in design
variables<-c('TimePoint', 'CellLine')
design.formula = build.formula(variables)

## build dge and dds
dge <- DGEList(cts.mat, samples=sample.meta)
design <- model.matrix(design.formula, data=dge$samples)
keep <- filterByExpr(dge, design, min.count=10, min.total.count=50)
dge.filter <- dge[keep, , keep.lib.sizes=F]
dim(dge.filter)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(dge.filter$counts)),
                              colData = dge.filter$samples,
                              design = design.formula)

## Sample Variance Section 
pdf(file='SampleVariance.pdf')
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
#head(assay(vsd),3)
#colData(vsd)
#head(assay(rld), 3)
dds <- estimateSizeFactors(dds) #required for log2 approach which need size factors to account for sequencing depth, and specify normalized=TRUE

plot.mds(dge.filter, variables)
plot.glmpca(dds, variables)
plot.vst(vsd, rld)
plot.sampleDist.vst(vsd, variables, '_')
plot.mds.vst(vsd, variables)
plot.pca.vst(vsd, variables)
plot.svaseq(dge.filter, 'Sample')
dev.off()


## DE
