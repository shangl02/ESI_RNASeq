source('src/R/util.load.cts.R')
source('src/R/util.load.metadata.R')
source('src/R/util.load.gtf.R')
source('src/R/util.verify.R')
source('src/R/util.design.R')
source('src/R/util.sampleVariance.R')
source('src/R/util.DE.R')
source('src/R/util.pathway.R')
source('src/R/util.DE.Nanostring.R')
source('src/R/util.pheatmap.R')


## Load parameter◘
#source('src/R/load.param.R')

## Load counts
## Option 1: load Sample counts per file
cts.mat = read_merge.cts(pattern, metric.vec, patName)
# set sample names

## Option 2: load count matrix
cts.mat = read.combined_cts(cts_file)

dim(cts.mat)

## Load genome annotation file
#gtf.file <- 'X:\\projects\\Ensembl\\release-98\\gtf\\Mus_musculus.GRCm38.98.chr.gtf'
#gene_anno = load.gtf(gtf.file)


## Load sample metadata
sample.meta <- load.meta(sample.meta.file)
stopifnot("Sample" %in% colnames(sample.meta))
sample.meta$mergeCond = pasteMultiCol(sample.meta, variables, '.')
dim(sample.meta)
# 
# ## add filter of based on sample metadata
# sample.meta=sample.meta[sample.meta$Quality=='good',]
# cts.mat = cts.mat[,colnames(cts.mat) %in% sample.meta$Sample]
# dim(cts.mat)
# dim(sample.meta)

## Verify count table and sample metadata
if (!verify.meta(cts.mat, sample.meta)) {
  print("Ordering columns of count matrix")
  cts.mat<-cts.mat[,unlist(sample.meta[1])]
  stopifnot(verify.meta(cts.mat, sample.meta))
}

## Read comparison file
compare_df = read.table(comparison.file,header=T, sep="\t")

## Output directory
dir.create(outdir)
setwd(outdir)

# Generate sample variance report
process.sampleVariance.all(cts.mat, sample.meta, variables, min_count, min_total_count, labelCol)

## DE analysis
if (exists("DE_algo") && DE_algo=="limma") {
  DE.limma(cts.mat, sample.meta, compare_df, min_count, min_total_count, outdir, lfc_cutoff, alpha, topN, species)
} else {
  DE.DESeq(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, topN, species)
}

## run GSEA
species = 'mouse'
category = 'C2' # chose GSE pathway type
deseq2_path = outdir
GSEA_wrapper(compare_df, sample.meta, deseq2_path, species)

## run D4C analysis
dir.create(d4c_out_dir)
run_d4c_for_raw_count(cts.mat, sample.meta, compare_df, domain, user, pwd)
run_d4c_for_deseq2_results(compare_df, domain, user, pwd)