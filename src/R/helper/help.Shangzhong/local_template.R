source('src/R/util.load.cts.R')
source('src/R/util.load.gtf.R')
source('src/R/util.verify.R')
source('src/R/util.design.R')
source('src/R/util.sampleVariance.R')
source('src/R/util.DE.R')
source('src/R/util.pathway.R')
source('src/R/util.pheatmap.R')
source('src/R/util.expression.boxplot.R')
source('src/R/util.D4C.R')
source('src/R/util.load.metadata.R')
## Load parameter◘
source('src/R/helper/help.Shangzhong/local_param.R')


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
topN = 1e5
DE.DESeq(cts.mat, sample.meta, compare_df, min_total_count, outdir, lfc_cutoff, alpha, topN, species)

## run D4C analysis
dir.create(d4c_out_dir)
run_d4c_for_raw_count(cts.mat, sample.meta, compare_df, domain, user, pwd)
run_d4c_for_deseq2_results(compare_df, domain, user, pwd)




#=============================================
#        pheatmap for subset of genes
#=============================================
gene_fn = '/home/rstudio/p001_RNASeq/CD38/genes_nad+.txt'
outFig = '/home/rstudio/p001_RNASeq/CD38/genes_nad+.png'
norm_cts_fn = '/home/rstudio/p001_RNASeq/CD38/DESeq2/norm_count.tsv'
sample_meta_fn = '/home/rstudio/p001_RNASeq/CD38/20220617.CD38.Mouse.Lung.Metadata_sub.tsv'
variables = c('Condition')
logTrans = T
z_score = F
show_rownames = T
species = 'mouse'
plot_order = c('YoungWT:Saline', 'AgingWT:Saline', 'AgingWT:Bleomycin')
# plot_order = c('Normoxia','IsotypeControl_3mg','Ab1076_0.3mg','Ab1076_3mg',
#                'Ab732_0.3mg','Ab732_3mg','ActRIIA-Fc_2.1mg')
plot.pheatmap(norm_cts_fn, sample_meta_fn, variables, species, gene_fn, outFig, 
              logTrans, show_rownames, z_score, plot_order)

deseq2_fn = '/home/rstudio/p001_RNASeq/CD38/DESeq2/AgingWT_Bleomycin_VS_AgingWT_Saline/AgingWT_Bleomycin_VS_AgingWT_Saline.result.csv'
plot_pheatmap_group_mean(norm_cts_fn, sample_meta_fn, variables, species, gene_fn, deseq2_fn, 
                                    outFig, logTrans,show_rownames, z_score, plot_order)
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
norm_ctx_fn = '/home/rstudio/p001_RNASeq/CD38/DESeq2/norm_count.tsv'
gene_fn = '/home/rstudio/p001_RNASeq/CD38/genes_nad.txt'
sample_meta_fn = '/home/rstudio/p001_RNASeq/CD38/20220617.CD38.Mouse.Lung.Metadata_sub.tsv'
figure = '/home/rstudio/p001_RNASeq/CD38/genes_nad_boxplot.png'
species = 'mouse'
logTrans = T
variables = c('Condition','Gender')
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

