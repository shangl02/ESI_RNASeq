## Parameter
## Count file
cts_file = "/media/C3aR/raw_counts.tsv"


## Load counts
pattern = "X:\\projects\\p041_PMS1\\AlignmentNCounting\\S*\\salmon_aln\\quant.genes.sf"
metric.vec = c(1,5) # gene + count or TMP (depend on the counter)
patName = 'S\\d+'

## sample meta file
sample.meta.file = "/media/C3aR/condition.tsv"
comparison.file = "/media/C3aR/compare.tsv"

## Species
species='human'

## gene filter
min_count=5
min_total_count=30

## variables in sample meta
variables<-c('condition')  ## Variables in design

## output dir
outdir <- '/media/C3aR//DESeq2'

## cutoff for DE
lfc_cutoff=1
alpha=0.05

