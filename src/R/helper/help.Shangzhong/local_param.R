## Parameter
## Count file
cts_file = "/media/CD38/raw_count.tsv"


## Load counts
pattern = "X:\\projects\\p041_PMS1\\AlignmentNCounting\\S*\\salmon_aln\\quant.genes.sf"
metric.vec = c(1,5) # gene + count or TMP (depend on the counter)
patName = 'S\\d+'

## sample meta file
sample.meta.file = "/media/CD38/condition.tsv"
comparison.file = "/media/CD38/compare.tsv"

## Species
species='mouse'

## gene filter
min_count=5
min_total_count=30

## variables in sample meta
# variables<-c('Cytokine','Nucleotide','Treatment')  ## Variables in design
variables = c('Strain','Diet')

## output dir
outdir <- '/media/CD38/DESeq2'

## cutoff for DE
lfc_cutoff=0.58
alpha=0.05

