## Parameter
## Count file
cts_file = "X:\\projects\\p029_SAMHD1\\RNASeq.Mouse.12032021\\AlignmentNCounting\\ORD_Pipeline\\RNASeq_140.mm10.12_07_2021_10_46_07\\RSEM\\RSEM.RAW.table"


## Load counts
pattern = "X:\\projects\\p041_PMS1\\AlignmentNCounting\\S*\\salmon_aln\\quant.genes.sf"
metric.vec = c(1,5) # gene + count or TMP (depend on the counter)
patName = 'S\\d+'

## sample meta file
sample.meta.file = "X:\\projects\\p041_PMS1\\Metadata\\PMS1.SampleMeta.txt"
comparison.file = "X:\\projects\\p041_PMS1\\Metadata\\Comparisons.txt"

## Species
species='mouse'

## gene filter
min_count=5
min_total_count=30

## variables in sample meta
variables<-c('Sex', 'Treatment')  ## Variables in design

## output dir
outdir <- 'X:\\projects\\p041_PMS1\\Output'

## cutoff for DE
lfc_cutoff=1
alpha=0.05

