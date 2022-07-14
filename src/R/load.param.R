## Parameter
## Count file
cts_file = "X:\\projects\\p041_PMS1\\REDA\\PMS1.allSamples.count.txt"

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

## label column used in PCA
labelCol='Sample'

## variables in sample meta
variables<-c('Sex', 'Treatment')  ## Variables in design

## output dir
outdir <- 'X:\\projects\\p041_PMS1\\Output'

## cutoff for DE
lfc_cutoff=1
alpha=0.05
topN=100


## d4c parameters 
domain = 'CD38'
user = 'mavershang'
pwd = 'TargetScience202206$'
d4c_code_path = '/home/rstudio/Code/ESI_RNASeq/src/python/'
