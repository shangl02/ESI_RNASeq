##
## ES&I, Target Science, Genetics group 
## Function: convert the JSON-like header and row names from D4C expression table
## Date: 02/25/2022
##

suppressPackageStartupMessages({
  require(stringr)
})

## parameter setting (change it based on your file structure)
id="GSE159940"
wd=paste0('X:/projects/p049_PIAS1/counts/', id)
exp.file = paste0('Gene_Expression.', id, '.tbl')
skip_line=1  ## many exp files from D4C has a data description in the 1st line, set it to 0 if it is not applicable

## load expression table
df_exp = read.table(file.path(wd, exp.file), skip = skip_line, row.names=1, header=TRUE, sep='\t')

## Parse gene from column names
## change it accordingly. 
extract.gene = function(v) {
  ##str_replace(v, '\\S+type\\.\\.gene\\.\\.name\\.\\.(\\w+)\\.\\.label\\.\\.\\S+','\\1')
  str_replace(v, '\\S+name\\.\\.(\\w+)\\.\\.type\\.\\.gene\\.\\.\\.label\\.\\.\\S+','\\1')
}
genes<-extract.gene(colnames(df_exp))
head(genes)

## Parse sample from row names
extract.sample = function(v) {
  #str_replace(v, '.+\\,\\s*label\\:\\s*GSM\\s*(\\w+)\\}', '\\1')
  str_replace(v, '.+label\\:\\s*GSM\\s*(\\w+).+', '\\1')
}
samples<-extract.sample(row.names(df_exp))
head(samples)

## Change row/column names and transpose
colnames(df_exp)<-genes
row.names(df_exp)<-samples
df_exp_t<-t(df_exp)

## write out
write.table(df_exp_t, file.path(wd, paste0(exp.file, ".parsed.txt")), row.names = TRUE, col.names = TRUE, sep='\t')


