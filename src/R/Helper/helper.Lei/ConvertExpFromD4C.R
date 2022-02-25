library(stringr)

id="GSE74201"
wd=paste0('X:/projects/p049_PIAS1/counts/', id)
exp.file = paste0('Gene_Expression.', id, '.tbl')

## many exp files from D4C has a data description in the 1st line, set it to 0 if it is not applicable
skip_line=1  

df_exp = read.table(file.path(wd, exp.file), skip = skip_line, row.names=1, header=TRUE, sep='\t')

## Parse gene from column names
genes<-extract.gene(colnames(df_exp))

## Parse sample from row names
samples<-extract.sample(row.names(df_exp))

## 
colnames(df_exp)<-genes
row.names(df_exp)<-samples
df_exp_t<-t(df_exp)

write.table(df_exp_t, file.path(wd, paste0(exp.file, ".parsed.txt")), row.names = TRUE, col.names = TRUE, sep='\t')

## change it accordingly
extract.gene = function(v) {
  ##str_replace(v, '\\S+type\\.\\.gene\\.\\.name\\.\\.(\\w+)\\.\\.label\\.\\.\\S+','\\1')
  str_replace(v, '\\S+name\\.\\.(\\w+)\\.\\.type\\.\\.gene\\.\\.label\\.\\.\\S+','\\1')
  
}

extract.sample = function(v) {
  str_replace(v, '.+\\,\\s*label\\:\\s*GSM\\s*(\\w+)\\}', '\\1')
}
