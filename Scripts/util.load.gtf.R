
load.gtf = function(gtf.file) {
  gtf.df <- as.data.frame(rtracklayer::import(gtf.file))
  gene_annotations <- gtf.df%>%filter(type=="gene")%>%
    select(ensembl_ID=gene_id, gene_type=gene_biotype, symbol=gene_name, chr=seqnames)
}