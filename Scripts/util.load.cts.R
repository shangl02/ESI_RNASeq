
## Count file has metrics for genes per line
## Example: Geneid	Chr	Start	End	Strand	Length	Count
## v: c(1,7) indicates geneID and count
load.cts = function(pattern, v) {
  stopifnot(length(v)==2)
  
  fnames <- Sys.glob(pattern)
  all.counts.list <- lapply(fnames, function(x) read.delim(x, comment.char="#")[, v])
}

verify.cts = function(cts_list) {
  all(sapply(cts_list, function(x) all(x[,1]==cts_list[[1]][, 1])))
}

merge.cts = function(cts_list) {
  all.counts.mat <- do.call("cbind", lapply(cts_list, function(x) x[, 2]))
}

read_merge.cts = function(pattern, v) {
  ## load all counts
  cts.list <- load.cts(pattern, v)
  
  ## verify gene names
  stopifnot(verify.cts(cts.list))
  
  ## merge count matrix
  cts.all.mat <- merge.cts(cts.list)
  
  ## set row names (genes)
  rownames(cts.all.mat) <- cts.list[[1]][, 1]
  return(cts.all.mat)
}



## Combined counts
read.combined_cts = function(file) {
  cts.list = read.table(file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  genes <- cts.list[,1]
  cts.list[1] = NULL
  cts.mat = matrix(unlist(cts.list), nrow=length(genes), ncol = length(cts.list))
  rownames(cts.mat) <- genes
  colnames(cts.mat) <- colnames(cts.list)
  return(cts.mat)
}