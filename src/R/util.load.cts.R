
## Count file has metrics for genes per line
## Example: Geneid	Chr	Start	End	Strand	Length	Count
## v: c(1,7) indicates geneID and count
load.cts = function(pattern, v, patName) {
  stopifnot(length(v)==2)
  
  fnames <- Sys.glob(pattern)
  
  all.counts.list <- lapply(fnames, function(x) read.delim(x, comment.char="#")[, v])
  all.names <- lapply(fnames, function(x) str_extract_all(x,patName)[[1]][2])
  return(list(all.counts.list, all.names))
}

verify.cts = function(cts_list) {
  all(sapply(cts_list, function(x) all(x[,1]==cts_list[[1]][, 1])))
}

merge.cts = function(cts_list) {
  all.counts.mat <- do.call("cbind", lapply(cts_list, function(x) x[, 2]))
}

read_merge.cts = function(pattern, v, patName) {
  ## load all counts
  l <- load.cts(pattern, v, patName)
  cts.list <- l[[1]]
  samples <- l[[2]]

  ## verify gene names
  stopifnot(verify.cts(cts.list))
  
  ## merge count matrix
  cts.all.mat <- merge.cts(cts.list)
  
  ## set row names (genes)
  rownames(cts.all.mat) <- cts.list[[1]][, 1]
  
  ## set column names (samples)
  colnames(cts.all.mat) = samples
  
  return(cts.all.mat)
}



## Combined counts
read.combined_cts = function(file) {
  delim=ifelse(tolower(tools::file_ext(file))=='csv', ',','\t')
  cts.list = read.table(file, sep=delim, header=TRUE, stringsAsFactors = FALSE)
  genes <- cts.list[,1]
  cts.list[1] = NULL
  cts.mat = matrix(unlist(cts.list), nrow=length(genes), ncol = length(cts.list))
  rownames(cts.mat) <- genes
  colnames(cts.mat) <- colnames(cts.list)
  return(cts.mat)
}


## Combined Nanostring counts
read.combined_nano_cts = function(file) {
  cts.list = read.table(file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
  genes = paste(cts.list$Code.Class, cts.list$Name, cts.list$Accession, sep='_')
  cts.mat = cts.list[,4:length(cts.list)]
  rownames(cts.mat) = genes
  colnames(cts.mat) = colnames(cts.list)[4:length(cts.list)]
  return(cts.mat)
}
