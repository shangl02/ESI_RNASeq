
verify.meta = function(cts.mat, sample.meta) {
  return(all(colnames(cts.mat)==sample.meta[,1]))
}
