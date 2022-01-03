build.formula = function(v) {
  as.formula(paste0("~ ", paste(v, collapse=" + ")))
}