library("tools")
library("readxl")

load.meta<-function(file) {
 ext=file_ext(file)
 if (ext=='xls' | ext=='xlsx') {
   sample.meta <- read_excel(file, sheet = 1)
 } else {
   sample.meta <- read.csv(sample.meta.file, sep='\t', header=TRUE, stringsAsFactors = FALSE)
 }
}