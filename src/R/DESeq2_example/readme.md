* raw_count.tsv: table that have raw count data of RNA-Seq.
* condition.tsv: table defining the condition for each sample.
* DE_pair.tsv: two columns defining control and test, each line represents one DE comparison.

Command to run the DE analysis
------------------------------

**Rscript util.DESeq2.R raw_count.tsv condition.tsv DE_pair.tsv outPath**