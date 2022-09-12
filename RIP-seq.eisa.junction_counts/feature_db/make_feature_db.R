suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))

txdb0 <- makeTxDbFromGFF('./resources/dm6_ensembl96.gtf')

# intronic_db = intronsByTranscript(txdb0, use.names = TRUE)
exonic_db = exonsBy(txdb0, by = 'tx', use.names = TRUE)

exonic_db = exonic_db %>% as.data.frame() %>% 
  group_by(group_name) %>%
  filter(exon_rank != 1) %>% 
  makeGRangesListFromDataFrame(names.field = 'group_name') %>% 
  unlist %>% reduce

export.bed(exonic_db, './dm6_ensembl96.exonic_database.bed')
export.bed(genes(txdb0), './dm6_ensembl96.genes.bed')
# export.bed(intronic_db, './dm6_ensembl96.intronic_database.bed')