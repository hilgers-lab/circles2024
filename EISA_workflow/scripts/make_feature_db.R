suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))

# Arguments: <input_gtf>  <outdir>
args0 = commandArgs(TRUE)
# txdb0 <- makeTxDbFromGFF('./resources/dm6_ensembl96.gtf')
txdb0 <- makeTxDbFromGFF(args0[1])

# intronic_db = intronsByTranscript(txdb0, use.names = TRUE)
exonic_db = exonsBy(txdb0, by = 'tx', use.names = TRUE)

exonic_db = exonic_db %>% as.data.frame() %>% 
  group_by(group_name) %>%
  filter(exon_rank != 1) %>% 
  makeGRangesListFromDataFrame(names.field = 'group_name') %>% 
  unlist %>% reduce

outprefix=paste0(args0[2], '/', gsub('(.*).gtf','\\1', basename(args0[1])))
export.bed(exonic_db, paste0(outprefix,'.exonic_database.bed'))
export.bed(genes(txdb0), paste0(outprefix,'.genes.bed'))

# export.bed(intronic_db, './dm6_ensembl96.intronic_database.bed')
