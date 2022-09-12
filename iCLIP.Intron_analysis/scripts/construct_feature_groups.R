suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-m", "--feature_map"), type='character', default = NULL,
              help="Create up/downstream introns by intron map of circRNA to intron mapping (tsv)"),
  
  make_option(c("-i", "--introns"), type='character',
              help="Intron database (gtf)"),
  
  make_option(c('-o','--outprefix'), type='character',
              help="outprefix for intron sets")
)

opt = parse_args(OptionParser(option_list=option_list))
# opt$introns = 'feature_db/dm6_ensembl96.introns.bed'
# opt$feature_map = 'Features/ciri_ref/circRNA.intron_mapping.tsv'

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

intron.db <- import.bed(opt$introns)
feature.map <- read_tsv(opt$feature_map)

groups.set = list()
t = intron.db[feature.map$upstream.idx]
mcols(t)$circ_rna_id = mcols(t)$feature_id = feature.map$circ_rna_id
groups.set[['upstream']] = t
t = intron.db[feature.map$downstream.idx]
mcols(t)$circ_rna_id = mcols(t)$feature_id = feature.map$circ_rna_id
groups.set[['downstream']] = t
t = subset(intron.db[-c(feature.map$downstream.idx,feature.map$upstream.idx),], name %in% feature.map$gene_id)
mcols(t)$feature_id = data.frame(gene_id = mcols(t)$name) %>%
  group_by(gene_id) %>% 
  mutate(feature_id = paste0(gene_id,'.intron',1:n())) %>% 
  ungroup() %>% 
  select(feature_id) %>% unlist %>% as.vector()
groups.set[['notflanking']] = t

for(n0 in names(groups.set)){
  export.gff3(groups.set[[n0]], paste0(opt$outprefix,'_', n0,'.gff'))
}


