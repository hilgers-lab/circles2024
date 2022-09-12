
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-g", "--gene_set"), type='character', default = NULL,
              help="Features for gene set"),
  make_option(c("-f", "--format"), type = 'character', default = 'tsv',
               help = 'Format of gene_set. Either \'gff\' or \'tsv\' (default)'),

  make_option(c("-i", "--introns"), type='character',
              help="Intron database (bed)"),

  make_option(c('-o','--outprefix'), type='character',
              help="directory for intron sets")
)

opt = parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

void = assertthat::assert_that(opt$format %in% c('tsv','gff'), msg = 'File format is not supported. Either \'gff\' or \'tsv\'')

intron.db <- import.bed(opt$introns)

if(opt$format == 'tsv'){
  gene_set <- read_tsv(opt$gene_set, show_col_types = FALSE)
} else if (opt$format == 'gff'){
  gene_set <- import.gff(opt$gene_set)
}
groups.set = list()
introns.gene_set = intron.db[sapply(strsplit(mcols(intron.db)$name,','), function(x, y) any(x%in%y), unique(strsplit(gene_set$gene_id,',')))]

mcols(introns.gene_set)$feature_id = data.frame(gene_id = mcols(introns.gene_set)$name) %>%
  group_by(gene_id) %>%
  mutate(feature_id = paste0(gene_id,'.intron',1:n())) %>%
  ungroup() %>%
  select(feature_id) %>% unlist %>% as.vector()
groups.set[['geneset']] = introns.gene_set

for(n0 in names(groups.set)){
  export.gff3(groups.set[[n0]], paste0(opt$outprefix, '_', n0,'.gff'))
}
