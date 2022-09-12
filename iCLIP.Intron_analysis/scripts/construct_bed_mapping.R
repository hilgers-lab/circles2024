suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-f", "--feature"), type='character',
              help="Reference database (bed)"),
  
  make_option(c("-g", "--gtf"), type='character',
              help="genome annotation (gtf)"),
  make_option(c("-i", "--introns"), type='character',
              help="Intron database (gtf)"),
  
  make_option(c('-o','--outprefix'), type = 'character',
               help="outprefix for circRNA-intron mapping")
)
# TODO:
#  [open] check if strsplit is necessary here, line 120

opt = parse_args(OptionParser(option_list=option_list))

# opt$ciri.ref="input/circRNA.reference.cutoff2.neuronal.ciri"
# opt$gtf = "resources/dm6_ensembl96.gtf"
# opt$introns="feature_db/dm6_ensembl96.introns.bed"
# opt$circRNA.offset = 25

void = assertthat::assert_that(file.exists(opt$feature), msg = 'Feature file (--feature) does not exist. Exit')
void = assertthat::assert_that(file.exists(opt$gtf), msg = 'genome annotation (--gtf) does not exist. Exit')
void = assertthat::assert_that(file.exists(opt$introns), msg = 'intron annotation (--introns) does not exist. Exit')

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))


# ciri.ref_df <- read_tsv('./resources/circRNA.reference.cutoff2.ciri')
# ciri.ref_df <- read_tsv(opt$ciri.ref)
# circid2jnc_count <- ciri.ref_df %>% dplyr::select(circ_rna_id, number_junction_reads) %>% deframe
# circid2non_jnc_count <- ciri.ref_df %>% dplyr::select(circ_rna_id, number_non_junction_reads) %>% deframe
# ciri.ref_all <- makeGRangesFromDataFrame(ciri.ref_df, keep.extra.columns = TRUE)
# circID2geneid <- melt(ciri_gene_set(ciri.ref_all)) %>%
#   dplyr::select(L1, value) %>% 
#   deframe

features <- import.bed(opt$feature)

dm6.genes <- import.gff(opt$gtf, feature.type = 'gene')
intronic.database <- import.bed(opt$introns)
mcols(intronic.database)$gene_id <- mcols(intronic.database)$name

# Identify flanking introns
circRNA.introns = list()
circRNA.introns[['precede']] = rep(NA, length(features))
circRNA.introns[['follow']] = rep(NA, length(features))


p <- precede(features, intronic.database)
intronic.precede <- intronic.database[na.omit(p)]
names(intronic.precede) = mcols(ciri.ref)$circ_rna_id[!is.na(p)]
mcols(intronic.precede)$circ_rna_id <- gsub('\\|','-',names(intronic.precede))

q <- follow(features, intronic.database)
intronic.follow <- intronic.database[na.omit(q)]
names(intronic.follow) = mcols(ciri.ref)$circ_rna_id[!is.na(q)]
mcols(intronic.follow)$circ_rna_id <- gsub('\\|','-',names(intronic.follow))


featureID.common = intersect(names(intronic.follow), names(intronic.precede))
ddf <- data.frame(feature_id = mcols(features)$name, 
                  feature_id.follow = sapply(mcols(intronic.followfeatureID.common)$gene_id, paste0, collapse = ','), 
                  feature_id.precede = sapply(mcols(intronic.precede[featureID.common])$gene_id,paste0,collapse = ','))

# Find gene overlap
circRNA.gene_mapping = data.frame(feature_id = mcols(feature)$name,
                                  intron_db.downstream = p[mcols(feature)$name %in% featureID.common],
                                  intron_db.upstream = q[mcols(feature)$name %in% featureID.common],
                                  do.call('rbind', 
                                          mapply(function(x, y, z){
                                            c(precede = any(x %in% y), 
                                              follow = any(x %in% z))
                                            }, 
                                            ddf$gene_id, # check if strsplit is necessary here
                                            strsplit(ddf$gene_id.precede, ','), 
                                            strsplit(ddf$gene_id.follow,','), 
                                          SIMPLIFY = FALSE)), row.names = NULL)


circRNA.gene_mapping_clean <- circRNA.gene_mapping %>% filter(!is.na(seqnames) & precede & follow)


outprefix=opt$outprefix
write_tsv(circRNA.gene_mapping, paste0(outprefix,'.intron_mapping.unfiltered.tsv'))
write_tsv(circRNA.gene_mapping_clean, paste0(outprefix,'.intron_mapping.tsv'))
write_tsv(circRNA.gene_mapping_clean %>% dplyr::select(gene_id), paste0(outprefix,'.intron_mapping.gene_set.tsv'))

# export.gff3(gr0,paste0(outprefix,'.igv_only.gff'))
# export.gff3(intronic.database[circRNA.gene_mapping_clean$intron_db.downstream],paste0(outprefix,'.introns_downstream.gff'))
# export.gff3(intronic.database[circRNA.gene_mapping_clean$intron_db.upstream],paste0(outprefix,'.introns_upstream.gff'))
