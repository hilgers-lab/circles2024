suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-c", "--ciri.ref"), type='character',
              help="Reference database (ciri2)"),
  
  make_option(c("-g", "--gtf"), type='character',
              help="genome annotation (gtf)"),
  make_option(c("-i", "--introns"), type='character',
              help="Intron database (gtf)"),
  
  make_option(c("-s", "--circRNA.offset"), default = 25, type = 'integer',
              help="Merge circRNA with maximum distance offset between boundaries (default %default)"),
  
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

void = assertthat::assert_that(file.exists(opt$ciri.ref), msg = 'reference circRNA (--ciri.ref) does not exist. Exit')
void = assertthat::assert_that(file.exists(opt$gtf), msg = 'genome annotation (--gtf) does not exist. Exit')
void = assertthat::assert_that(file.exists(opt$introns), msg = 'intron annotation (--introns) does not exist. Exit')

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))

ciri_gene_set <- function(ciri){
  t <- strsplit(mcols(ciri)$gene_id, split = ',')
  names(t) <- mcols(ciri)$circ_rna_id
  t
}

# circ.class <- read_tsv('./data/circRNA.reference.cutoff2_overlaps.tsv')
# subset(circ.class, deltaelav.male_neurons.pooled == 1)

circRNA_offset_thrs <- opt$circRNA.offset 

# ciri.ref_df <- read_tsv('./resources/circRNA.reference.cutoff2.ciri')
ciri.ref_df <- read_tsv(opt$ciri.ref)
circid2jnc_count <- ciri.ref_df %>% dplyr::select(circ_rna_id, number_junction_reads) %>% deframe
circid2non_jnc_count <- ciri.ref_df %>% dplyr::select(circ_rna_id, number_non_junction_reads) %>% deframe
ciri.ref_all <- makeGRangesFromDataFrame(ciri.ref_df, keep.extra.columns = TRUE)
circID2geneid <- melt(ciri_gene_set(ciri.ref_all)) %>%
  dplyr::select(L1, value) %>% 
  deframe


dm6.genes <- import.gff(opt$gtf, feature.type = 'gene')
intronic.database <- import.bed(opt$introns)
mcols(intronic.database)$gene_id <- mcols(intronic.database)$name


# Collapse circRNA where the boundaries almost match
q.self <- findOverlaps(ciri.ref_all, drop.self = TRUE, drop.redundant = TRUE)

d.tss <- distance(resize(ciri.ref_all, width = 1, fix = 'start')[queryHits(q.self)], resize(ciri.ref_all, width = 1, fix = 'start')[subjectHits(q.self)])
d.tes <- distance(resize(ciri.ref_all, width = 1, fix = 'end')[queryHits(q.self)], resize(ciri.ref_all, width = 1, fix = 'end')[subjectHits(q.self)])

dist.df <- data.frame(circ_rna_id = mcols(ciri.ref_all)$circ_rna_id[queryHits(q.self)],
                      coord = gsub('\\|','-',mcols(ciri.ref_all)$circ_rna_id[queryHits(q.self)]), 
                      dist_tss = d.tss,
                      dist_tes = d.tes)

max.offset <- circRNA_offset_thrs
q.to_collapse = q.self[apply(dist.df[,3:4], 1, function(v, c) all(0 <= v & v < c), max.offset),'circ_rna_id']

i.discard = apply(as.data.frame(q.to_collapse),1, function(ii, ref) ii[which.min(c(ref[ii[1]], ref[ii[2]]))], width(ciri.ref_all))

ciri.ref <- ciri.ref_all
if(length(i.discard) > 0)
  ciri.ref <- ciri.ref[-i.discard]

names(ciri.ref) <- mcols(ciri.ref)$circ_rna_id

# Identify flanking introns
circRNA.introns[['precede']] = rep(NA, length(ciri.ref))
circRNA.introns[['follow']] = rep(NA, length(ciri.ref))

p <- precede(ciri.ref, intronic.database)
intronic.precede <- intronic.database[na.omit(p)]
names(intronic.precede) = mcols(ciri.ref)$circ_rna_id[!is.na(p)]
mcols(intronic.precede)$circ_rna_id <- gsub('\\|','-',names(intronic.precede))

q <- follow(ciri.ref, intronic.database)
intronic.follow <- intronic.database[na.omit(q)]
names(intronic.follow) = mcols(ciri.ref)$circ_rna_id[!is.na(q)]
mcols(intronic.follow)$circ_rna_id <- gsub('\\|','-',names(intronic.follow))


circID.common = intersect(names(intronic.follow), names(intronic.precede))
ddf <- data.frame(circ_rna_id = circID.common, 
                  gene_id = circID2geneid[circID.common],
                  gene_id.follow = sapply(mcols(intronic.follow[circID.common])$gene_id, paste0, collapse = ','), 
                  gene_id.precede = sapply(mcols(intronic.precede[circID.common])$gene_id,paste0,collapse = ','))

# Find gene overlap
circRNA.gene_mapping = data.frame(circ_rna_id = ddf$circ_rna_id, 
                                  jnc_count = circid2jnc_count[ddf$circ_rna_id],
                                  non_jnc_count = circid2non_jnc_count[ddf$circ_rna_id],
                                  data.frame(dm6.genes, row.names = mcols(dm6.genes)$gene_id)[ddf$gene_id,c('seqnames','start','end','strand','gene_id','gene_name')],
                                  intron_db.downstream = p[mcols(ciri.ref)$circ_rna_id %in% circID.common],
                                  intron_db.upstream = q[mcols(ciri.ref)$circ_rna_id %in% circID.common],
                                  do.call('rbind', 
                                          mapply(function(x, y, z){
                                            c(precede = any(x %in% y), 
                                              follow = any(x %in% z))
                                            }, 
                                            ddf$gene_id, # check if strsplit is necessary here
                                            strsplit(ddf$gene_id.precede, ','), 
                                            strsplit(ddf$gene_id.follow,','), 
                                          SIMPLIFY = FALSE)), row.names = NULL)

# flip, if intron is too short - has been moved to the upstream function
# b.precede <- width(intronic.database)[circRNA.gene_mapping$intron_db.upstream] > opt$intron.size.threshold
# b.follow <- width(intronic.database)[circRNA.gene_mapping$intron_db.downstream] > opt$intron.size.threshold

# cat('Fraction too short:\n');print(c(upstream = sum(!b.precede), downstream = sum(!b.follow))/nrow(circRNA.gene_mapping))

# circRNA.gene_mapping$precede[!b.precede] = FALSE
# circRNA.gene_mapping$follow[!b.follow] = FALSE

circRNA.gene_mapping_clean <- circRNA.gene_mapping %>% filter(!is.na(seqnames) & precede & follow)

outprefix=opt$outprefix
write_tsv(circRNA.gene_mapping, paste0(outprefix,'.intron_mapping.unfiltered.tsv'))
write_tsv(circRNA.gene_mapping_clean, paste0(outprefix,'.intron_mapping.tsv'))
write_tsv(circRNA.gene_mapping_clean %>% dplyr::select(gene_id), paste0(outprefix,'.intron_mapping.gene_set.tsv'))

# 
# gr0 <- circRNA.gene_mapping_clean %>%
#   dplyr::select(-seqnames,-start,-end,strand) %>% 
#   left_join(ciri.ref_df %>% dplyr::select(circ_rna_id, chr, circ_rna_start, circ_rna_end, strand)) %>% 
#   dplyr::select(chr, circ_rna_start, circ_rna_end, strand, everything()) %>% 
#   filter(!is.na(chr)) %>% 
#   makeGRangesFromDataFrame(keep.extra.columns = TRUE)
#
# export.gff3(gr0,paste0(outprefix,'.igv_only.gff'))
# export.gff3(intronic.database[circRNA.gene_mapping_clean$intron_db.downstream],paste0(outprefix,'.introns_downstream.gff'))
# export.gff3(intronic.database[circRNA.gene_mapping_clean$intron_db.upstream],paste0(outprefix,'.introns_upstream.gff'))
