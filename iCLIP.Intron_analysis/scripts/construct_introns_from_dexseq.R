suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-c", "--dexseq.targets"), type='character',
              help="circRNA targets (ciri2)"),

  make_option(c("-g", "--genes"), type='character',
              help="Prepared gene database (bed)"),
  make_option(c("-i", "--introns"), type='character',
              help="Prepared Intron database (bed)"),

  # make_option(c( "--collapse_circRNA"), default = 10, type = 'integer', # apply 10 nt
  #             help="Collapse overlapping circRNA, if boundary distances are below <int> (default %default)"),
  make_option(c( "--exon_boundaries"), default = 10, type = 'integer', # apply 10 nt
              help="Discard exons with distance to exon boundaries greater than <int> (default %default)"),
  make_option(c("--intron_maxgap"), default = 10, type = 'integer', # apply 10 nt
              help="Assigned introns with max distance <int> (default %default)"),

  make_option(c('-o','--outprefix'), type = 'character',
              help="outprefix for circRNA-intron mapping")
)
# defined output:
# - intron index for feature set
# - intron gff (feature sorted) for upstream/downstream/non-flanking (mind the group of circRNA!)

opt = parse_args(OptionParser(option_list=option_list))

void = assertthat::assert_that(file.exists(opt$dexseq.targets), msg = 'reference dexseq exons (--dexseq.targets) does not exist. Exit')
void = assertthat::assert_that(file.exists(opt$genes) & file.exists(opt$introns), msg = 'Gene, exon or intron database is missing. Exit')

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(GenomicFeatures))

# 1) Select circRNA matching exon boundaries, max offset

# input: ciri2, intron database
# return: ciri2-intron-map <upstream>/<downstream>
create_intron_map <- function(exons, introns, maxgap = 0){
  cat('>>> create_intron_map\n')
  q.upstream <- findOverlaps(resize(exons,1,'start'), resize(introns,1,'end'), maxgap = maxgap)
  q.downstream <- findOverlaps(resize(exons,1,'end'), resize(introns,1,'start'), maxgap = maxgap)

  void <- assertthat::assert_that(all(elementNROWS(split(subjectHits(q.upstream), queryHits(q.upstream) )) %>% range == 1),
                                  msg = '[debug] found more than one intron per circRNA')
  void <- assertthat::assert_that(all(elementNROWS(split(subjectHits(q.downstream), queryHits(q.downstream) )) %>% range == 1),
                                  msg = '[debug] found more than one intron per circRNA')


  ## Perform gene_id mapping
  intron.map <- data.frame(group_id = mcols(exons)$group_id,
                           feature_id = mcols(exons)$feature_id,
                           upstream.idx = NA,
                           downstream.idx = NA)
  intron.map$upstream.idx[queryHits(q.upstream)] = subjectHits(q.upstream)
  intron.map$downstream.idx[queryHits(q.downstream)] = subjectHits(q.downstream)

  # Find gene overlap
  gene_id.map <- data.frame(group_id = mcols(exons)$group_id,
                            feature_id = mcols(exons)$feature_id,
                            gene_id = mcols(exons)$group_id,
                            upstream = NA,
                            downstream = NA)

  gene_id.map$upstream = mcols(introns)$gene_id[intron.map$upstream]
  gene_id.map$downstream = mcols(introns)$gene_id[intron.map$downstream]
  # match by set of genes
  gene_id.map$upstream.match = mapply(function(x,y) any(x %in% y), strsplit(gene_id.map$gene_id,'\\+'),strsplit(gene_id.map$upstream,'\\+'))
  gene_id.map$downstream.match = mapply(function(x,y) any(x %in% y), strsplit(gene_id.map$gene_id,'\\+'),strsplit(gene_id.map$downstream,'\\+'))

  # compose
  exon.intron_map = data.frame(exons) %>%
    left_join(intron.map, by = c('group_id','feature_id')) %>%
    left_join(gene_id.map %>%
                dplyr::select(-upstream, -downstream), by = c('group_id','feature_id'))
  exon.intron_map
}

get_flanking_introns <- function(exons, introns, intron.map){
  cat('>>> get_flanking_introns\n')
  b.is_single_flanking <- is.na(intron.map$upstream.idx) | is.na(intron.map$downstream.idx)
  warning('> exons with single intron are discarded. Discarding #exons: ', sum(b.is_single_flanking))
  intron.map <- subset(intron.map, !b.is_single_flanking)

  groups.set = list()
  t = introns[intron.map$upstream.idx]
  mcols(t)$group_id =  intron.map$group_id
  mcols(t)$feature_id = intron.map$feature_id
  groups.set[['upstream']] = t

  t = introns[intron.map$downstream.idx]
  mcols(t)$group_id =  intron.map$group_id
  mcols(t)$feature_id = intron.map$feature_id
  groups.set[['downstream']] = t


  t = subset(introns[-c(intron.map$downstream.idx,intron.map$upstream.idx),], name %in% unlist(strsplit(intron.map$group_id,'\\+')))
  mcols(t)$feature_id = data.frame(gene_id = mcols(t)$name) %>%
    group_by(gene_id) %>%
    mutate(feature_name = paste0(gene_id,'.intron',1:n()),
           feature_id = 1:n()) %>%
    ungroup() %>%
    dplyr::select(feature_id) %>% unlist %>% as.vector()
  t <- t[!t %over% exons]
  groups.set[['notflanking']] = t

  list(group_id = intron.map$group_id,
       feature_id = intron.map$feature_id,
       intron.set = groups.set,
       exons_single = exons[b.is_single_flanking])
}

# MAIN
# import exons
exon.targets <- read_tsv(opt$dexseq.targets) %>% janitor::clean_names()
exon.ref_all <- makeGRangesFromDataFrame(exon.targets, keep.extra.columns = TRUE)
mcols(exon.ref_all)$gene_id = gsub('\\+',',',mcols(exon.ref_all)$group_id)
mcols(exon.ref_all)$feature_id = paste0(gsub('\\+',',',mcols(exon.ref_all)$group_id), ":", mcols(exon.ref_all)$feature_id)
# Import data bases
## genes
organism.genes <- import.bed(opt$genes)
mcols(organism.genes)$gene_id <- mcols(organism.genes)$name
organism.introns <- import.bed(opt$introns)
mcols(organism.introns)$gene_id <- mcols(organism.introns)$name

exon.ref <- exon.ref_all

cat('intron maxgap:', opt$intron_maxgap,'\n')
intron.index <- create_intron_map(exon.ref, organism.introns, as.integer(opt$intron_maxgap))

introns.per_exon = get_flanking_introns(exon.ref, organism.introns, intron.index)

intron.set <- introns.per_exon$intron.set

exon.set <- subset(exon.ref,paste0(group_id, feature_id) %in% paste0(introns.per_exon$group_id, introns.per_exon$feature_id))

# select corresponding exons
for(n0 in names(intron.set)){
  export.gff3(intron.set[[n0]], paste0(opt$outprefix,'_dexseq_', n0,'.gff'))
}

export.gff3(exon.set, paste0(opt$outprefix, '_dexseq.gff'))
export.gff3(introns.per_exon$exons_single, paste0(opt$outprefix, '_dexseq_discarded.gff'))
