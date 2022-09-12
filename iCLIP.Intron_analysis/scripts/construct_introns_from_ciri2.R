suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-c", "--circRNA.targets"), type='character',
              help="circRNA targets (ciri2)"),

  make_option(c("-g", "--genes"), type='character',
              help="Prepared gene database (bed)"),
  make_option(c("-e", "--exons"), type='character',
              help="Prepared exon database (bed)"),
  make_option(c("-i", "--introns"), type='character',
              help="Prepared Intron database (bed)"),

  make_option(c( "--internal_introns"), default = FALSE, type = 'logical', action = 'store_true', 
              help="extract internal introns (default: %default)"),
  make_option(c( "--collapse_circRNA"), default = 10, type = 'integer', # apply 10 nt
              help="Collapse overlapping circRNA, if boundary distances are below <int> (default %default)"),
  make_option(c( "--exon_boundaries"), default = 10, type = 'integer', # apply 10 nt
              help="Discard circRNA with distance to exon boundaries greater than <int> (default %default)"),
  make_option(c("--intron_maxgap"), default = 10, type = 'integer', # apply 10 nt
              help="Assigned introns with max distance <int> (default %default)"),

  make_option(c('-o','--outprefix'), type = 'character',
               help="outprefix for circRNA-intron mapping")
)
# defined output:
# - intron index for feature set
# - intron gff (feature sorted) for upstream/downstream/non-flanking (mind the group of circRNA!)

opt = parse_args(OptionParser(option_list=option_list))

# TODO:
# opt$circRNA.targets = '../iCLIP.introns_analysis.elav_spliced/input/circRNA.neuronal_targets.ciri'
# opt$genes = '../iCLIP.introns_analysis.elav_spliced//feature_db/dm6_ensembl96.genes.bed'
# opt$exons = '../iCLIP.introns_analysis.elav_spliced//feature_db/dm6_ensembl96.exons.bed'
# opt$introns = '../iCLIP.introns_analysis.elav_spliced//feature_db/dm6_ensembl96.introns.bed'
# opt$internal=TRUE

void = assertthat::assert_that(file.exists(opt$circRNA.targets), msg = 'reference circRNA (--circRNA.targets) does not exist. Exit')
void = assertthat::assert_that(file.exists(opt$genes) & file.exists(opt$exons) & file.exists(opt$introns), msg = 'Gene, exon or intron database is missing. Exit')

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(GenomicFeatures))

ciri_gene_set <- function(ciri){
  cat('>>> ciri_gene_set\n')
  t <- strsplit(mcols(ciri)$gene_id, split = ',')
  names(t) <- mcols(ciri)$circ_rna_id
  t
}

# 1) Select circRNA matching exon boundaries, max offset
make_circ_rna_id <- function(ciri.gr){
  cat('>>> make_circ_rna_id\n')
  paste0(seqnames(ciri.gr),':',start(ciri.gr),'|',end(ciri.gr))
}

is_at_exon_boundary <- function(ciri2, exonset, maxgap){
  cat('>>> is_at_exon_boundary\n')
  b.start <- overlapsAny(ciri2, exonset, type='start',maxgap = maxgap)
  b.end <- overlapsAny(ciri2, exonset, type='end', maxgap = maxgap)

  b.start & b.end
}

# <check> boundary
# i.start = unique(queryHits(subset(distanceToNearest(resize(ciri.ref_all,1,'start'), resize(organism.exons,1,'start')), distance > exon_boundary_window)))
# i.end = unique(queryHits(subset(distanceToNearest(resize(ciri.ref_all,1,'end'), resize(organism.exons,1,'end')), distance > exon_boundary_window)))
#
# ciri.discard = NULL
# if(length(union(i.start, i.end)) > 0){
#   ciri.discard = ciri.ref_all[union(i.start,i.end)];ciri.discard
# }
# void <- assertthat::assert_that(length(ciri.ref) + length(ciri.discard) == length(ciri.ref_all), msg = '[debug] circRNA count doesnt add up. Check code')
# /<check>

# 2) Pool self-overlaps with max <exon_boundary_window>
pool_ambiguous_ciri2 <- function(ciri2, maxgap){
  cat('>>> pool_ambiguous_ciri2\n')
  qq <- findOverlaps(ciri2, drop.self = TRUE, drop.redundant = TRUE, type = 'equal', maxgap = maxgap)
  if(length(qq) > 0){
    ciri.collapsed <- reduce(split(ciri2[queryHits(qq)], subjectHits(qq))) %>% as.data.frame() %>% makeGRangesFromDataFrame()
    ciri.set <- sort(c(granges(ciri.ref)[-unlist(as.data.frame(qq))], ciri.collapsed), ignore.strand = TRUE)
    mcols(ciri.set)$circ_rna_id <- make_circ_rna_id(ciri.set)
  } else {
    ciri.set = granges(ciri.ref)
    mcols(ciri.set)$circ_rna_id <- make_circ_rna_id(ciri.set)
  }
  return(ciri.set)
}

# input: ciri2, intron database
# return: ciri2-intron-map <upstream>/<downstream>
create_intron_map <- function(ciri2, introns, internal=FALSE, maxgap = 0){
  cat('>>> create_intron_map\n')
  cat('>>> Internal:',internal,'\n')
  if(internal) {
    # 1) find introns per circRNA
    # 2) select first and last intron
    qi <- findOverlaps(introns, ciri2, type = 'within')
    
    # sort.GRangesList is strand independent
    qi.grouped <- split(queryHits(qi), subjectHits(qi))
    print(length(qi.grouped))
    
    qi.mat <- lapply(qi.grouped, function(ii, ref){
      if(length(ii) == 1) 
        return(c(upstream = ii, downstream = ii))
      else {
        gr = ref[ii]
        assertthat::assert_that(length(unique(strand(gr))) == 1, msg = 'internal introns: Found introns with contradictory strands')
        
        if(all(strand(gr) == "-")) 
          return(c(upstream = ii[length(ii)], downstream = ii[1]))
        else 
          return(c(upstream = ii[1], downstream = ii[length(ii)]))
      }
    }, introns) %>% do.call(what = 'rbind') %>% as.data.frame
    
    intron.map = data.frame(circ_rna_id = mcols(ciri2)$circ_rna_id,
                            upstream.idx = NA,
                            downstream.idx = NA)
    
    intron.map$upstream.idx[as.integer(names(qi.grouped))] = qi.mat$upstream
    intron.map$downstream.idx[as.integer(names(qi.grouped))] = qi.mat$downstream
    
  } else {
    q.upstream <- as.data.frame(findOverlaps(resize(ciri2,1,'start'), resize(introns,1,'end'), maxgap = maxgap))
    q.downstream <- as.data.frame(findOverlaps(resize(ciri2,1,'end'), resize(introns,1,'start'), maxgap = maxgap))
    
    void <- assertthat::assert_that(all(elementNROWS(split(q.upstream$subjectHits, q.upstream$queryHits )) %>% max == 1), msg = '[debug] found more than one intron per circRNA')
    void <- assertthat::assert_that(all(elementNROWS(split(q.downstream$subjectHits, q.downstream$queryHits )) %>% max == 1), msg = '[debug] found more than one intron per circRNA')
    
    
    ## Perform gene_id mapping
    intron.map <- data.frame(circ_rna_id = mcols(ciri2)$circ_rna_id,
                             upstream.idx = NA,
                             downstream.idx = NA)
    intron.map$upstream.idx[q.upstream$queryHits] = q.upstream$subjectHits
    intron.map$downstream.idx[q.downstream$queryHits] = q.downstream$subjectHits
  }
  
  # Find gene overlap
  gene_id.map <- data.frame(circ_rna_id = mcols(ciri2)$circ_rna_id,
                            gene_id = mcols(ciri2)$gene_id,
                            upstream = NA,
                            downstream = NA)
  
  gene_id.map$upstream = mcols(introns)$gene_id[intron.map$upstream]
  gene_id.map$downstream = mcols(introns)$gene_id[intron.map$downstream]
  # match by set of genes
  gene_id.map$upstream.match = mapply(function(x,y) any(x %in% y), strsplit(gene_id.map$gene_id,','),strsplit(gene_id.map$upstream,','))
  gene_id.map$downstream.match = mapply(function(x,y) any(x %in% y), strsplit(gene_id.map$gene_id,','),strsplit(gene_id.map$downstream,','))
  
  # compose
  circRNA.intron_map = data.frame(circ_rna_id = mcols(ciri2)$circ_rna_id,ciri2) %>%
    left_join(intron.map, by = 'circ_rna_id') %>%
    left_join(gene_id.map %>%
                dplyr::select(-upstream, -downstream), by = c('circ_rna_id','gene_id'))
  circRNA.intron_map
}


make_intron_set <- function(ciri2, introns, intron.map, is.internal = FALSE){
  
  if(is.internal){
    
    cat('Single exon circRNA:', sum(is.na(intron.map$upstream.idx)),'\n')
    intron.map <- subset(intron.map, !is.na(upstream.idx))
    
    b.single = intron.map$upstream.idx == intron.map$downstream.idx
    
    # split single and multi introns
    groups.set = list()
    t = introns[intron.map[b.single,]$upstream.idx]
    mcols(t)$circ_rna_id = mcols(t)$feature_id = intron.map$circ_rna_id[b.single]
    groups.set[['internal_single']] = t
    
    t = introns[intron.map[!b.single,]$upstream.idx]
    mcols(t)$circ_rna_id = mcols(t)$feature_id = intron.map$circ_rna_id[!b.single]
    groups.set[['internal_upstream']] = t
    
    t = introns[intron.map[!b.single,]$downstream.idx]
    mcols(t)$circ_rna_id = mcols(t)$feature_id = intron.map$circ_rna_id[!b.single]
    groups.set[['internal_downstream']] = t
    
    t = subset(introns[-c(intron.map$downstream.idx,intron.map$upstream.idx),], name %in% intron.map$gene_id)
    mcols(t)$feature_id = data.frame(gene_id = mcols(t)$name) %>%
      group_by(gene_id) %>%
      mutate(feature_name = paste0(gene_id,'.intron',1:n()),
             feature_id = 1:n()) %>%
      ungroup() %>%
      dplyr::select(feature_id) %>% unlist %>% as.vector()
    t <- t[!t %over% ciri2]
    groups.set[['internal_notflanking']] = t
    
  } else {
    
    cat('>>> get_flanking_introns\n')
    b.is_single_flanking <- is.na(intron.map$upstream.idx) | is.na(intron.map$downstream.idx)
    warning('> circRNA with single intron are discarded. Discarding #circRNA: ', sum(b.is_single_flanking))
    intron.map <- subset(intron.map, !b.is_single_flanking)
    
    groups.set = list()
    t = introns[intron.map$upstream.idx]
    mcols(t)$circ_rna_id = mcols(t)$feature_id = intron.map$circ_rna_id
    groups.set[['upstream']] = t
    
    t = introns[intron.map$downstream.idx]
    mcols(t)$circ_rna_id = mcols(t)$feature_id = intron.map$circ_rna_id
    groups.set[['downstream']] = t
    
    
    t = subset(introns[-c(intron.map$downstream.idx,intron.map$upstream.idx),], name %in% intron.map$gene_id)
    mcols(t)$feature_id = data.frame(gene_id = mcols(t)$name) %>%
      group_by(gene_id) %>%
      mutate(feature_name = paste0(gene_id,'.intron',1:n()),
             feature_id = 1:n()) %>%
      ungroup() %>%
      dplyr::select(feature_id) %>% unlist %>% as.vector()
    t <- t[!t %over% ciri2]
    groups.set[['notflanking']] = t
  }
  
  list(circ_rna_id = intron.map$circ_rna_id, intron.set = groups.set)
}

# MAIN
# import circRNA
ciri2targets <- read_tsv(opt$circRNA.targets,show_col_types = FALSE) %>% janitor::clean_names()
circid2jnc_count <- ciri2targets %>% dplyr::select(circ_rna_id, number_junction_reads)
circid2non_jnc_count <- ciri2targets %>% dplyr::select(circ_rna_id, number_non_junction_reads)
ciri.ref_all <- makeGRangesFromDataFrame(ciri2targets, keep.extra.columns = TRUE)
circID2geneid <- melt(ciri_gene_set(ciri.ref_all)) %>%
  dplyr::select(L1, value) %>%
  dplyr::rename(c('circ_rna_id'='L1','gene_id'='value')) %>% unique

# Import data bases
## genes
organism.genes <- import.bed(opt$genes)
mcols(organism.genes)$gene_id <- mcols(organism.genes)$name
organism.exons <- import.bed(opt$exons)
mcols(organism.exons)$gene_id <- mcols(organism.exons)$name
organism.introns <- import.bed(opt$introns)
mcols(organism.introns)$gene_id <- mcols(organism.introns)$name

ciri.ref <- ciri.ref_all
if(!is.null(opt$exon_boundaries)){
  b.exon_boundaries = is_at_exon_boundary(ciri.ref_all, organism.exons, as.integer(opt$exon_boundaries))
  ciri.ref <- ciri.ref_all[b.exon_boundaries];length(ciri.ref)
}

if(!is.null(opt$collapse_circRNA))
  ciri.ref <- pool_ambiguous_ciri2(ciri.ref, as.integer(opt$collapse_circRNA));length(ciri.ref)


mcols(ciri.ref) <- mcols(ciri.ref) %>% as.data.frame() %>% 
  left_join(circID2geneid %>% group_by(circ_rna_id) %>% summarize(gene_id = paste(gene_id,collapse = ',')), by = 'circ_rna_id')

cat(opt$intron_maxgap,'\n')
intron.index <- create_intron_map(ciri.ref, organism.introns, internal = opt$internal_introns, as.integer(opt$intron_maxgap))
introns.per_ciri2 = make_intron_set(ciri.ref, organism.introns, intron.index,  is.internal = opt$internal_introns)

intron.set <- introns.per_ciri2$intron.set

circ.set <- subset(ciri.ref, circ_rna_id %in% introns.per_ciri2$circ_rna_id)

# select corresponding circRNA
# cat("Export data\n")
print(sapply(intron.set, length))
for(n0 in names(intron.set)){
  export.gff3(intron.set[[n0]], paste0(opt$outprefix,'_ciri2_', n0,'.gff'))
}

export.gff3(circ.set, paste0(opt$outprefix, ifelse(opt$internal_introns, '_internal_ciri2.gff','_ciri2.gff')))
