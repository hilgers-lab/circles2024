suppressPackageStartupMessages(library(optparse))


option_list = list(
  make_option(c("-g", "--gtf"), type='character',
              help="genome annotation (gtf)"),

  make_option(c("--transcript_table"), type='character', default = NULL,
              help="Table with target transcript_ids. Requires column \'transcript_id\' (default %default)"),

  make_option(c("--biotypes"), type='character', default = 'protein_coding',
              help="permitted biotypes, tag: \'gene_biotype\', comma-separated (default %default)"),
  make_option(c("--exon_size_min"), type='integer', default = 0,
              help="Mnimimum span of an exon (default %default)"),

  make_option(c('--out_genes'), type = 'character',
              help="outprefix for circRNA-intron mapping"),
  make_option(c('--out_exons'), type = 'character',
              help="outprefix for circRNA-intron mapping")
)


opt = parse_args(OptionParser(option_list=option_list))
# opt$gtf = '../resources/dm6_ensembl96.gtf'
# opt$transcript_table = '../input/RNAseq.tx_expressed.tsv'
# opt$biotypes = NULL
# opt$exon_size_min = 0

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(tibble))

biotypes.valid <- opt$biotypes
if(!is.null(biotypes.valid))
  biotypes.valid <- strsplit(opt$biotypes,',')[[1]]

if(!is.null(opt$transcript_table))
  tx.expressed <- read_tsv(opt$transcript_table)

genes <- import.gff(opt$gtf, feature.type = 'gene')
void <- assertthat::assert_that(all(biotypes.valid %in% mcols(genes)$gene_biotype),
                                msg = paste0("\ngene_biotype not found: \'", paste(biotypes.valid[!biotypes.valid %in% mcols(genes)$gene_biotype],collapse = ', '),'\'',
                                             '\n Available biotypes: ', paste(sort(unique(mcols(genes)$gene_biotype)),collapse = ', ')))

gid2biotype <- (mcols(genes)[,c('gene_id','gene_biotype')]) %>% as.data.frame()

transcripts <- import.gff(opt$gtf, feature.type = 'transcript') 
tx2gid <- mcols(transcripts)[c('transcript_id','gene_id')] %>% deframe

txdb0 <- makeTxDbFromGFF(opt$gtf)

## extract exons
exonic_db = unlist(GRangesList(exonsBy(txdb0, by = 'tx', use.names = TRUE)))
mcols(exonic_db)$transcript_id = names(exonic_db)
intronic_db <- intronsByTranscript(txdb0, use.names = TRUE)
tss <- resize(transcripts, fix = 'start', width = 1) %>% unique
tes <- resize(transcripts, fix = 'end', width = 1) %>% unique

# Expressed Transcript filter
if(!is.null(opt$transcript_table)){
  cat('>>> Subsetting by expressed isoforms\n')
  
  tx.tab <- transcripts %>% as.data.frame()
  
  stats = list()
  stats[['transcripts']] = c('total' = nrow(tx.tab),
                           tx.tab %>% summarize(expressed = sum(transcript_id %in% tx.expressed$transcript_id)))
  stats[['genes']] <- c('total' = length(unique(tx.tab$gene_id)),
                      'expressed' = tx.tab %>% filter(transcript_id %in% tx.expressed$transcript_id) %>% dplyr::select(gene_id) %>% unique %>% summarize(genes = n()))
  stats[['exons']] = c('total' = length(exonic_db), 'expressed' = length(exonic_db[names(exonic_db) %in% tx.expressed$transcript_id,]))
  stats[['tss']] = c('total'= length(tss), tss %>% as.data.frame() %>% summarize(expressed = sum(transcript_id %in% tx.expressed$transcript_id)))
  stats[['tes']] = c('total' = length(tes), tes %>% as.data.frame() %>% summarize(expressed = sum(transcript_id %in% tx.expressed$transcript_id)))
  print(do.call(rbind, stats)[,c(2,1)])
  
  transcripts <- transcripts %>% subset(transcript_id %in% tx.expressed$transcript_id)
  tx2gid <- mcols(transcripts)[c('transcript_id','gene_id')] %>% deframe
  
  exonic_db <- subset(exonic_db, transcript_id %in% tx.expressed$transcript_id)
  intronic_db <- intronic_db[tx.expressed$transcript_id]
  tss <- tss %>% subset(transcript_id %in% tx.expressed$transcript_id)
  tes <- tes %>% subset(transcript_id %in% tx.expressed$transcript_id)
}

## intronic utr-filter ##
mcols(exonic_db) <- mcols(exonic_db) %>% as.data.frame() %>%
  mutate(transcript_id = names(exonic_db)) %>%
  mutate(gene_id = tx2gid[transcript_id]) %>%
  left_join(gid2biotype, by = 'gene_id')

# exclude first TSS, first TES
aTSS <- tss
aTSS <- aTSS[!aTSS %over% resize(genes, width = 1, fix = 'start')]
aTES <- tes
aTES <- aTES[!aTES %over% resize(genes, width = 1, fix = 'end')]

# discard exons overlapping UTRs and intron
b.feature_select <- !((exonic_db %over% aTSS & exonic_db %within% intronic_db) | 
                        (exonic_db  %over% aTES & exonic_db %within% intronic_db)) # ;sum(b.feature_select)/length(b.feature_select)

exonic_db.filtered <- subset(exonic_db, b.feature_select & width >= opt$exon_size_min)

if(!is.null(biotypes.valid)){
  exonic_db.filtered <- subset(exonic_db.filtered, gene_biotype %in% biotypes.valid)
}
exonic_db.reduced <- exonic_db.filtered %>% reduce

## Gene definition to be max extend of _expressed_ transcript, instead of all annotated tx
gr.set = unlist(reduce(split(transcripts, mcols(transcripts)$gene_id)))
mcols(gr.set)$gene_id = names(gr.set)
genes.set = subset(gr.set, gene_id %in% mcols(exonic_db.filtered)$gene_id)

## Annotate reduced exon set exonic_db.reduced
qi = findOverlaps(exonic_db.reduced, exonic_db)
gid.set = sapply(split(mcols(exonic_db)$gene_id[subjectHits(qi)], queryHits(qi)), function(v) paste0(sort(unique(v)),collapse = ','))
names(exonic_db.reduced)[as.integer(names(gid.set))] = gid.set

### Others
# genes.set = subset(genes(txdb0), gene_id %in% mcols(exonic_db.filtered)$gene_id)
#####
## debug output
# report exonic_db a) selected and b) discard
# dir.create('exon_database.debug')
# export.gff(exonic_db.reduced, 'exon_database.debug/exons_reduced.selected.gff')
# export.gff(exonic_db, 'exon_database.debug/exons.expressed.gff')
# export.gff(exonic_db[!b.feature_select], 'exon_database.debug/exons.discarded.gff')
# export.bed(genes.set, paste0('exon_database.debug/genes.selected.bed'))
# Problematic X   4075743-4076426 
### 

export.bed(exonic_db.reduced, ifelse(endsWith(opt$out_exons,'.bed'), opt$out_exons, paste0(opt$out_exons,'.bed')))
export.bed(genes.set, ifelse(endsWith(opt$out_genes,'.bed'), opt$out_genes, paste0(opt$out_genes,'.bed')))
