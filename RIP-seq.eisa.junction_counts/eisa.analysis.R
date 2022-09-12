library(eisaR)
library(readr)
library(rtracklayer)
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)

parse_eisa <- function(eisa.obj, outprefix = NULL, map.gid2symbol = NULL) {
  #> filtering quantifyable genes...keeping 11759 from 20288 
  tab0 <- eisa.obj$tab.ExIn %>% rownames_to_column('gene_id') %>% 
    left_join(as.data.frame(eisa.obj$contrasts) %>% rownames_to_column('gene_id'))
  
  if(!is.null(map.gid2symbol)){
    tab0 <- tab0 %>% mutate(gene_name = map.gid2symbol[gene_id]) 
  }
  if(!is.null(outprefix)){
    write_tsv(tab0, paste0(outprefix,'.exon_intron.tsv'))
    write_rds(list('eisa.obj'=eisa.obj, 'tab.ExIn' = tab0), paste0(outprefix,'.rds'))
  }
  return(list('eisa.obj'=eisa.obj, 'tab.ExIn' = tab0))
}

dm6.genes <- import.gff('data/resources/dm6_ensembl96.gtf', feature.type = 'gene')
gid2symbol = mcols(dm6.genes)[c('gene_id','gene_name')] %>% deframe


outdir0='eisa.quantification.junction_counts/eisa_analysis.results'
assertthat::assert_that(dir.exists(outdir0))

cntEx_raw <- read_tsv('./Quantification.junctions/featureCounts/featureCounts.mode_paired.tsv.jcounts', comment = '#'); colnames(cntEx_raw) = basename(colnames(cntEx_raw))
cntEx <- cntEx_raw %>%
  dplyr::select(-SecondaryGenes,-Site1_chr,-Site1_location,-Site1_strand,-Site2_chr,-Site2_location,-Site2_strand) %>%
  melt  %>% 
  dcast(formula = PrimaryGene ~ variable, value.var = 'value', fun.aggregate = sum)
cntEx <- cntEx %>% filter(!is.na(PrimaryGene)) %>% as.data.frame() %>% column_to_rownames('PrimaryGene')

cntIn_raw <- read_tsv('./Quantification.junctions/featureCounts.eij/featureCounts.mode_paired.tsv', comment = '#'); colnames(cntIn_raw) = basename(colnames(cntIn_raw))
cntIn <- (cntIn_raw %>% column_to_rownames('Geneid')) %>% dplyr::select(-Chr,-Start,-End,-Strand,-Length)
cntGb_raw <- read_tsv('./Quantification.junctions/featureCounts/featureCounts.mode_paired.tsv', comment = '#'); colnames(cntGb_raw) = basename(colnames(cntGb_raw))
cntGb <- (cntGb_raw %>% column_to_rownames('Geneid')) %>% dplyr::select(-Chr,-Start,-End,-Strand,-Length)

# Minimum expression filter in Genebody
gid.keep = intersect(rownames(cntGb)[rowMeans(cntGb) > 10], rownames(cntEx))

# sum(keep)/length(keep)
cntEx = cntEx[gid.keep,]
cntIn = cntIn[gid.keep,]
cntGb = cntGb[gid.keep,]

## checkout exon-intron ratios
# dt0 <- data.frame(Exons = colSums(cntEx), Introns = colSums(cntIn), Genebody=colSums(cntGb)) %>% 
#   rownames_to_column('sample') %>% 
#   filter(!grepl('_HA_', sample))
# 
# dt0 <- dt0 %>% mutate(group = gsub('.*(elav_FLAG|heads|w1118_FLAG).*(input|IP).*','\\1',sample),
#                       ip = gsub('.*(elav_FLAG|heads|w1118_FLAG).*(input|IP).*','\\2',sample),
#                       rep = gsub('_|R','',gsub('.*(_[[:digit:]]|R[[:digit:]]).*','\\1', sample)))
# df0 <- melt(dt0)

## exon-intron ratios
# ggplot(dt0, aes(Introns/Genebody, Exons/Genebody)) + geom_point(aes(color = group, pch = ip), size = 3) +
#   coord_cartesian(xlim = c(0,.5), ylim = c(0,.5)) +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_bw()

# ggplot(dt0 %>% mutate(ratio = (Introns)/(Exons + Introns)),aes(group)) + geom_bar(aes(weight = ratio, fill = rep), position = 'dodge') +
#   coord_cartesian(ylim = c(0,1)) +
#   scale_fill_viridis_d() +
#   geom_hline(mapping = aes(yintercept = median(ratio))) +
#   ylab('Intron fraction') +
#   facet_wrap(~ip) +
#   theme_bw()


Rex <- cntEx
Rin <- cntIn


# create condition factor (contrast will be TN - ES)

f0 <- list.files('./Samplesheets/', pattern = 'samplesheet.tsv$', full.names = TRUE)
names(f0) <- gsub('(.*).samplesheet.tsv','\\1',basename(f0))
samplesheets = lapply(f0, read_tsv)

# method_set = list('Gaidatzis2015' = c(method = 'Gaidatzis2015'),
#                   'sizeFactor_intron' = c(sizeFactor = 'intron', method = NULL),
#                   'sizeFactor_exon' = c(sizeFactor = 'exon', method = NULL),
#                   'sizeFactor_individual' = c(sizeFactor = 'individual', method = NULL))
#
# res.set <- list()
# for(method0 in names(method_set)) {
#   for(sheet0 in names(samplesheets)){
#     dir.create(paste0(outdir0,'_',method0), showWarnings = FALSE)
#     outprefix0 = paste0(outdir0,'_',method0,'/','eisa.', sheet0)
#     
#     samplesheet = samplesheets[[sheet0]]
#     cond0 <- as.factor(samplesheet$type)
#     
#     cat('>> Outprefix:', outprefix0, '\n')
#     if(!grepl('sizeFactor',method0)){
#       eisa.obj <- runEISA(Rex[,samplesheet$sample], Rin[,samplesheet$sample], cond0, method = method_set[[method0]]['method'])
#     } else if (method0 == 'Gaidatzis2015'){
#       eisa.obj <- runEISA(Rex[,samplesheet$sample], Rin[,samplesheet$sample], cond0, sizeFactor = method_set[[method0]][['sizeFactor']])
#     } else {
#       eisa.obj <- runEISA(Rex[,samplesheet$sample], Rin[,samplesheet$sample], cond0, method = NULL, unlist(method_set[[method0]]))
#     }
#     res.set[[paste0(sheet0,'_',method0)]] = parse_eisa(eisa.obj, outprefix0, map.gid2symbol = gid2symbol)
#   }
# }
#


# modelSamples = FALSE,
# geneSelection = 'filterByExpr'
# statFramework = 'QLF'
# effects = 'predFC'
# pscnt = 2
# sizeFactor = 'individual', 
# recalcNormFactAfterFilt = TRUE
# recalcLibSizeAfterFilt = FALSE
colnames(Rin) <- gsub('.eij_CIGAR_M.bam','.filtered.bam', colnames(Rin))
res.set <- list()
method0 = 'libSizeIndividual' # libsize individual
for(sheet0 in names(samplesheets)){
  dir.create(paste0(outdir0,'_',method0), showWarnings = FALSE)
  outprefix0 = paste0(outdir0,'_',method0,'/','eisa.', sheet0)
  
  samplesheet = samplesheets[[sheet0]]
  cond0 <- as.factor(samplesheet$type)
  
  cat('>> Outprefix:', outprefix0, '\n')
  eisa.obj <- runEISA(Rex[,samplesheet$sample], Rin[,samplesheet$sample], cond0, 
                      modelSamples = FALSE,
                      geneSelection = 'filterByExpr',
                      statFramework = 'QLF',
                      effects = 'predFC',
                      pscnt = 2,
                      sizeFactor = 'individual', # use individual, comm fri221021
                      recalcNormFactAfterFilt = TRUE,
                      recalcLibSizeAfterFilt = FALSE
                      )
  res.set[[paste0(sheet0,'_',method0)]] = parse_eisa(eisa.obj, outprefix0, map.gid2symbol = gid2symbol)
}


