suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-b", "--iclip_bedgraph"), type='character',
              help="iCLIP bedgraph (iCount)"),
  make_option('--intron_split', type = 'logical', action = 'store_true', default = FALSE,
               help = "split introns into 5\' and 3\' half"),
  make_option('--intron_filter', type = 'integer', default = 0,
              help = "retain introns with at least X xl-events"),
  make_option(c("-g", "--intron_groups_dir"), type='character',
              help="Folder with prepared sets of introns"),
  
  make_option(c('-o','--outprefix'), type = 'character',
              help="outprefix. One file per intron group (tsv)")
)

opt = parse_args(OptionParser(option_list=option_list))
# Rscript iCLIP_enrichment.R -b data/elav_iclip.intron_strict.pooled_cDNA_peaks.bedgraph -g Intron_groups/ -i feature_db/dm6_ensembl96.introns.bed
opt$iclip_bedgraph = 'elav_iclip.intron_strict.pooled_cDNA_peaks.bedgraph'
opt$intron_groups_dir = 'intron_groups/'
opt$intron_filter = 0
# opt$outprefix

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

intron_filter <- function(xlsites, feature, cutoff){
  is.valid <- rep(F, length(feature))
  qq = findOverlaps(xlsites, feature)
  xl.introns <- sapply(split(score(xlsites)[queryHits(qq)], subjectHits(qq)), sum)
  v <- sapply(xl.introns, function(x, c) x > c, cutoff)
  
  is.valid[as.integer(names(v))] = v
  return(is.valid)
}


relative_position <- function(xlsites, feature, weighted = TRUE){
  
  qq = findOverlaps(xlsites, feature)
  d.set = distance(resize(xlsites, fix = 'start', width = 1)[queryHits(qq)], 
                   resize(feature,fix = 'start', width = 0)[subjectHits(qq)])
  
  rel.pos = data.frame(distance = d.set/(width(feature)[subjectHits(qq)]),
                       score = score(xlsites)[queryHits(qq)])
  if(weighted){
    rel.pos_median = sapply(split(rel.pos, subjectHits(qq)), function(tab)
      (weighted.mean(tab$distance, tab$score))) # /sum(tab$score))))
  } else {
    rel.pos_median = sapply(split(rel.pos$distance, subjectHits(qq)), median, na.rm = TRUE)
  }
  
  xlsite.events = sapply(split(queryHits(qq), subjectHits(qq)), length)
  xlsite.count = sapply(split(abs(score(xlsites)[queryHits(qq)]), subjectHits(qq)), sum)
  
  df = matrix(0, ncol = 4, nrow = length(feature))
  
  df[unique(subjectHits(qq)),1] = xlsite.events
  df[unique(subjectHits(qq)),2] = xlsite.count
  df[,3] = width(feature)
  df[unique(subjectHits(qq)),4] = rel.pos_median
  df[!1:nrow(df) %in% subjectHits(qq),4] = NA
  
  df = as.data.frame(df)
  colnames(df) = c('xl.event','xl.count','feature.length','relative_position')
  df = cbind(feature_id = mcols(feature)$feature_id, df)
  return(df)
}

f0 <- opt$iclip_bedgraph
xlsites0 = import.bedGraph(f0)
# require strand!
strand(xlsites0) <- ifelse(score(xlsites0) > 0, '+','-')
score(xlsites0) <- abs(score(xlsites0))

f0 = list.files(opt$intron_groups_dir, pattern = '.*(downstream|geneset|notflanking|upstream).gff$', full.names = TRUE)
names(f0) = gsub('(.*)_(downstream|geneset|notflanking|upstream).gff','\\1.\\2', basename(f0))
feature.set0 = lapply(f0, import.gff)


void <- assertthat::assert_that(all(sapply(lapply(feature.set0, mcols),function(x, ref) ref %in% colnames(x), 'feature_id')),
                        msg = 'feature_id in Intron group gff missing')


# Intron filter - make upsteam/downstream symmetric, isolated for other groups
intron_filter.xlcounts = 0 # opt$intron_filter 
{
  intron_valid.set <- lapply(feature.set0, intron_filter, xlsites = xlsites0, intron_filter.xlcounts)
  
  # overwrite for upstream/downstream with 'any'
  # pair upstream/downstream introns. If any, take it
  intron_pair.labels = data.frame(label = names(intron_valid.set),
                                  set = gsub('(.*).(geneset|downstream|upstream|notflanking)','\\1',names(intron_valid.set)),
                                  feature = gsub('(.*).(geneset|downstream|upstream|notflanking)','\\2',names(intron_valid.set))) %>%
    filter(feature %in% c('downstream','upstream')) 
  
  # manipulate vectors directly!
  void <- lapply(split(intron_pair.labels$label, intron_pair.labels$set), function(labels, valid) {
    v = apply(do.call('cbind', valid[labels]), 1, function(x) any(x)) 
    for(lab in labels){
      valid[[lab]] = v
    }
    valid[labels]
  }, intron_valid.set)
  
  for(n0 in names(void)){
    for(n1 in names(void[[n0]])){
      intron_valid.set[[n1]] = void[[n0]][[n1]]
    }
  }
}

# intron split operation, do this for all intron sets
feature.set1 <- mapply(function(ref, b) {ref[b]}, feature.set0, intron_valid.set)

d0 = lapply(feature.set1, relative_position, xlsites = subset(xlsites0, score > -1), weighted = TRUE)
dt0 = melt(d0, id.vars = colnames(d0[[1]]))
dt0 <- dt0 %>% mutate(feature = gsub('(.*)\\.(upstream|downstream|notflanking|geneset)','\\1', L1),
               group = gsub('(.*)\\.(upstream|downstream|notflanking|geneset)','\\2', L1))

dt.set <- split(dt0 %>% select(-L1), dt0$L1)

# for(n0 in names(dt.set)){
#   write_tsv(dt.set[[n0]], gsub('__+','_',paste0(opt$outprefix,'_',gsub('.gff$','',n0),'.tsv')))
# }

library(ggplot2)


dt0 <- melt(dt.set, id.vars = colnames(dt.set[[1]])) %>% mutate(group = factor(group, c('upstream','downstream','notflanking','geneset')))

p.pos_distr = ggplot(dt0, aes(x = relative_position)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50) +
  geom_density(aes(relative_position, color = group), show.legend = FALSE) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,4)) + # , ylim = c(0,3)) +
  xlab('relative position (5\' intron)') + 
  geom_vline(xintercept = c(0,1), color = 'darkgrey') +
  ggtitle('Centrality - iCLIP in intron') +
  theme_bw() +
  facet_grid(feature~group) 
print(p.pos_distr)




# Intron splitting analysis
feature.set5prime = lapply(feature.set1, function(ref) {resize(ref, width = width(ref)/2, fix = 'start')})
feature.set3prime = lapply(feature.set1, function(ref) {resize(ref, width = width(ref)/2, fix = 'end')})

feature.list = list('set5prime' = feature.set5prime, 'set3prime' = feature.set3prime)
for(set_name in names(feature.list)){
  set = feature.list[[set_name]]
  set0 <- set[!grepl('geneset|notflanking',names(set)) ]
  
  d0 = lapply(set0, relative_position, xlsites = subset(xlsites0, score > -1), weighted = FALSE)
  dt0 = melt(d0, id.vars = colnames(d0[[1]]))
  dt0 <- dt0 %>% mutate(feature = gsub('(.*)\\.(upstream|downstream|notflanking|geneset)','\\1', L1),
                        group = gsub('(.*)\\.(upstream|downstream|notflanking|geneset)','\\2', L1))
  
  dt.set <- split(dt0 %>% select(-L1), dt0$L1)
  
  # for(n0 in names(dt.set)){
  #   write_tsv(dt.set[[n0]], gsub('__+','_',paste0(opt$outprefix,'_',gsub('.gff$','',n0),'.tsv')))
  # }
  
  library(ggplot2)
  
  
  dt0 <- melt(dt.set, id.vars = colnames(dt.set[[1]])) %>% mutate(group = factor(group, c('upstream','downstream','notflanking','geneset')))
  
  p.pos_distr = ggplot(dt0, aes(x = relative_position)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50) +
    geom_density(aes(relative_position, color = group), show.legend = FALSE) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,4)) + # , ylim = c(0,3)) +
    xlab('relative position (5\' intron)') + 
    geom_vline(xintercept = c(0,1), color = 'darkgrey') +
    ggtitle(paste0('Centrality - iCLIP in intron', '-',set_name)) +
    theme_bw() +
    facet_grid(feature~group) 
  print(p.pos_distr)
}
