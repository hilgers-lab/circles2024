suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-b", "--iclip_bedgraph"), type='character',
              help="iCLIP bedgraph (iCount)"),
  
  make_option(c("-g", "--intron_set"), type='character',
              help='mapped introns'),
  
  make_option(c('-o','--outfile'), type = 'character',
              help="outfile (tsv)")
)

opt = parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))


centrality <- function(xlsites, feature){
  
  qq = findOverlaps(xlsites, feature)
  d.set = distance(resize(xlsites, fix = 'start', width = 1)[queryHits(qq)], 
                   resize(feature, fix = 'start', width = 0)[subjectHits(qq)])
  rel.pos = d.set/(width(feature)[subjectHits(qq)])
  
  rel.pos_median = sapply(split(rel.pos, subjectHits(qq)), median)
  
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

feature.set0 = import.gff(opt$intron_set)

void <- assertthat::assert_that(any(colnames(mcols(feature.set0)) %in% 'feature_id'),
                        msg = 'feature_id in Intron group gff missing')

d0 = centrality(xlsites = xlsites0, feature.set0)

write_tsv(d0, paste0(gsub('.tsv$','',opt$outfile),'.tsv'))

