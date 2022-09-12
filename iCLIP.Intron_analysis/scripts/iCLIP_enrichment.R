suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-b", "--iclip_bedgraph"), type='character',
              help="iCLIP bedgraph (iCount)"),
  
  make_option(c("-g", "--intron_set"), type='character',
              help='mapped introns'),

  make_option(c("-i", "--intron_db"), type='character',
              help="Intron database (gtf)"),
  
  make_option(c("-p", "--pseudo_count"), type='logical', action = 'store_true', default=FALSE,
              help="Add pseudocount"),
  
  make_option(c('-o','--outtable'), type = 'character',
              help="Output table name (tsv)")
)

opt = parse_args(OptionParser(option_list=option_list))
# Rscript iCLIP_enrichment.R -b data/elav_iclip.intron_strict.pooled_cDNA_peaks.bedgraph -g Intron_groups/ -i feature_db/dm6_ensembl96.introns.bed
# opt$iclip_bedgraph = '../input/elav_iclip.intron_strict.pooled_cDNA_peaks.bedgraph'
# opt$intron_set = 'Intron_groups/ectopic_ciri2_downstream.gff'
# opt$intron_db = 'feature_db/dm6_ensembl96.introns.bed'


suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

enrichment <- function(xlsites, features, feature_space = features, pseudocount = FALSE){
  
  # template: https://github.com/tomazc/iCount/blob/2.0.0/iCount/analysis/summary.py
  
  # make sure xlsites are stranded and scores are positive
  assertthat::assert_that(all(strand(xlsites) != '*'), msg = '[debug] missing strand for xlsites')
  assertthat::assert_that(all(score(xlsites) > 0), msg = '[debug] require positive scores for xlsites')
  
  # feature set must be of unique items
  features <- unique(features)
  
  qq <- findOverlaps(xlsites, features) 
  
  # xlsites within feature universe. exclude the rest. is it ?
  total.xlsites <- sum(score(xlsites)) + (as.integer(pseudocount) * length(features))
  
  total.feature_length = sum(width(feature_space))# *2, we count strand specific, so no need for the double-strand correction
  
  xlsites.summed <- sum(score(xlsites)[queryHits(qq)]) + (as.integer(pseudocount) * length(features))
  feature.length <- sum(width(features))
  
  # enrichment[unique(subjectHits(qq)),'enrichment'] = (xlsites.summed/total.xlsites) /
  #   (((feature.length) / (total.feature_length*2)))
  # [ todo ] find significance test
  enrichment = data.frame(xlsites.count = xlsites.summed, 
                 xlsites.total = total.xlsites, 
                 feature.length = feature.length, 
                 total_length = total.feature_length,
                 enrichment = (xlsites.summed/total.xlsites) / (sum(feature.length) / (total.feature_length)),
                 feature.count = length(features),
                 empty.features = length(features[-subjectHits(qq)]))
  
  return(enrichment)
}

f0 <- opt$iclip_bedgraph
xlsites0 = import.bedGraph(f0)
# require strand!
strand(xlsites0) <- ifelse(score(xlsites0) > 0, '+','-')
score(xlsites0) <- abs(score(xlsites0))

features0 = import.gff(opt$intron_set)

Intron_db = import.bed(opt$intron_db)

# Test this function on feature:exons and entire bedgraph, instead of intronic bedgraph
enrichment0 <- data.frame(feature.set = gsub('.gff$','',basename(opt$intron_set)), 
                          enrichment(features0,   
                                     xlsites = xlsites0,
                                     feature_space = Intron_db,
                                     pseudocount = opt$pseudo_count))

write_tsv(enrichment0, paste0(gsub('.tsv$','',opt$outtable),'.tsv'))





