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


profile_tab <- function(xlsites, feature){
  
  qq = findOverlaps(xlsites, feature)
  d.set = distance(resize(xlsites, fix = 'start', width = 1)[queryHits(qq)], 
                   resize(feature, fix = 'start', width = 0)[subjectHits(qq)])
  rel.pos = data.frame(relative_position = d.set/(width(feature)[subjectHits(qq)]), 
                       xl.counts = score(xlsites)[queryHits(qq)],
                       subjectHits = as.character(subjectHits(qq))) %>% 
    group_by(subjectHits) %>% 
    mutate(i = 1:n(),
           xl.log_density = -log10(xl.counts/sum(xl.counts))) %>% 
    ungroup
  list(table = rel.pos, 
       features_empty = feature[-subjectHits(qq)])
}

# opt$iclip_bedgraph = '../input/elav_iclip.intron_strict.pooled_cDNA_peaks.bedgraph'
# opt$intron_set = '../iCLIP.introns_analysis.elav_spliced/Intron_groups/elav_down_dexseq_downstream.gff'
# opt$intron_set = '../iCLIP.introns_analysis.elav_spliced/Intron_groups/elav_down_dexseq_upstream.gff'
# opt$intron_set = '../iCLIP.introns_analysis.elav_spliced/Intron_groups/elav_up_dexseq_downstream.gff'
# opt$intron_set = '../iCLIP.introns_analysis.elav_spliced/Intron_groups/elav_up_dexseq_upstream.gff'

f0 <- opt$iclip_bedgraph
xlsites0 = import.bedGraph(f0)
# require strand!
strand(xlsites0) <- ifelse(score(xlsites0) > 0, '+','-')
score(xlsites0) <- abs(score(xlsites0))

feature.set0 = import.gff(opt$intron_set)

void <- assertthat::assert_that(any(colnames(mcols(feature.set0)) %in% 'feature_id'),
                                msg = 'feature_id in Intron group gff missing')

dobj = profile_tab(xlsites = xlsites0, feature.set0)
d0 = dobj$table
d.empty = dobj$features_empty

write_tsv(d0, paste0(gsub('.tsv$','',opt$outfile),'.tsv'))
export.gff3(d.empty, paste0(gsub('.tsv$','',opt$outfile),'.empty_features.gff'))

