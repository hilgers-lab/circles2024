suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))


args <- commandArgs(T)
# args  = c(5e3, '~/Projects/circCenters/output2/visualize.circRNA/STAR_chimeric/STAR_chimeric/library/w1118_embryos_14_16_R1_RR.Chimeric.out.junction')

max.length = args[1]
## Test data
# STAR junction file

if(tolower(max.length) %in% 'inf'){
  max.length = Inf
} else {
  max.length = as.integer(max.length)
}
tab0 <- suppressMessages(read_tsv(args[2], col_names = FALSE))

tab1 = tab0 %>% 
  # discard fusion transcripts
  filter(X1 == X4 & X3 == X6) %>% 
  # fix coordinates
  mutate(start = ifelse(X3 == '-', X2, X5), end = ifelse(X3 =='-',X5, X2)) %>% 
  # discard inverse mappings after fixing start/stop
  filter(start < end) %>% 
  # prepare GRanges 
  mutate(seqnames = X1, strand = X3, width = X5 - X2)


if(!is.infinite(max.length))
  tab1 <- tab1 %>% filter(abs(width) <= max.length)

# tab.sum 
tab2 <- tab1 %>%
  group_by(seqnames, start, end, strand) %>% summarize(score = n()) %>%
  mutate(name = '.') %>% select(seqnames, start, end, name, score, strand)

cat(paste0("track name=",basename(args[2])," graphType=junctions\n",format_tsv(tab2 %>% arrange(seqnames,start), col_names = FALSE)))
