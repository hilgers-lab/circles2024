# group 1 (split): introns, flanking circRNA. Split into upstream and downstream
# group 2: genome-wide. Note, Discard short introns
# group 3: Introns of circRNA genes, not flanking the circRNA exon
# group 4: introns of functional targets (carrasco2020)
# [open] group 5: differential Exons (DEXseq, carrasco2020), upstream/downstream

suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-m", "--intron_map"), type='character', default = NULL,
              help="Intron map of circRNA to intron mapping (tsv)"),
  
  make_option(c("-i", "--introns"), type='character',
              help="genome annotation (gtf)"),
  
  make_option(c("-c", "--intron_size_thrs"), type='integer', default = 0,
              help="minimum size of an intron required"),
  
  make_option(c('-o','--outdir'), type='character',
              help="directory for intron sets")
)

opt = parse_args(OptionParser(option_list=option_list))


# opt$intron_map = './Features/dm6_ensembl96.circRNA_intron_mapping.tsv'
# opt$introns = './feature_db/dm6_ensembl96.introns.bed'
intron.length.thrs <- opt$intron_size_thrs

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
intron.db <- import.bed(opt$introns)
intron.map <- read_tsv(opt$intron_map)

### Implemented side effects! ###
suppressPackageStartupMessages(library(openxlsx))
## make use of differential table (CIRIquant), instead of circRNA! 
# ciri.ref <- read_tsv('./resources/circRNA.reference.cutoff2.ciri')
# ciriquant.elav <- read_csv('./resources/mutNeuron_NeuronFemale.differential.tsv')
# ciriquant.neurons <- read_csv('./resources/mutNeuron_NeuronFemale.differential.tsv')

targets.apa <- read.xlsx('resources/Table_S1_Carrasco.xlsx')
targets.as <- read.xlsx('resources/Table_S2_Carrasco.xlsx')
### /Implemented side effects! ###

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

fromList2 <- function (input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  rownames(data) <- elements
  return(data)
}


t0 <- fromList2(list(APA = targets.apa$gene_id, AS = targets.as$gene_id.uniq))
t1 <- apply(t0, 1, paste0, collapse = '')
func.map = setNames(c('11' = 'APA/AS', '10' = 'APA','01' = 'AS')[t1], nm = names(t1))


groups.set = list()
# group 1a: introns, flanking circRNA, all circRNA
t = intron.db[intron.map$intron_db.upstream]
mcols(t)$circ_rna_id = intron.map$circ_rna_id
groups.set[['group1a_upstream']] = t
t = intron.db[intron.map$intron_db.downstream]
mcols(t)$circ_rna_id = intron.map$circ_rna_id
groups.set[['group1a_downstream']] = t

# group 1b: introns, flanking circRNA, elav-dependent
t = intron.db[subset(intron.map, circ_rna_id %in% circID.elav)$intron_db.upstream]
mcols(t)$circ_rna_id = subset(intron.map, circ_rna_id %in% circID.elav)$circ_rna_id
groups.set[['group1b_upstream']] = t

t = intron.db[subset(intron.map, circ_rna_id %in% circID.elav)$intron_db.downstream]
mcols(t)$circ_rna_id = subset(intron.map, circ_rna_id %in% circID.elav)$circ_rna_id
groups.set[['group1b_downstream']] = t

# group 2: genome-wide. Note, Discard short introns
intron.db_filt <- subset(intron.db, width > intron.length.thrs)
groups.set[['group2']] = intron.db_filt

# group 3: Introns of circRNA genes, not flanking the exon
introns.circRNAgene.notFlanking = subset(intron.db[-c(intron.map$intron_db.downstream,intron.map$intron_db.upstream),], name %in% intron.map$gene_id)
introns.circRNAgene.notFlanking <- subset(introns.circRNAgene.notFlanking, width > intron.length.thrs)
groups.set[['group3']] = introns.circRNAgene.notFlanking

# group 4: introns of functional targets (carrasco2020)
introns.func_targets = subset(intron.db, name %in% names(func.map) & width > intron.length.thrs)
groups.set[['group4']] = introns.func_targets

print(elementNROWS(groups.set))

for(n0 in names(groups.set)){
  export.gff3(groups.set[[n0]], paste0(opt$outdir,'/','Intronset.', n0,'.gff'))
}

# library(ggplot2)
# library(reshape2)
# library(ggpubr)
# 
# dt0 <- melt(lapply(c(groups.set, 'Intron.db' = intron.db), width))
# 
# p0 <- ggplot(dt0 %>% filter(value < 150), aes(value)) + 
#   geom_histogram(bins = 50) + 
#   geom_vline(xintercept = intron.length.thrs, color = 'grey') +
#   # scale_x_log10() + 
#   xlab('Intron length [nt]') +
#   coord_cartesian(ylim = c(0,0.5e3)) +
#   theme_bw() + 
#   facet_wrap(~L1) 
# 
# p1 <- ggplot(dt0, aes(value)) + 
#   geom_histogram(bins = 50) + 
#   geom_vline(xintercept = intron.length.thrs, color = 'grey') +
#   # scale_x_log10() + 
#   xlab('Intron length [nt]') +
#   # coord_cartesian(ylim = c(0,0.5e3)) +
#   scale_x_log10() +
#   theme_bw() + 
#   facet_wrap(~L1) 
# ggarrange(plotlist = list(p0,p1), common.legend = TRUE, ncol = 1)
# 
# 
# 
# 
# subset(groups.set$group1_upstream, width < 50)
# subset(groups.set$group1_downstream, width < 50) %>% as.data.frame() %>% View
