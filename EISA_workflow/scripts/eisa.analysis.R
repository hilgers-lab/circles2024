suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-e", "--eej_count"), type='character', default = NULL,
              help="count matrix for exons, feature counts junctions files"),
  
  make_option(c("-i", "--eij_count"), type='character',
              help="count matrix for introns, feature counts tsv"),
  
  make_option(c("-g", "--gene_count"), type='character',
              help="count matrix for genes, feature counts tsv"),

  make_option(c("--gtf"), type='character',
              help="genome annotation (gtf)"),
  
  make_option(c("-s", "--samplesheet"), type='character',
              help="Samplesheet for eisa (tsv)"),
  
  make_option(c('-o','--outprefix'), type='character',
              help="outprefix for intron sets")
)

opt = parse_args(OptionParser(option_list=option_list))


suppressPackageStartupMessages(library(eisaR))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

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
  # return(list('eisa.obj'=eisa.obj, 'tab.ExIn' = tab0))
  return(paste0(outprefix,'.exon_intron.tsv'))
}

# dm6.genes <- import.gff('data/resources/dm6_ensembl96.gtf', feature.type = 'gene')


dm6.genes <- import.gff(opt$gtf, feature.type = 'gene')
gid2symbol = mcols(dm6.genes)[c('gene_id','gene_name')] %>% deframe

# case 1 - exon count table, intron required
opt$eej_count 
opt$eij_count 
opt$gene_count 


# cntEx_raw <- read_tsv('./Quantification.junctions/featureCounts/featureCounts.mode_paired.tsv.jcounts', comment = '#'); colnames(cntEx_raw) = basename(colnames(cntEx_raw))
cntEx_raw <- read_tsv(opt$eej_count, comment = '#', show_col_types = FALSE); 
colnames(cntEx_raw) = basename(colnames(cntEx_raw))
cntEx <- cntEx_raw %>%
  dplyr::select(-SecondaryGenes,-Site1_chr,-Site1_location,-Site1_strand,-Site2_chr,-Site2_location,-Site2_strand) %>%
  melt  %>% 
  dcast(formula = PrimaryGene ~ variable, value.var = 'value', fun.aggregate = sum)
cntEx <- cntEx %>% filter(!is.na(PrimaryGene)) %>% as.data.frame() %>% column_to_rownames('PrimaryGene')

# cntIn_raw <- read_tsv('./Quantification.junctions/featureCounts.eij/featureCounts.mode_paired.tsv', comment = '#'); colnames(cntIn_raw) = basename(colnames(cntIn_raw))
cntIn_raw <- read_tsv(opt$eij_count, comment = '#', show_col_types = FALSE); 
colnames(cntIn_raw) = basename(colnames(cntIn_raw))
cntIn <- (cntIn_raw %>% column_to_rownames('Geneid')) %>% dplyr::select(-Chr,-Start,-End,-Strand,-Length)
# cntGb_raw <- read_tsv('./Quantification.junctions/featureCounts/featureCounts.mode_paired.tsv', comment = '#'); colnames(cntGb_raw) = basename(colnames(cntGb_raw))

cntGb_raw <- read_tsv(opt$gene_count, comment = '#', show_col_types = FALSE); 
colnames(cntGb_raw) = basename(colnames(cntGb_raw))
cntGb <- (cntGb_raw %>% column_to_rownames('Geneid')) %>% dplyr::select(-Chr,-Start,-End,-Strand,-Length)

# Minimum expression filter in Genebody
gid.keep = intersect(rownames(cntGb)[rowMeans(cntGb) > 10], rownames(cntEx))

# sum(keep)/length(keep)
cntEx = cntEx[gid.keep,]
cntIn = cntIn[gid.keep,]
cntGb = cntGb[gid.keep,]

Rex <- cntEx
Rin <- cntIn


# create condition factor (contrast will be TN - ES)
# f0 <- list.files('./Samplesheets/', pattern = 'samplesheet.tsv$', full.names = TRUE)
# names(f0) <- gsub('(.*).samplesheet.tsv','\\1',basename(f0))
# samplesheets = lapply(f0, read_tsv)

sheet0 = read_tsv(opt$samplesheet, show_col_types = FALSE)
print(sheet0)
colnames(Rin) <- gsub('.eij_CIGAR_M.bam','.filtered.bam', colnames(Rin))
method0 = 'libSizeIndividual' # libsize individual
cond0 <- factor(sheet0$type, levels = unique(sheet0$type))
print(cond0)

# modelSamples = FALSE,
# geneSelection = 'filterByExpr'
# statFramework = 'QLF'
# effects = 'predFC'
# pscnt = 2
# sizeFactor = 'individual', 
# recalcNormFactAfterFilt = TRUE
# recalcLibSizeAfterFilt = FALSE
eisa.obj <- runEISA(Rex[,sheet0$sample], Rin[,sheet0$sample], cond0, 
                    modelSamples = FALSE,
                    geneSelection = 'filterByExpr',
                    statFramework = 'QLF',
                    effects = 'predFC',
                    pscnt = 2,
                    sizeFactor = 'individual', # use individual, comm fri221021
                    recalcNormFactAfterFilt = TRUE,
                    recalcLibSizeAfterFilt = FALSE
)

outprefix0 = opt$outprefix
dir.create(dirname(outprefix0), recursive = TRUE, showWarnings = FALSE)
cat('>> Outprefix:\n')
parse_eisa(eisa.obj, outprefix0, map.gid2symbol = gid2symbol)


