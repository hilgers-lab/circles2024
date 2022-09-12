library(readr)
library(dplyr)
library(optparse)


option_list <- list(
  make_option(c("-s", "--samplesheet"), type = 'character',
              help="samplesheet à la CIRIquant, without paths"),
  make_option(c("-p", "--path"), type = 'character',
              help="source path to get gtf from"),
  make_option(c("-o", "--outFile"), type="character",
              help="name of output table, i.e. CIRIquant samplesheet"),
  make_option("--genes", action='store_true', default = FALSE,
              help = "Sample sheet should be prepared for gene expression"),
  # make_option("--use_label", action='store_true', default = FALSE,
  #             help = "Report \'label\' column instead of sample \'sample\'"),
  make_option("--use_batch", action='store_true', default = FALSE,
              help = "Use \'batch\' column for subject columns. Ignore if \'--genes\' is set")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# opt = list()
# opt$samplesheet = '../../elav_mutants.circRNA/CIRIquant.pooled.libtype2/SAMPLESHEETS/NeuronAll_nonNeuron.tsv'
# opt$path = '../../elav_mutants.circRNA/CIRIquant.pooled.libtype2/CIRIquant/circles_w1118_heads/'
# opt$genes = TRUE
# opt$use_batch = TRUE

samplesheet = read_tsv(opt$samplesheet)
proj.path = opt$path
outtab_name = opt$outFile

pattern0 = '.gtf$'
if(opt$genes){
  pattern0 = paste0('_out',pattern0)
  proj.path = file.path(proj.path, 'gene')
}

fset = normalizePath(list.files(proj.path, pattern = pattern0, full.names = TRUE))

if(length(fset) == 0)
  stop("No files found. Please check the path")

fileset = tibble(sample = gsub(pattern0,'', basename(fset)), path = fset)
if(opt$use_batch){
  b = 'batch' %in% colnames(samplesheet)
  cat("Using batch:", b, '\n')
  opt$use_batch = b
  if(!b){
    cat("Skipping batch. No \'batch\' column available\n")
  }
}

  
if(opt$use_batch){
  tabout = samplesheet %>% left_join(fileset, by = c('sample')) %>%
    mutate(subject = as.numeric(as.factor(batch))) %>%
    select(sample, path, condition, subject)
} else {
  tabout = samplesheet %>% left_join(fileset, by = c('sample')) %>%
    select(sample, path, condition)
}
tabout = tabout %>% mutate(condition = c('C','T')[as.integer(relevel(as.factor(condition),unique(condition)[1]))])

if(opt$genes){
  tabout = tabout %>%
    select(sample, path)
}
if(any(is.na(tabout$path)))
  stop('Paths not found. samplenames must match between samplesheet and CIRIquant directory')

write_tsv(tabout, path = outtab_name, col_names = FALSE)
