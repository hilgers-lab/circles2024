######load packages needed##########
packages.toLoad <- c("GenomicRanges", "BSgenome.Dmelanogaster.UCSC.dm6", "plyranges", "TxDb.Dmelanogaster.UCSC.dm6.ensGene", "AnnotationDbi", "readxl", "readr", 
        "dplyr", "stringr", "stats", "tidyverse", "shiny", "Biobase")
loaded <- (.packages())
load_all <- function(list) {
  while (any(packages.toLoad %in% loaded==FALSE)) {
    lapply(list[!list %in% loaded], require, character.only = T)
    loaded <- (.packages())
  }
}
load_all(packages.toLoad)
loaded <- (.packages())
all(packages.toLoad %in% loaded)
rm(loaded, packages.toLoad)

#########check what is the bed-refFlat file############
setwd("/data/hilgers/group2/Shi/2024_circSplice")
example <- read.table("doc/USCS-refFlat_dm6.txt", sep = "\t")

# refFlat_dm6 was downloaded from UCSC 
path_ref <- "/data/hilgers/group2/Shi/2024_circSplice/doc/bed-refFlat_dm6.txt"
dm6_ref <- read.table(path_ref, sep = "\t")
dm6_ref_circSplice <- dm6_ref %>%
  select(-c(V1, V12, V15, V16)) %>% #kick out unwanted column
  rename(TranscriptID = V2, Chrom = V3, strand = V4, Start = V5, End = V6, 
         cdsStart = V7, cdsEnd = V8, exonCount = V9, exonStarts = V10, exonEnds = V11, geneID = V13, type = V14) %>%
  # mutate(Chrom = gsub("chr", "", Chrom)) %>%
  filter(!grepl("FBti", geneID)) %>% #kick out transposon 
  mutate(startMinus2 = Start -2, endPlus2 = End + 2,
         type = ifelse(type == "cmpl", "mRNA", "lncRNA"),
         TxID = TranscriptID, chrom = Chrom, Strand = strand,
         numExon = exonCount) %>%
  select(c(Chrom, startMinus2, endPlus2, TranscriptID, exonCount, strand, geneID, TxID, chrom, Strand, Start, End, cdsStart, cdsEnd,
           numExon, exonStarts, exonEnds, type)) %>%
  mutate(startMinus2 = ifelse(startMinus2 < 0, 0, startMinus2))

write.table(dm6_ref_circSplice, "doc/bed_refFlat_dm6_ucsc.txt", sep = "\t", col.names = F, row.names = F, quote = F)


  