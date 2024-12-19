Identify and count circRNA using CircSplice.

workfolw:
1. make refFlat file required  for CircSplice from UCSC genome using make_refFlat.R.
2. download snakemake file, modified CircSplice.pl and .pl scripts from https://github.com/GeneFeng/CircSplice.
3. Excute by bash run.sh
4. get count table. 

require packages:
fastqc
trim-galore
STAR2.7.0
bedtools 
samtools
perl

Small modification in CircSplice.pl is to make it excutable in snakemake. 
