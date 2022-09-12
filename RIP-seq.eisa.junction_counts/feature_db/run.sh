module load R/4.0.3 bedtools2

Rscript make_feature_db.R # make sure it's R>=4.0.0

bedtools subtract -s -a dm6_ensembl96.genes.bed -b dm6_ensembl96.exonic_database.bed > dm6_ensembl96.intronic_database.bed
