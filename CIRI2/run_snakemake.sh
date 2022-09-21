# Requirements
# module load bwa/0.7.17 seqtk/1.2 cutadapt/2.10

read directory config params <<< "$@"

mkdir -p $directory/cluster_log
snakemake -d $directory --configfile $config $params --keep-going
