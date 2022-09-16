set -e
# Requirements
# module load samtools/1.12 deeptools/3.5.0 subread/2.0.0 

read directory config params <<< "$@"

snakemake $params --configfile $config -d $directory --latency-wait 180 --keep-going
