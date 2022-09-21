set -e
# Requirements
module load samtools/1.12 deeptools/3.5.0 subread/2.0.0 R/4.0.3 slurm bedtools2/2.27.0 

read directory config params <<< "$@"

snakemake $params --configfile $config -d $directory --latency-wait 180 --keep-going  --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log"
