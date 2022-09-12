set -e
module load samtools/1.12 deeptools/3.5.0 subread/2.0.0 slurm

read directory config params <<< "$@"

snakemake $params --configfile $config -d $directory --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log" --latency-wait 180 --keep-going
