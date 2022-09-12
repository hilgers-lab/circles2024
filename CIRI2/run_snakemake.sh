module load slurm bwa/0.7.17 seqtk/1.2 cutadapt/2.10

read directory config params <<< "$@"

mkdir -p $directory/cluster_log
snakemake -d $directory --configfile $config $params --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log" --latency-wait 500 --keep-going
