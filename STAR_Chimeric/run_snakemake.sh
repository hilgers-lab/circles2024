# conda activate ciriquant w/ STAR/2.7.7a
module load slurm snakemake R/4.0.0 sambamba

read workdir configyaml params <<< "$@"

snakemake -d ${workdir} ${params} --configfile $configyaml --cluster "SlurmEasy -n {rule} -t {threads} -m 8G -l cluster_log" --latency-wait 180 --keep-going
