set -e
module load slurm R/4.0.3 bedtools2 snakemake bedops/2.4.39 #slurm

# Example run: bash run_snamekake <config> <workdir> -n -j 12
read yaml directory params <<< "$@"

snakemake --snakefile Snakefile --configfile $yaml -d $directory $params --latency-wait 60 --cluster "SlurmEasy -n {rule} -t {threads} -m 8G -l ./cluster_log"  --keep-going
