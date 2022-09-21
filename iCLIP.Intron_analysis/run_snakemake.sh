set -e
# Requirements
# module load R/4.0.3 bedtools2 snakemake bedops/2.4.39

# Example run: bash run_snamekake <config> <workdir> -n -j 12
read directory yaml params <<< "$@"

snakemake --snakefile Snakefile --configfile $yaml -d $directory $params --latency-wait 60 --keep-going
