# Requirements
# module load slurm snakemake R/4.0.0 sambamba
# conda activate ciriquant w/ STAR/2.7.7a
#
read workdir configyaml params <<< "$@"

snakemake -d ${workdir} ${params} --configfile $configyaml
