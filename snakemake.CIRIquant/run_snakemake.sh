# conda activate ciriquant
set -e
module load slurm R/4.0.3

read workdir circles fastq_dir samplelist_dir params <<< "$@"

[ ! -d "$(realpath $fastq_dir)" ] && { echo "Required input fastq dir. Exit..."; exit; }

conda="--use-conda --conda-prefix <conda prefix>"
cluster=" --latency-wait 180 --keep-going"
config="samplelist_dir=$(realpath $samplelist_dir) fastq_directory=$(realpath $fastq_dir) circles=$(realpath $circles)"
snakemakeconfig="$config $params $conda $cluster"
snakemake --snakefile Snakefile --directory $workdir --configfile "genome.yaml" --config $snakemakeconfig --cluster "SlurmEasy -n {rule} -t {threads} -m 8G -l cluster_log"
