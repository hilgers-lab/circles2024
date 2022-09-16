# conda activate ciriquant
set -e

read workdir circles fastq_dir samplelist_dir params <<< "$@"

[ ! -d "$(realpath $fastq_dir)" ] && { echo "Required input fastq dir. Exit..."; exit; }

# Update your <conda prefix>
conda="--use-conda --conda-prefix <conda prefix>"
echo "Update your <conda prefix>"
exit

config="samplelist_dir=$(realpath $samplelist_dir) fastq_directory=$(realpath $fastq_dir) circles=$(realpath $circles)"
snakemakeconfig="$config $params $conda"
snakemake --snakefile Snakefile --directory $workdir --configfile "genome.yaml" --config $snakemakeconfig
