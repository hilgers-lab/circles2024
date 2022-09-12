set -e
read name fasta gtf bwa_index hisat2_index <<< "$@"

bwa_cmd="$(which bwa)"
hisat2_cmd="$(which hisat2)"
stringtie_cmd="$(which stringtie)"
samtools_cmd="$(which samtools)"

bwa_index=$(realpath $(dirname $bwa_index))/$(basename $bwa_index)
hisat2_index=$(realpath $(dirname $hisat2_index))/$(basename $hisat2_index)

echo "name: $name"
echo "tools:"
echo -e "  bwa: $bwa_cmd"
echo -e "  hisat2: $hisat2_cmd"
echo -e "  stringtie: $stringtie_cmd"
echo -e "  samtools: $samtools_cmd"
echo "reference:"
echo -e "  fasta: $(realpath $fasta)"
echo -e "  gtf: $(realpath $gtf)"
echo -e "  bwa_index: $bwa_index"
echo -e "  hisat_index: $hisat2_index"
