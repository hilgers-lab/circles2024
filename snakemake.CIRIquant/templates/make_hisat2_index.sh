# conda activate <environment>
read genomefa rest <<< "$@"
mkdir hisat2_index
hisat2-build -p 24 $genomefa hisat2_index/$(basename $genomefa)
