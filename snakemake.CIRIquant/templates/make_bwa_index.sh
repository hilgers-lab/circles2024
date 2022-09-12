# conda activate <environment>

read name genomefa rest <<< "$@"

mkdir -p bwa_index
bwa index -p bwa_index/$name $genomefa 
