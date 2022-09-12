set -e
module load R/4.0.3 bedtools2/2.27.0

read introns_out_bed intron_size_thrs genes exons <<< "$@"

while read line;
do 
  gid=$(echo -e $line | awk '{print($4)}')
  bedtools subtract -s -a <(echo -e "$line") -b <(grep "$gid" "$exons")
done < <(cat $genes) | bedtools sort | bedtools merge -s -c 4,5,6 -o distinct 1> ${introns_out_bed%.bed}.unfiltered.bed

cat ${introns_out_bed%.bed}.unfiltered.bed | awk -v intron_size_thrs=$intron_size_thrs '$3 -$2 >= intron_size_thrs' > $introns_out_bed
