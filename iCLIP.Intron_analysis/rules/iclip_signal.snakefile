### iCLIP section
rule feature_introns_iclip:
    input:
        # circRNA.integration.elav_iclip/feature_db/dm6_ensembl96.introns.bed
        introns=os.path.join("feature_db","dm6_ensembl96.introns.bed"),
        iclip=resources['iclip']
    output:
        introns=os.path.join('iCLIP_coverage','dm6_ensembl96.introns.coverage.bed')
    shell:
        "bedtools map -a {input.introns} -b {input.iclip} -c 4 -o sum | awk -v OFS='\\t' '{{ $5=sqrt($7**2); $7=\"\"; print($0) }}' > {output.introns}"

rule flanks_iclip:
    input:
        # circRNA.integration.elav_iclip/feature_db/dm6_ensembl96.introns.bed
        introns=os.path.join("Intron_groups","{featureA}_{type}_{flank}.gff"),
        iclip=resources['iclip']
    output:
        introns=os.path.join('iCLIP_coverage','dm6_ensembl96.introns.{featureA}_{type}_{flank}.coverage.bed')
    log:
        os.path.join("iCLIP_coverage","log","dm6_ensembl96.introns.{featureA}_{type}_{flank}.coverage.log")
    shell:
        "bedtools map -bed -a <(gff2bed < {input.introns} | cut -f 1-6) -b {input.iclip} -c 4 -o sum | awk -v OFS='\\t' '{{ $5=sqrt($7**2); $7=\"\"; print($0) }}' 1> {output.introns} 2> {log}"
