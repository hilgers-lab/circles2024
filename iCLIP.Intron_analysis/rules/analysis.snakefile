rule enrichment:
    input:
        iclip_bedgraph = resources['iclip'],
        intron_db=rules.feature_databases.output.introns,
        intron_set=lambda wildcard: os.path.join('Intron_groups', '{group}_{type}_{set}.gff'),
        script=os.path.join(maindir, "scripts","iCLIP_enrichment.R")
    output:
        outtable=os.path.join('iCLIP_enrichment','{group}.{type}_{set}.{pseudocnt}.iclip_enrichment.tsv')
    params:
        param='--pseudo_count' if "{pseudocnt}" == 'pseudocount' else ''
    log: os.path.join('iCLIP_enrichment','log','{group}.{type}_{set}.{pseudocnt}.iclip_enrichment.log')
    shell:
        "Rscript {params} {input.script} -b {input.iclip_bedgraph} -g {input.intron_set} -i {input.intron_db} -o {output.outtable} &> {log}"

# change centrality script from aggregation to one per feature
rule centrality:
    input:
        iclip_bedgraph = resources['iclip'],
        # neuronal-down_ciri2_geneset.gff
        intron_set=lambda wildcard: os.path.join('Intron_groups', '{group}_{type}_{set}.gff'),
        script=os.path.join(maindir, "scripts","iCLIP_centrality.R")
    output:
        # group: feature_group, set {up/dowstrean, not_flanking,geneset}, type {dexseq,ciri2}
        outfile=os.path.join('iCLIP_centrality','{group}.{type}_{set}.centrality.tsv')
    params:
    log: os.path.join('iCLIP_centrality','log','{group}.{type}_{set}.centrality.log')
    shell:
        "Rscript {input.script} -b {input.iclip_bedgraph} -g {input.intron_set} -o {output.outfile} &> {log}"

rule iclip_profile:
    input:
        iclip_bedgraph = resources['iclip'],
        # neuronal-down_ciri2_geneset.gff
        intron_set=lambda wildcard: os.path.join('Intron_groups', '{group}_{type}_{set}.gff'),
        script=os.path.join(maindir, "scripts","iCLIP_profile.R")
    output:
        # group: feature_group, set {up/dowstrean, not_flanking,geneset}, type {dexseq,ciri2}
        outfile=os.path.join('iCLIP_profile','{group}.{type}_{set}.profile.tsv')
    params:
    log: os.path.join('iCLIP_profile','log','{group}.{type}_{set}.profile.log')
    shell:
        "Rscript {input.script} -b {input.iclip_bedgraph} -g {input.intron_set} -o {output.outfile} &> {log}"
