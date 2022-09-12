set_ciri2_variables = ['upstream', 'internal_upstream','downstream', 'internal_downstream','internal_single','notflanking','internal_notflanking','geneset']

set_dexseq_variables = ['upstream','downstream','notflanking','geneset']

rule iCLIP_report_ciri2:
    input:
        expand(os.path.join('iCLIP_centrality','{group}.ciri2_{set}.centrality.tsv'),
            group=circ_set, set = set_ciri2_variables),
        expand(os.path.join('iCLIP_enrichment','{group}.ciri2_{set}.default.iclip_enrichment.tsv'),
            group=circ_set, set = set_ciri2_variables),
        expand(os.path.join('iCLIP_profile','{group}.ciri2_{set}.profile.tsv'),
            group=circ_set, set = set_ciri2_variables)
    output:
        os.path.join('iCLIP.introns_analysis.ciri2.html')
    params:
        centrality_set=','.join(expand(os.path.join('iCLIP_centrality','{group}.ciri2_{set}.centrality.tsv'),
            group=circ_set, set = set_ciri2_variables)),
        enrichment_set=','.join(expand(os.path.join('iCLIP_enrichment','{group}.ciri2_{set}.default.iclip_enrichment.tsv'),
            group=circ_set, set = set_ciri2_variables)),
        profile_set=','.join(expand(os.path.join('iCLIP_profile','{group}.ciri2_{set}.profile.tsv'),
            group=circ_set, set = set_ciri2_variables))
    script:
        os.path.join(maindir,"scripts/iCLIP.introns_analysis.ciri2.Rmd")

rule iCLIP_report_dexseq:
    input:
        expand(os.path.join('iCLIP_centrality','{group}.dexseq_{set}.centrality.tsv'),
            group=exon_set, set = set_dexseq_variables),
        expand(os.path.join('iCLIP_enrichment','{group}.dexseq_{set}.default.iclip_enrichment.tsv'),
            group=exon_set, set = set_dexseq_variables),
        expand(os.path.join('iCLIP_profile','{group}.dexseq_{set}.profile.tsv'),
            group=exon_set, set = set_dexseq_variables),
    output:
        os.path.join('iCLIP.introns_analysis.dexseq.html')
    params:
        centrality_set=','.join(expand(os.path.join('iCLIP_centrality','{group}.dexseq_{set}.centrality.tsv'),
            group=exon_set, set = set_dexseq_variables)),
        enrichment_set=','.join(expand(os.path.join('iCLIP_enrichment','{group}.dexseq_{set}.default.iclip_enrichment.tsv'),
            group=exon_set, set = set_dexseq_variables)),
        profile_set=','.join(expand(os.path.join('iCLIP_profile','{group}.dexseq_{set}.profile.tsv'),
            group=exon_set, set = set_dexseq_variables))
    script:
        os.path.join(maindir,"scripts/iCLIP.introns_analysis.dexseq.Rmd")

rule iCLIP_report_genesets:
    input:
        expand(os.path.join('iCLIP_centrality','{group}.genes_{set}.centrality.tsv'),
            group=gene_set, set = ['geneset']),
        expand(os.path.join('iCLIP_enrichment','{group}.genes_{set}.default.iclip_enrichment.tsv'),
            group=gene_set, set = ['geneset']),
        expand(os.path.join('iCLIP_profile','{group}.genes_{set}.profile.tsv'),
            group=gene_set, set = ['geneset'])
    output:
        os.path.join('iCLIP.introns_analysis.genes.html')
    params:
        centrality_set=','.join(expand(os.path.join('iCLIP_centrality','{group}.genes_{set}.centrality.tsv'),
            group=gene_set, set = ['geneset'])),
        enrichment_set=','.join(expand(os.path.join('iCLIP_enrichment','{group}.genes_{set}.default.iclip_enrichment.tsv'),
            group=gene_set, set = ['geneset'])),
        profile_set=','.join(expand(os.path.join('iCLIP_profile','{group}.genes_{set}.profile.tsv'),
            group=gene_set, set = ['geneset']))
    script:
        os.path.join(maindir,"scripts/iCLIP.introns_analysis.genesets.Rmd")
