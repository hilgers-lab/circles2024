
rule feature_databases:
    input:
        gtf=resources["gtf"],
        script_r=os.path.join(maindir, "scripts", "exon_database.R"),
        script_sh=os.path.join(maindir, "scripts", "intron_database.sh")
    output:
        genes=os.path.join("feature_db","dm6_ensembl96.genes.bed"),
        exons=os.path.join("feature_db","dm6_ensembl96.exons.bed"),
        introns=os.path.join("feature_db","dm6_ensembl96.introns.bed")
    params:
        script_r=os.path.join(maindir, "scripts", "exon_database.R"),
        script_sh=os.path.join(maindir, "scripts", "intron_database.sh"),
        biotypes=config["gene_biotypes"],
        intron_size_thrs=config['intron_size_thrs']
    log:
        R=os.path.join("feature_db","log","exon_databases.log"),
        bedtools=os.path.join("feature_db","log","prepare_features.bedtools_introns.log"),
    shell:
        # add intron size filter using e.g. awk
        """
        Rscript {params.script_r} --gtf {input.gtf} --biotypes {params.biotypes} --out_genes {output.genes} --out_exons {output.exons} 2> {log.R} &&
            bash {params.script_sh} {output.introns} {params.intron_size_thrs} {output.genes} {output.exons} &> {log.bedtools}
        """
### Map introns to feature

rule introns_dexseq:
    input:
        exons_set=lambda wildcard: exon_targets[wildcard.featureA],
        genes=rules.feature_databases.output.genes,
        introns=rules.feature_databases.output.introns,
        script=os.path.join(maindir, "scripts", "construct_introns_from_dexseq.R")
    output:
        outupstream=os.path.join("Intron_groups","{featureA}_dexseq_upstream.gff"),
        outdownstream=os.path.join("Intron_groups","{featureA}_dexseq_downstream.gff"),
        outothers=os.path.join("Intron_groups","{featureA}_dexseq_notflanking.gff"),
        outexons=os.path.join("Intron_groups","{featureA}_dexseq.gff")
    params:
        params=" --exon_boundaries {exon_boundaries} --intron_maxgap {intron_maxgap} --outprefix {outpref}".format(
            exon_boundaries=config["exon_boundaries"],
            intron_maxgap=config['intron_maxgap'],
            outpref=os.path.join("Intron_groups","{featureA}"))
    log:
        os.path.join('Intron_groups','log',"introns_exons.{featureA}_dexseq.log")
    shell:
        "Rscript {input.script} {params.params} --dexseq.targets {input.exons_set} -g {input.genes} -i {input.introns} &> {log} "

# create intron set, related to circRNA
rule extract_ciri2_targets:
    input:
        circRNA_tab=lambda wildcards: circRNA_targets[wildcards.featureA],
        ciri2_database=resources['ciri2_reference'],
        script=os.path.join(maindir, "scripts", "extract_ciri2.sh")
    output:
        ciri2_set=os.path.join("input","circRNA.{featureA}.ciri")
    # params:
        # script=os.path.join(maindir, "scripts", "extract_ciri2.sh")
    log: os.path.join("input","log","extract_ciri2_targets.{featureA}.log")
    shell:
        "bash {input.script} {input.circRNA_tab} {input.ciri2_database} 1> {output.ciri2_set} 2> {log}"

rule introns_cir2:
    input:
        ciri2_set=rules.extract_ciri2_targets.output.ciri2_set,
        genes=rules.feature_databases.output.genes,
        exons=rules.feature_databases.output.exons,
        introns=rules.feature_databases.output.introns,
        script=os.path.join(maindir, "scripts", "construct_introns_from_ciri2.R")
    output:
        outupstream=os.path.join("Intron_groups","{featureA}_ciri2_upstream.gff"),
        outdownstream=os.path.join("Intron_groups","{featureA}_ciri2_downstream.gff"),
        outothers=os.path.join("Intron_groups","{featureA}_ciri2_notflanking.gff"),
        # outinternal=os.path.join("Intron_groups","{featureA}_internal.gff"),
        outcircrna=os.path.join("Intron_groups","{featureA}_ciri2.gff")
    params:
        script=os.path.join(maindir, "scripts", "construct_introns_from_ciri2.R"),
        params="--collapse_circRNA {collapse_circRNA} --exon_boundaries {exon_boundaries} --intron_maxgap {intron_maxgap} --outprefix {outpref}".format(
            collapse_circRNA=config["circRNA_merge_distance"],
            exon_boundaries=config["exon_boundaries"],
            intron_maxgap=config['intron_maxgap'],
            outpref=os.path.join("Intron_groups","{featureA}"))
    log: os.path.join("Intron_groups","log","introns_circRNA.{featureA}_ciri2.log")
    shell:
        "Rscript {params.script} {params.params} -c {input.ciri2_set} -g {input.genes} -e {input.exons} -i {input.introns} &> {log}"


rule introns_internal_cir2:
    input:
        ciri2_set=rules.extract_ciri2_targets.output.ciri2_set,
        genes=rules.feature_databases.output.genes,
        exons=rules.feature_databases.output.exons,
        introns=rules.feature_databases.output.introns,
        script=os.path.join(maindir, "scripts", "construct_introns_from_ciri2.R")
    output:
        outupstream=os.path.join("Intron_groups","{featureA}_ciri2_internal_upstream.gff"),
        outdownstream=os.path.join("Intron_groups","{featureA}_ciri2_internal_downstream.gff"),
        outsingle=os.path.join("Intron_groups","{featureA}_ciri2_internal_single.gff"),
        outothers=os.path.join("Intron_groups","{featureA}_ciri2_internal_notflanking.gff"),
        # outinternal=os.path.join("Intron_groups","{featureA}_internal.gff"),
        outcircrna=os.path.join("Intron_groups","{featureA}_internal_ciri2.gff")
    params:
        script=os.path.join(maindir, "scripts", "construct_introns_from_ciri2.R"),
        params="--internal_introns --collapse_circRNA {collapse_circRNA} --exon_boundaries {exon_boundaries} --intron_maxgap {intron_maxgap} --outprefix {outpref}".format(
            collapse_circRNA=config["circRNA_merge_distance"],
            exon_boundaries=config["exon_boundaries"],
            intron_maxgap=config['intron_maxgap'],
            outpref=os.path.join("Intron_groups","{featureA}"))
    log: os.path.join("Intron_groups","log","introns_circRNA.{featureA}_ciri2.log")
    shell:
        "Rscript {params.script} {params.params} -c {input.ciri2_set} -g {input.genes} -e {input.exons} -i {input.introns} &> {log}"

# create intron setfor genes. related to *_{exons,circles}.gff
rule introns_feature_geneset:
    input:
        # use output from rukes introns_*
        geneset=os.path.join("Intron_groups","{group}_{type}.gff"),
        introns=rules.feature_databases.output.introns,
        script=os.path.join(maindir, 'scripts','construct_introns_from_geneset.R')
    output:
        outset=os.path.join("Intron_groups","{group}_{type}_geneset.gff")
    params:
        script=os.path.join(maindir, 'scripts','construct_introns_from_geneset.R'),
        outprefix=os.path.join("Intron_groups","{group}_{type}"),
    log: os.path.join("Intron_groups","log", "introns_geneset.{group}_{type}.log")
    shell:
        "Rscript {params.script} --gene_set {input.geneset} -f gff --introns {input.introns} --outprefix {params.outprefix} &> {log}"

rule introns_genes_geneset:
    input:
        # use output from rukes introns_*
        geneset=lambda wildcard: gene_targets[wildcard.group],
        introns=rules.feature_databases.output.introns,
        script=os.path.join(maindir, 'scripts','construct_introns_from_geneset.R')
    output:
        outset=os.path.join("Intron_groups","{group}_genes_geneset.gff")
    params:
        script=os.path.join(maindir, 'scripts','construct_introns_from_geneset.R'),
        outprefix=os.path.join("Intron_groups","{group}_genes"),
    log: os.path.join("Intron_groups","log", "introns_geneset.{group}_genes.log")
    shell:
        "Rscript {params.script} --gene_set {input.geneset} --introns {input.introns} --outprefix {params.outprefix} &> {log}"
