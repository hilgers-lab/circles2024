## get genomic coordinates out of it - and convert to IGV-junction Files
## Also map the entire library
localrules: chimeric2igv_junctions
rule prep_star:
    input:
        fasta=reference_fasta
    output:
        index=os.path.join("STAR_index","Genome")
    params:
        outdir=os.path.join("STAR_index")
    log: os.path.join("STAR_index","log","STAR_index.log")
    threads: 16
    shell:
        "STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {params.outdir} --genomeFastaFiles {input.fasta} &> {log}"


rule star_chip_mapping_library:
    input:
        R1=os.path.join("FASTQ","{samplename}_R1.fastq.gz"),
        R2=os.path.join("FASTQ","{samplename}_R2.fastq.gz"),
        index=rules.prep_star.output.index,
        genome_gtf=reference_gtf
    output:
        bam=os.path.join("STAR_chimeric","{library}","{samplename}.Aligned.out.bam"),
        chim=os.path.join("STAR_chimeric","{library}","{samplename}.Chimeric.out.junction")
    params:
        index_dir=rules.prep_star.params.outdir,
        outprefix=os.path.join("STAR_chimeric","{library}","{samplename}.")
    log: os.path.join("STAR_chimeric","{library}","log","STAR_{samplename}.log")
    threads: 12
    shell:
        """
        STAR \
        --genomeDir {params.index_dir} \
        --readFilesIn {input.R1} {input.R2} \
        --sjdbGTFfile {input.genome_gtf} --sjdbOverhang 99 \
        --runThreadN {threads} --outReadsUnmapped Fastx --quantMode GeneCounts \
        --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --chimOutType Junctions SeparateSAMold \
        --outFileNamePrefix {params.outprefix} &> {log}

        """

rule bamsort:
    input:
        bam=os.path.join("STAR_chimeric","{library}","{samplename}.Aligned.out.bam")
    output:
        bam=os.path.join("STAR_chimeric","{library}","{samplename}.bam"),
        bai=os.path.join("STAR_chimeric","{library}","{samplename}.bam.bai")
    threads: 8
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && samtools index -@ {threads} {output.bam}"

# requires a rule to convert STAR junctions to IGV junctions
rule chimeric2igv_junctions:
    input:
        junc=os.path.join("STAR_chimeric","{library}","{samplename}.Chimeric.out.junction"),
    output:
        junc=os.path.join("STAR_chimeric","{library}_IGV_junctions","{samplename}.Chimeric.out.junction.bed")
    params:
        span_cutoff = 150000,
        script = os.path.join(maindir, "rules", "scripts", "chimeric2igv_junction_bed.R")
    shell:
        "Rscript {params.script} {params.span_cutoff} {input.junc} > {output.junc}"
