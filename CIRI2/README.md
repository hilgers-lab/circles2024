# circRNA detection with CIRI2

A snakemake workflow to identify circRNA using CIRI2. Additionally, the
workflow allows for automated downsampling of FASTQ files by total read counts
(not fractions),in order to check for detection limits (Saturation analysis).

# Execution

`bash run_snakemake <workdir> <config> <snakemake parameters>`


# Setup
The config file, yaml:

```
genome_gtf:
genome_fasta:
```

In order to run to use downsampling, please execute as such:

`bash run_snakemake <workdir> <config> <snakemake parameters> --config downsample=<num>`

`<num>` can either be a) fraction of reads or b) total number of reads. See `seqtk sample` for
more information.

# Output directory

Output of the snakemake instructions

```
rnaser.full_set/
├── CIRI2             # CIRI2 output
├── FASTQ             # symlinks to original fastq files
└── FASTQ_Cutadapt    # trimmed fastq files
```

The `bwa mem` alignments are considered temporary results and thus discarded
after successful execution of the workflow.
