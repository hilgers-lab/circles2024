# EISA analysis based on exon-exon/intron junctions

RIP-seq analysis based on exon-exon/intron junction quantification

# TODO

Create eisa.analysis.R CLI

# Execution

Step 1. `bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
Step 2. `Rscript eisa.analysis.R <parameters>`

# Initialization

## EEJ/EIJ Quantification

Run:
`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
- `<workdir>`: project directory with results files
- `<config yaml>`: configuration containing `genome_gtf` and `genome_fasta`
- `<snakemake parameters>`: snakemake parameters with e.g. number of cores (`-j n`), dryrun (`-n`) etc.

# Output directory

Output of the snakemake instructions

```
<workdir>/
├── bamCoverage         # Signal tracks
│   └── log
├── cluster_log
├── featureCounts       # featureCounts on exon-exon junctions
│   └── log
├── featureCounts.eij   # featureCounts on exon-intron junctions
│   └── log
└── tmp
```
