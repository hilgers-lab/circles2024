# STAR aligned back-splice junction tracks

This workflow uses STAR chimeric to identify back-splice junctions tracks and
produce a format to display it on genome viewers (IGV tested)

# Initialization

Run:
`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
- `<workdir>`: project directory with results files
- `<config yaml>`: configuration containing `genome_gtf` and `genome_fasta`
- `<snakemake parameters>`: snakemake parameters with e.g. number of cores (`-j n`), dryrun (`-n`) etc.

# Output directory

```
<workdir>/
├── FASTQ                      # symlinks to original fastq files
├── STAR_chimeric
│   ├── library                # alignments
│   └── library_IGV_junctions  # junction tracks bed files
└── STAR_index     # STAR index for alignment
```
