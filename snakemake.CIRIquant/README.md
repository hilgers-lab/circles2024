# CIRIquant snakemake

snakemake'd version for the CIRIquant workflow. 

# Initialization

Run:
`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
- `<workdir>`: project directory with results files
- `<config yaml>`: configuration containing `genome_gtf` and `genome_fasta`
- `<snakemake parameters>`: snakemake parameters with e.g. number of cores (`-j n`), dryrun (`-n`) etc.

# Output directory

```
<workdir>/
CIRIquant/
├── CIRCLES                     # set of reference circRNA
├── CIRIquant                   # Quantification by CIRIquant
│   ├── [...]
│   └── circles_<circRNA set>
│       ├── circ                ## Support files for circRNA quantification
│       ├── gene                ## Support files for gene quantification
│       └── log
├── CIRIquant.dge               # Differential circRNA analysis
│   ├── [...]
│   └── circles_<circRNA set>
│       ├── differential        ## Support files for differential analysis
│       │   └── log
│       └── log
├── CIRIquant.resources         # CIRIquant resources
│   ├── bwa_index
│   └── hisat2_index
├── FASTQ                       # SYMLINKS to fastq files
└── SAMPLESHEETS                # Samplesheets for differential circRNA
```
