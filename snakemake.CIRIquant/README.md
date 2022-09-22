# CIRIquant snakemake

snakemake'd version for the CIRIquant workflow.


# Initialization

Run:
`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
- `<workdir>`: project directory with results files
- `<config yaml>`: configuration containing `genome_gtf` and `genome_fasta`
- `<snakemake parameters>`: snakemake parameters with e.g. number of cores (`-j n`), dryrun (`-n`) etc.

# Setup

1. Install conda environment using the `conda.ciriquant.yaml`.
2. Configure ciriQuant, see github doc at https://github.com/bioinfo-biols/CIRIquant
3. Prepare `genome.yaml` with following items:
  - `name`: <project name>
  - `gtf`: genome annotation file, gtf format
  - `fasta`: genome sequences, fasta format
4. place `genome.yaml` in the same folder as the `run_snakemake.sh`

# Output directory

```
<workdir>/
CIRIquant/
├── CIRCLES                     # set of reference circRNA
├── CIRIquant                   # Quantification by CIRIquant
<<<<<<< HEAD
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
=======
│   ├── [...]                   # one per <circRNA set>
│   └── circles_<circRNA set>
│       ├── circ                ## Support files for circRNA quantification
│       └── gene                ## Support files for gene quantification
├── CIRIquant.dge               # Differential circRNA analysis
│   ├── [...]                   # one per <circRNA set>
│   └── circles_<circRNA set>
│       └── differential        ## Support files for differential analysis
>>>>>>> 49995272bc42fc9695b748656ae7e728e2758dd8
├── CIRIquant.resources         # CIRIquant resources
│   ├── bwa_index
│   └── hisat2_index
├── FASTQ                       # SYMLINKS to fastq files
└── SAMPLESHEETS                # Samplesheets for differential circRNA
```
