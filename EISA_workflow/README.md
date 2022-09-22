# EISA analysis based on exon-exon/intron junctions

RIP-seq analysis based on exon-exon/intron junction quantification as well as
EISA default approach.

# Execution

`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`

# Initialization

## EEJ/EIJ Quantification

Run:
`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
- `<workdir>`: project directory with results files
- `<config yaml>`: configuration containing `genome_gtf` and `genome_fasta`
- `<snakemake parameters>`: snakemake parameters with e.g. number of cores (`-j n`), dryrun (`-n`) etc.

## Notes

Currently, it is not possible from user input side to provide your own exon or intron database.
This is an interesting feature that could be added, but it is out of scope.

# Output directory

Output of the snakemake instructions

```
<workdir>
├── bamCoverage                     # signal tracks from EIJ
├── eisa_analysis_default.paired    # eisa standard analysis, any read in feature  
├── eisa_analysis_junc.paired       # eisa junction-based analysis
├── featureCounts.classic_paired    # featureCounts for eisa standard analysis
├── featureCounts.eej_paired        # featureCounts; exon-exon juntions
├── featureCounts.eij_paired        # featureCounts; exon-intron junctions
└── feature_db                      # gene, exon, intron databases
```
