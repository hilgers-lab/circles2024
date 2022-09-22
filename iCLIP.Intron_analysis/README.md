# iCLIP intron analysis

The iCLIP.intron_analysis workflow quantifies the enrichment of iCLIP in intron of
circRNAs. Here, first an intron database is constructed, please check the scripts for details.
Then introns are assigned to the feature sets. Finally,
cross-lining profiles and enrichments are computed.

- Cross-link profiles in scaled introns
- Enrichment follows the idea in iCount, Curk et al. (2019) iCount.

# Method details

## Config file variables for filtering

The yaml files allows folllowing filtering:
- `gene_biotypes`: `gene_biotype` as provided in the genome gtf, not if empty
- `intron_maxgap`: <int>
- `exon_boundaries`: <int>
- `circRNA_merge_distance`: <int>
- `intron_size_thrs`: <int>

See `iCLIP.introns_analysis.template.yaml` for further details.

## circRNA preprocessing

circRNA that that overlap in their boundaries within `circRNA_merge_distance`
are collapsed.

## Intron database

The intron database is constructed per gene by subtracting the exonic coordinates from the gene body.
For the construction of the exonic database,
alternative TSS and TES,as well as isoform specific exons are removed.

Filters can be change in `scripts/exon_database.R`.

Minimum Intron sizes are defined by `intron_size_thrs`

## Intron-to-feature assignment

Introns are assigned to their feature while allowing for a maximum offset between
feature and intron boundaries by `intron_max_gap`. This overcomes distances that
might occor at alternative splice/donor sites.

# Initialization

Run:
`bash run_snakemake <workdir> <config yaml> <snakemake parameters>`
- `<workdir>`: project directory with results files
- `<config yaml>`: cmp `iCLIP.introns_analysis.template.yaml`
- `<snakemake parameters>`: snakemake parameters with e.g. number of cores (`-j n`), dryrun (`-n`) etc.

# Output directory

```
<workdir>/
├── feature_db           # constructed feature database
├── iCLIP_centrality     # Disabled. Can be activated by uncommenting in "rule all"
├── iCLIP_coverage       # cross-link count per feature
├── iCLIP_enrichment     # enrichment tables
├── iCLIP_profile        # iCLIP cross-link profiles
├── input                # circRNA; ciri2 format as defined from the circRNA tables
└── Intron_groups        # intron sets by feature (ciri2, dexseq, genes)

```
