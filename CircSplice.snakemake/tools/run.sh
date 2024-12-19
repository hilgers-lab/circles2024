#!/bin/bash

module load snakemake/5.16.0 slurm

snakemake -s /data/hilgers/group2/Shi/2024_circSplice/code/CircSplice_snakefile3.0 --use-conda --cluster "SlurmEasy -t 20" -j 30 --rerun-incomplete --latency-wait 60
