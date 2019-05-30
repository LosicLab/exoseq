#!/bin/bash

# Run directory. Pipeline outputs to $rundir/preprocessing
rundir=$PWD"/testdata/tiny"

# Reference genome to use [Default: GRCh38; other refs not yet supported]
ref="hg38"

# Path to directory containing the pipeline to run
pipeline='/sc/orga/projects/losicb01a/common_folder/nextflow-pipelines/sandbox/exoseq'
#mkdir -p $rundir
cd $rundir

module purge
module load openssl
module load anaconda
module load nextflow/0.30.2

inputBAM="bwamem-Align/tiny*.bam"

module load nextflow/0.30.2

nextflow run $pipeline/bamQC.nf --outdir $rundir --bam "$inputBAM" --genome $ref --saveAlignedIntermediates "true" --multiLane "false" -resume -profile "chimera_local"
