#!/bin/bash

# The purpose of this script is to run the LosicLab's custom fork of nextflow's whole exome sequencing processing pipeline for fastq->vcf.

nID='n' # identifier for normal sample
tID='t' # identifier for tumor sample
       # (leave blank & comment out lines using this variable if no somatic mutation calling)

nbam="GATK_ApplyBQSR/tiny*"$nID"*.bam"
tbam="GATK_ApplyBQSR/tiny*"$tID"*.bam"

# Run directory. Pipeline outputs to $rundir/preprocessing
rundir=$PWD"/testdata/tiny"
cd $rundir

# Reference genome to use [Default: GRCh38; other refs not yet supported]
ref="hg38"

# Path to directory containing the pipeline to run
pipeline='/sc/orga/projects/losicb01a/common_folder/nextflow-pipelines/exoseq'

module purge
module load openssl
module load anaconda
module load nextflow
nextflow run $pipeline/variantCalling.nf --nbam "$nbam" --tbam "$tbam" --genome $ref -resume -profile chimera_local
