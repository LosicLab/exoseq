#!/bin/bash

# The purpose of this script is to run the LosicLab's custom fork of nextflow's wholeexome sequencing processing pipeline for fastq->vcf.


# specify account for lsf runs
account='acc_JanssenIBD' # do not leave this to default unless running interactive profile
job_queue='premium' #default queue is alloc

# Run directory. Pipeline outputs to $rundir/bamQC
rundir='/sc/orga/projects/losicb01a/common_folder/nextflow-pipelines/wxs-test'

# Reference genome to use [Default: GRCh38; other refs not yet supported]
ref="GRCh38"

# Path to directory containing the pipeline to run
pipeline='/path/to/pipeline/repo'
mkdir -p $rundir
cd $rundir

inputBAM="$rundir/preprocessing/bwamem-Align/Sample*"

module load nextflow/0.30.2

nextflow run $pipeline/bamQC.nf \
--bam "$inputBAM" \
--genome $ref \
--saveAlignedIntermediates true \
--minerva_account $account \
# do you have multiple lanes per sample that need to be merged?
--multiLane false \ 
--job_queue $job_queue \
# if this is a re-run, pick up where ya left off
-resume \
# local | minerva
-profile local 