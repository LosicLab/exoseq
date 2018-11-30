#!/usr/bin/env nextflow


/*
========================================================================================
               nf-core/ E X O S E Q    B E S T    P R A C T I C E
========================================================================================
 #### Homepage / Documentation
 https://github.com/scilifelab/NGI-ExoSeq
 #### Authors
 Senthilkumar Panneerselvam @senthil10 <senthilkumar.panneerselvam@scilifelab.se>
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Alex Peltzer @alex_peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
 Marie Gauder <marie.gauder@student.uni-tuebingen.de>

 Some code parts were lent from other NGI-Pipelines (e.g. CAW), specifically the error
 handling, logging messages. Thanks for that @CAW guys.
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - alignment (BWA)
 - recalibration (GATK)
 - realignment (GATK)
 - variant calling (GATK)
 - variant evaluation (SnpEff)
*/

// Help message
helpMessage = """
===============================================================================
LosicLab/ExoSeq - Bam QC: Exome/Targeted sequence capture best practice analysis v${params.version}
===============================================================================

Usage: nextflow </path/to/pipelines/exoseq/bamQC.nf> --bam '/path/to/*_L{1,2,3,4}.bam' --multilane true --genome GRCh38

This is a typical usage where the required parameters (with no defaults) were
given. It is recommended to process one sample per nextflow run.

The available paramaters are listed below based on category

Required parameters:
    --bam                         Absolute path to input bam file for sample (MUST use regex '*' in name)
    --genome                      Name of Genome reference, [Default: 'GRCh38']
    --multiLane                   specify if you have mutliple lanes per sample that need to be concatenated. Assumes 4 lanes [Default: false]

Required if running minerva profile:
    --minerva_account             name of fund for minerva to charge lsf jobs to
    --job_queue                   job queue to submit to [Default: 'alloc']

Output:
    --outdir                      Path where the results to be saved [Default: ./bamQC ]
    -w/--work-dir                 The temporary directory where intermediate data will be saved

Optional Parameters:
    -profile                      load custom profile ['standard' | 'minerva' ; Default: 'standard']
    --name                        custom run name?
    --help                        show help message & exit
    --singleEnd                   input reads are single-end [Default: paired-end]

    --saveAlignedIntermediates    [Default: false]
    --saveIntermediateVariants    [Default: false]

For more detailed information regarding the parameters and usage refer to package
documentation at https://github.com/nf-core/ExoSeq""".stripIndent()


// Output configuration
params.outdir = "./bamQC"
params.saveAlignedIntermediates = false

// Check blocks for certain required parameters, to see they are given and exist
if (!params.bam || !params.genome){
    exit 1, "Parameters '--bam' and '--genome' are required to run the pipeline"
}

// Create a channel for input files
inputBam = Channel
    .fromFilePairs(params.bam, size: params.multiLane ? 4 : 1)
    .ifEmpty { exit 1, "Cannot find any or enough bams matching: ${params.bam}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }


// start pipeline! 


// merge multi-lane bams if necessary
if(! params.multiLane){
    mergedBamForMkDup = inputBam

}else{
    
    process mergeBamFiles {
        tag "${name}"

        publishDir "${params.outdir}/mergedBamFiles", mode: 'symlink',
            saveAs: { filename -> filename.indexOf("*") > 0 ? filename : null }

        input:
        set val(name), file(bams) from inputBam

        output:
        set val(name), file("${name}_merged.sorted.bam"), file("${name}_merged.sorted.bai") into mergedBamResults, mergedBamForMkDup

        script:
        def avail_mem = task.memory ? "${task.memory.toMega().intdiv(task.cpus)}M" : ''

        """
        samtools merge -@ ${task.cpus} ${name}.merged.bam $bams
        samtools sort ${name}.merged.bam -m ${avail_mem} -@ ${task.cpus} -O bam -T ${name} > ${name}.merged.sorted.bam
        samtools index ${name}.merged.sorted.bam
        """
    }

}

// mark duplicate reads with Picard
process markDuplicates {
    tag "${name}"

    publishDir "${params.outdir}/Picard_Markduplicates/metrics", mode: 'symlink',
        saveAs: { filename -> filename.indexOf("*") > 0 ? filename : null }

    input:
    set val(name), file(merged_bam) from mergedBamForMkDup

    output:
    set val(name), file("${name}_merged.sorted.markdup.bam"), file("${name}_merged.sorted.markdup.bai") into samples_markdup_bam, samples_for_applyBQSR, mkdupResults
    file("${name}.dup_metrics") into markdup_results
    file '.command.log' into markDuplicates_stdout

    script:
    """
        java -Xmx${task.memory.toGiga()}g -jar $PICARD MarkDuplicates \\
        INPUT=$merged_bam \\
        OUTPUT=${name}_merged.sorted.markdup.bam \\
        METRICS_FILE=${name}.dup_metrics \\
        REMOVE_DUPLICATES=false \\
        CREATE_INDEX=true \\
        TMP_DIR=`pwd`/tmp
       
    """
}



// Recalibrate BAM file with known variants and BaseRecalibrator

process recalibrateBam {
    tag "${name}"
    publishDir "${params.outdir}/GATK_Recalibration", mode: 'symlink',
        saveAs: { filename -> filename.indexOf("*") > 0 ? filename : null }


    input:
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam

    output:
    set val(name), file("${name}_table.recal") into samples_recal_reports, recalReportResults
    file '.command.log' into gatk_stdout
    file '.command.log' into gatk_base_recalibration_results

    script:
    if(params.exome){
    """
    gatk BaseRecalibrator \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        -L $params.target \\
        --known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else {
    """
    gatk BaseRecalibrator \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        --known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    }
}

process applyBQSR {
    tag "${name}"
    publishDir "${params.outdir}/GATK_ApplyBQSR", mode: 'symlink',
        saveAs: { filename -> filename.indexOf("*") > 0 ? filename : null }

    input:
    set val(name), file("${name}_table.recal") from samples_recal_reports
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_for_applyBQSR

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into bam_vcall, bam_metrics, bam_vanno, recalBamsResult

    script:
    if(params.exome){
    """
    gatk ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        -L $params.target \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else {
    """
    gatk ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    }

}

// calculate QC metrics for bam files with Picard
process collectMultiMetrics {
    tag "${name}"
    publishDir "${params.outdir}/picard_multimetrics", mode: 'symlink',
        saveAs: { filename -> filename.indexOf("*") > 0 ? filename : null }

    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_metrics

    output:
    file "${name}.multimetrics.*" into picard_multimetrics
    file '.command.log' into qc_stdout

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar $PICARD CollectMultipleMetrics \\
      I=$realign_bam \\
      O=${name}.multimetrics \\
      R=$params.gfasta \\
    """
}