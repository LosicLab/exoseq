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


// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

if (!params.refs){
    exit 1, "No Exome Metafiles specified!"
}

// Validate Input indices for BWA Mem and GATK
if(params.aligner == 'bwa' ){
    bwaId = Channel
        .fromPath("${params.gfasta}.bwt")
        .ifEmpty { exit 1, "BWA index not found: ${params.gfasta}.bwt" }
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Create a channel for input files
inputBam = Channel
    .fromFilePairs(params.bam, size: params.multiLane ? 4 : 1)
    .ifEmpty { exit 1, "Cannot find any or enough bams matching: ${params.bam}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }

// Create a summary for the logfile
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome Assembly']       = params.genome
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
if(params.aligner == 'bwa'){
    summary['Aligner'] = "BWA Mem"
    if(params.bwa_index)          summary['BWA Index']   = params.bwa_index
    else if(params.gfasta)          summary['Fasta Ref']    = params.gfasta
}
summary['Save Intermediate Aligned Files'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Save Intermediate Variant Files'] = params.saveIntermediateVariants ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

try {
    if( ! workflow.nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
        }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}




// start pipeline! 


// merge multi-lane bams if necessary
if(! params.multiLane){
    mergedBamForMkDup = inputBam

}else{
    
    process mergeBamFiles {
        tag "${name}"

        input:
        set val(name), file(bams) from inputBam

        output:
        set val(name), file("${name}.bam"), file("${name}.bam.bai") into mergedBamResults, mergedBamForMkDup

        script:
        def avail_mem = task.memory ? "${task.memory.toMega().intdiv(task.cpus)}M" : ''

        """
        samtools merge -@ ${task.cpus} ${name}.bam $bams 
        samtools index ${name}.bam
        """
    }

}

process editBamHeaders {
    tag "${name}"
 

}


// mark duplicate reads with Picard
process markDuplicates {
    tag "${name}"

    publishDir "${params.outdir}/${name}/Picard_Markduplicates", mode: 'copy'

    input:
    set val(name), file(merged_bam), file(merged_bai) from mergedBamForMkDup

    output:
    set val(name), file("${name}.sorted.markdup.bam"), file("${name}.sorted.markdup.bai") into unrecal_bam_metrics, samples_markdup_bam, samples_for_applyBQSR, mkdupResults
    file("${name}.dup_metrics") into markdup_results
    file '.command.log' into markDuplicates_stdout

    script:
    """
        java -Xmx${task.memory.toGiga()}g -jar $PICARD MarkDuplicates \\
        INPUT=$merged_bam \\
        OUTPUT=${name}.sorted.markdup.bam \\
        METRICS_FILE=${name}.dup_metrics \\
        REMOVE_DUPLICATES=false \\
        CREATE_INDEX=true \\

       
    """
}

process collectMultiMetrics_unrecal {
    tag "${name}"
    publishDir "${params.outdir}/${name}/picard_multimetrics/", mode: 'copy'

    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from unrecal_bam_metrics

    output:
    file "${name}.multimetrics.*" into picard_multimetrics_unrecal
    file '.command.log' into qc_stdout

    script:
    """
    java -Xmx${task.memory.toGiga()}g -jar $PICARD CollectMultipleMetrics \\
      I=$realign_bam \\
      O=${name}.multimetrics \\
      R=$params.gfasta \\
    """
}

// Recalibrate BAM file with known variants and BaseRecalibrator
/*
process recalibrateBam {
    tag "${name}"
    publishDir "${params.outdir}/GATK_Recalibration", mode: 'symlink',
        saveAs: { filename -> filename.indexOf("*") > 0 ? filename : null }


    input:
    set val(name), file(markdup_bam), file(markdup_bai) from samples_markdup_bam

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
    set val(name), file("${name}.sorted.markdup.recal.bam"), file("${name}.sorted.markdup.recal.bai") into recal_bam_metrics, recalBamsResult

    script:
    if(params.exome){
    """
    gatk ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.sorted.markdup.recal.bam \\
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
        -O ${name}.sorted.markdup.recal.bam \\
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
    set val(name), file(realign_bam), file(realign_bam_ind) from recal_bam_metrics

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
*/


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def exoMessage() {
  // Display ExoSeq message
  log.info "nf-core/ExoSeq ANALYSIS WORKFLOW ~ ${params.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outdir
  log.info "Genome      : " + params.gfasta
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def versionMessage() {
  // Display version message
  log.info "nf-core/ExoSeq ANALYSIS WORKFLOW"
  log.info "  version   : " + version
  log.info workflow.commitId ? "Git info    : ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]" : "  revision  : " + this.grabRevision()
}


workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.exoMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.exoMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
