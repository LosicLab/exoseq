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
nf-core/ExoSeq : Exome/Targeted sequence capture best practice analysis v${params.version}
===============================================================================

Usage: nextflow nf-core/ExoSeq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 --kitfiles 'kitpath' --metafiles 'metapath'

This is a typical usage where the required parameters (with no defaults) were
given. The available paramaters are listed below based on category

Required parameters:
    --reads                       Absolute path to project directory
    --genome                      Name of iGenomes reference

Output:
    --outdir                      Path where the results to be saved [Default: './results']
    -w/--work-dir                 The temporary directory where intermediate data will be saved

Kit files:
    --kitfiles                    Path to kitfiles defined in metafiles.config
    --metafiles                   Path to metafiles defined in metafiles.config
    --kit                         Kit used to prep samples [Default: 'agilent_v5']

AWSBatch options:
    --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
    --awsregion                   The AWS Region for your AWS Batch job to run on

For more detailed information regarding the parameters and usage refer to package
documentation at https://github.com/nf-core/ExoSeq""".stripIndent()

// Variables and defaults
params.name = false
params.help = false
params.reads = false
params.singleEnd = false
params.run_id = false
//params.aligner = 'bwa' //Default, but stay tuned for later ;-)
params.saveReference = true
params.exome = true


// Output configuration
params.outdir = './results'
params.saveAlignedIntermediates = false
params.saveIntermediateVariants = false


// Clipping options
params.notrim = false
params.saveTrimmed = false
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0


// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Check blocks for certain required parameters, to see they are given and exist
if (!params.reads || !params.genome){
    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
}
if (!params.kitfiles){
    exit 1, "No Exome Capturing Kit specified!"
}
if (!params.metafiles){
    exit 1, "No Exome Metafiles specified!"
}
//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Create a channel for input files

Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }


// Validate Input indices for BWA Mem and GATK
if(params.aligner == 'bwa' ){
    bwaId = Channel
        .fromPath("${params.gfasta}.bwt")
        .ifEmpty { exit 1, "BWA index not found: ${params.gfasta}.bwt" }
}


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
summary['Container']      = workflow.container
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

// Build BWA Index if this is required

if(! params.bwa_index){
    // Create Channels
    fasta_for_bwa_index = Channel
        .fromPath("${params.gfasta}")
    fasta_for_samtools_index = Channel
        .fromPath("${params.gfasta}")
    // Create a BWA index for non-indexed genomes
    process makeBWAIndex {


        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_bwa_index

        output:
        file "*.{amb,ann,bwt,pac,sa}" into bwa_index

        script:
        """
        bwa index $fasta
        """
    }
    // Create a FastA index for non-indexed genomes
    process makeFastaIndex {

        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_samtools_index

        output:
        file "*.fai" into samtools_index

        script:
        """
        samtools faidx $fasta
        """
    }
} else {
    bwa_index = file("${params.bwa_index}")
}

/*
 * 
 * STEP 0 - FastQC 
 * 
*/

process fastqc {

    tag "$name"
        publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

        input:
        set val(name), file(reads) from read_files_fastqc

        output:
        file '*_fastqc.{zip,html}' into fastqc_results
        file '.command.out' into fastqc_stdout

        script:
        """
        fastqc -q $reads
        """
}
/*
 * STEP 1 - trim with trim galore
 */

if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
    trimgalore_logs = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy', 
            saveAs: {filename -> 
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }
        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file(reads) into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results, trimgalore_logs


        script:
        single = reads instanceof Path
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            trim_galore --gzip $c_r1 $tpc_r1 $reads --fastqc
            """
        } else {
            """
            trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads --fastqc
            """
        }
    }
}


/*
 * STEP 2 - Map with BWA Mem
 */

process bwamem {
    tag "$name"

    input:
    set val(name), file(reads) from trimmed_reads
    file(bwa_index) from bwa_index

    output:
    set val(name), file("${name}_bwa.bam") into samples_sorted_bam
    file '.command.log' into bwa_stdout


    script:
    def avail_mem = task.memory ? "-m ${task.memory.toMega().intdiv(task.cpus)}M" : ''
    rg="\'@RG\\tID:${params.run_id}\\tSM:${params.run_id}\\tPL:illumina\'"

    """
    bwa mem \\
    -R $rg \\
    -t ${task.cpus} \\
    $params.gfasta \\
    $reads \\
    | samtools sort ${avail_mem} -O bam -T - >${name}_bwa.bam
    """
}


/*
*  STEP 4 - Mark PCR duplicates in sorted BAM file
*/

process markDuplicates {

    tag "${name}"
    publishDir "${params.outdir}/Picard_Markduplicates/metrics", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".dup_metrics") > 0 ? filename : null }

    input:
    set val(name), file(sorted_bam) from samples_sorted_bam

    output:
    set val(name), file("${name}_markdup.bam"), file("${name}_markdup.bai") into samples_markdup_bam, samples_for_applyBQSR
    file("${name}.dup_metrics") into markdup_results
    file '.command.log' into markDuplicates_stdout

    script:
    """
        mkdir `pwd`/tmp
        java -Xmx${task.memory.toGiga()}g -jar $PICARD MarkDuplicates \\
        INPUT=$sorted_bam \\
        OUTPUT=${name}_markdup.bam \\
        METRICS_FILE=${name}.dup_metrics \\
        REMOVE_DUPLICATES=false \\
        CREATE_INDEX=true \\
        TMP_DIR=`pwd`/tmp
       
    """
}




/*
 * Step 5 - Recalibrate BAM file with known variants and BaseRecalibrator
 *
*/
process recal_bam_files {
    tag "${name}"
    publishDir "${params.outdir}/GATK_Recalibration", mode: 'copy'


    input:
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam

    output:
    set val(name), file("${name}_table.recal") into samples_recal_reports
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
    publishDir "${params.outdir}/GATK_ApplyBQSR", mode: 'copy'

    input:
    set val(name), file("${name}_table.recal") from samples_recal_reports
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_for_applyBQSR

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into bam_vcall, bam_metrics

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

/*
 * Step 7 - Determine quality metrics of mapped BAM files using Picard
 *
*/


process CollectMultiMetrics {
    tag "${name}"
    publishDir "${params.outdir}/picard_multimetrics", mode: 'copy'

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


/*
 * Step 8 - Call Variants with HaplotypeCaller in GVCF mode (differentiate between exome and whole genome data here)
 *
*/

process variantCall {
    tag "${name}"
    publishDir "${params.outdir}/GATK_VariantCalling/", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/variants/, "raw_variants")}

    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(name), file("${name}_variants.vcf"), file("${name}_variants.vcf.idx") into raw_variants

    script:
    if(params.exome){
    """
    gatk HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        -L $params.target \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else { // We have a winner (genome)
    """
    gatk HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    }
}

/*
 * Step 9 - Annotate gVCFs using GATK VariantAnnotator
 * 
*/ 

process annotateVariants {
    tag "${name}"
    publishDir "${params.outdir}/GATK_VariantAnnotation", mode: 'copy', 
 //   saveAs: {filename -> params.saveIntermediateVariants ? "$filename" : null }

    input:
    set val(name), file(raw_vcf), file(raw_vcf_idx) from raw_variants

    output:
    set val(name), file("${name}_annotated_variants.vcf"), file("${name}_annodated_variants.vcf.idx") into ann_raw_variants


// TODO !!
    script:
    """
    
    
    
    """


}



/*
 * Step 10 - Create separate files for SNPs and Indels 
 * 
*/ 

process variantSelect {
    tag "${name}"
    publishDir "${params.outdir}/GATK_VariantSelection", mode: 'copy', 
    saveAs: {filename -> params.saveIntermediateVariants ? "$filename" : null }

    input:
    set val(name), file(ann_raw_vcf), file(ann_raw_vcf_idx) from ann_raw_variants

    output:
    set val(name), file("${name}_snp.vcf"), file("${name}_snp.vcf.idx") into raw_snp
    set val(name), file("${name}_indels.vcf"), file("${name}_indels.vcf.idx") into raw_indels

    script:
    """
    gatk SelectVariants \\
        -R $params.gfasta \\
        --variant $raw_vcf \\
        --output ${name}_snp.vcf \\
        --select-type-to-include SNP \\
        --java-options -Xmx${task.memory.toGiga()}g

    gatk SelectVariants \\
        -R $params.gfasta \\
        --variant $raw_vcf \\
        --output ${name}_indels.vcf \\
        --select-type-to-include INDEL \\
        --select-type-to-include MIXED \\
        --select-type-to-include MNP \\
        --select-type-to-include SYMBOLIC \\
        --select-type-to-include NO_VARIATION \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}



/*
 * Step 11 - Recalibrate SNPs using Omni, 1000G and DBSNP databases 
 * 
*/ 

process recalSNPs {
    tag "${name}"
    publishDir "${params.outdir}/GATK_RecalibrateSNPs/", mode: 'copy', 
    saveAs: {filename -> params.saveIntermediateVariants ? "$filename" : null }

    input:
    set val(name), file(raw_snp), file(raw_snp_idx) from raw_snp

    output:
    set val(name), file("${sample}_filtered_snp.vcf"), file("${sample}_filtered_snp.vcf.idx") into filtered_snp

    script:
    """
    gatk VariantRecalibrator \\
        -R $params.gfasta \\
        --variant $raw_snp \\
        --output ${name}_snp.recal \\
        --max-gaussians 4 \\
        --tranches-file ${name}_snp.tranches \\
        --resource omni,known=false,training=true,truth=true,prior=12.0:$params.omni \\
        --resource 1000G,known=false,training=true,truth=false,prior=10.0:$params.thousandg \\
        --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$params.dbsnp \\
        --mode SNP \\
        -an QD \\
        -an FS \\
        -an MQ \\
        --java-options -Xmx${task.memory.toGiga()}g

    gatk ApplyVQSR \\
        -R $params.gfasta \\
        --output ${name}_filtered_snp.vcf \\
        --variant $raw_snp \\
        --mode SNP \\
        --tranches-file ${name}_snp.tranches \\
        --recal-file ${name}_snp.recal \\
        --truth-sensitivity-filter-level 99.0 \\
        -mode SNP \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}


/*
 * Step 12 - Recalibrate INDELS using the Mills golden dataset 
 * 
*/ 

process recalIndels {
    tag "${name}"
    publishDir "${params.outdir}/GATK_RecalibrateIndels", mode: 'copy', 
    saveAs: {filename -> params.saveIntermediateVariants ? "$filename" : null }

    input:
    set val(name), file(raw_indel), file(raw_indel_idx) from raw_indels

    output:
    set val(name), file("${name}_filtered_indels.vcf"), file("${name}_filtered_indels.vcf.idx") into filtered_indels

    script:
    """
    gatk VariantRecalibrator \\
        -R $params.gfasta \\
        --variant $raw_indel \\
        --output ${name}_indel.recal \\
        --max-gaussians 4 \\
        --tranches-file ${name}_indel.tranches \\
        --resource mills,known=false,training=true,truth=true,prior=12.0:$params.mills \\
        --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$params.dbsnp \\
        -an QD -an DP -an FS -an SOR \\
        -mode INDEL \\
        --java-options -Xmx${task.memory.toGiga()}g

    gatk ApplyVQSR \\
        -R $params.gfasta \\
        --output ${name}_filtered_indels.vcf \\
        --variant $raw_indel \\
        --mode SNP \\
        --tranches-file ${name}_indel.tranches \\
        --recal-file ${name}_indel.recal \\
        --truth-sensitivity-filter-level 99.0 \\
        -mode INDEL \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}

/*
 * Step 13 - Combine recalibrated files again
 * 
*/ 

filtered_snp
    .cross(filtered_indels)
    .map{ it -> [it[0][0], it[0][1], it[0][2], it[1][1], it[1][2]] }
    .set{ variants_filtered }




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
