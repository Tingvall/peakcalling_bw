#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run peakcalling_bw.nf --samples samples.txt --outdir outdir

  """.stripIndent()
  }


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * VALIDATE
 */

if (params.samples)     { ch_samples = Channel.fromPath(params.samples, checkIfExists: true) } else { exit 1, 'Samples not specified' }
    ch_samples
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample, row.count, row.rep, row.type, file(row.bam), file(row.bam_ctrl)] }
        .into { ch_samples_split1; ch_samples_split2}


println ("""
        ===========================================================================================
                                            Macs2 peakcallig & Pooled biwbig generation
        ===========================================================================================
        Sample file: ${params.samples}
        Macs2 q-value: ${params.macs_q}
        Generate consensus peaks: ${params.skip_consensus}
        Generate pooled bigwig: ${params.skip_bigwig}
        ===========================================================================================
        """)



/*
 * 1. Peak calling
 */
process PEAKCALLING {
    publishDir "${params.outdir}/${sample}/peaks", mode: 'copy', pattern: '*_peaks.*'

    input:
    set val(sample), val(rep), val(count), val(type), path(bam), path(bam_ctrl) from ch_samples_split1

    output:
    tuple val(sample), val(count), val(rep), path("${sample}_rep${rep}_peaks.${type}Peak") into ch_peaks


    script:
    """
    # peak calling for replicates
    macs2 callpeak -t ${bam} -c ${bam_ctrl} -f BAM -g ${params.genome_size} -n ${sample}_rep${rep} -B -q ${params.macs_q}  2> ${sample}_rep${rep}_macs2.log
    """
}

ch_peaks.groupTuple()
    .map { it -> [ it[0], it[1].max() as int, it[2].join(' '), it[3].join(' ')] }
    .set { ch_peaks_group}


/*
 * 2. Generate consensus peaks per sample
 */
 process INTERACTION_PEAK_INTERSECT {
   //publishDir "${params.outdir}/${sample}/peaks", mode: 'copy', pattern: '*_consensus_peaks.bed'
   publishDir "${params.outdir}/${sample}/peaks", mode: 'copy', pattern: '*.bed'

   when:
   !params.skip_consensus

   input:
   tuple val(sample), val(count), val(rep), val(peaks) from ch_peaks_group

   output:
   path "${sample}_consensus_peaks.bed" into ch_consensus_peaks


   script:
     """
     #Combine peaks per sample and count number of sample overlapping each region
     bedtools multiinter -i $peaks -names $rep > ${sample}_merge.bed

     #Filter based on region overlapping
     awk '\$4 >= $count' ${sample}_merge.bed > ${sample}_merge_filt.bed

     #Merge regions
     bedtools merge -i ${sample}_merge_filt.bed > ${sample}_consensus_peaks.bed
     """
 }


 ch_samples_split2.groupTuple()
     .map { it -> [ it[0], it[4].join(' '), it[5].join(' ')] }
     .set { ch_bam_group}


 /*
  * 3. Pool bams for merged bigwig generation
  */
process POOL_BAMS {

    when:
    !params.skip_bigwig

    input:
    set val(sample), val(bams), val(bam_ctrls) from ch_bam_group

    output:
    tuple val(sample), path("${sample}_pooled.bam"), path("${sample}_ctrl_pooled.bam") into ch_bam_pooled

    script:
    """
    # Merge bam
    samtools merge -u ${sample}_pooled.bam ${bams}

    # Merge ctrl bams
    samtools merge -u ${sample}_ctrl_pooled.bam ${bam_ctrls}
    """
}



/*
 * 11. Generate pooled bigwigs for samples and ctrls
 */
process GENERATE_POOLED_BIGWIGS {
    publishDir "${params.outdir}/${sample}/bigwigs", mode: 'copy', pattern: '*.bw'

    when:
    !params.skip_bigwig

    input:
    set val(sample), path(bam_pool), path(bam_ctrl_pool) from ch_bam_pooled

    output:
    tuple val(sample), path("${sample}_pooled.bw"), path("${sample}_Input_pooled.bw")  into ch_bigwig_pooled

    script:
    """
    # Generate bigwigs for pooled reps/inputs
    samtools index ${bam_pool}
    bamCoverage -b ${bam_pool} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample}_pooled.bw
    samtools index ${bam_ctrl_pool}
    bamCoverage -b ${bam_ctrl_pool} --effectiveGenomeSize ${params.genome_size} --centerReads --normalizeUsing RPGC -o ${sample}_Input_pooled.bw
    """
}
