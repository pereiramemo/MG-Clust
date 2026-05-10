#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
================================================================================
  mg-clust: metagenomics ORF clustering pipeline
================================================================================
  Modules:
    1 - De novo assembly (MEGAHIT) + read mapping (BWA-MEM + samtools)
    2 - ORF prediction (FragGeneScanRs) + coverage (bedtools)
    3 - Concatenate samples, filter ORFs by length (bbduk), create MMseqs2 DB
    4 - Cluster ORFs (MMseqs2) + generate abundance tables
================================================================================
*/

// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1: de novo assembly + read mapping
// Input:  paired-end reads (R1, R2)
// Output: assembly FASTA + sorted BAM
// ─────────────────────────────────────────────────────────────────────────────

process MODULE1 {
    tag "${sample_name}"
    container "ghcr.io/epereira/mg-clust/module-1:latest"

    publishDir "${params.output_dir}/module-1/${sample_name}", mode: "copy"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name),
          path("assembly/${sample_name}.contigs.fa"),
          path("${sample_name}_sorted.bam")

    script:
    def markdup_flag = params.markdup ? "--markdup" : ""
    """
    mg-clust-module-1.py \
        --reads1             ${reads[0]} \
        --reads2             ${reads[1]} \
        --sample_name        ${sample_name} \
        --output_dir         . \
        --nslots             ${params.nslots} \
        --assem_preset       ${params.assem_preset} \
        --min_contig_length  ${params.min_contig_len} \
        --min_seq            ${params.min_seq} \
        ${markdup_flag} \
        --overwrite
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2: ORF prediction + coverage estimation
// Input:  assembly FASTA + sorted BAM (per sample)
// Output: ORF protein FASTA + mean coverage + reads coverage tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE2 {
    tag "${sample_name}"
    container "ghcr.io/epereira/mg-clust/module-2:latest"

    publishDir "${params.output_dir}/module-2/${sample_name}", mode: "copy"

    input:
    tuple val(sample_name), path(assembly), path(bam)

    output:
    tuple val(sample_name),
          path("${sample_name}_orfs.faa"),
          path("${sample_name}_orfs_meancov.tsv"),
          path("${sample_name}_orfs_readscov.tsv")

    script:
    """
    mg-clust-module-2.py \
        --assembly_file  ${assembly} \
        --bam_file       ${bam} \
        --sample_name    ${sample_name} \
        --output_dir     . \
        --nslots         ${params.nslots} \
        --train_file_name ${params.train_file_name} \
        --overwrite
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// MODULE 4: cluster ORFs + generate abundance tables
// Input:  MMseqs2 DB + merged coverage tables
// Output: per-cluster abundance tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE4 {
    container "ghcr.io/epereira/mg-clust/module-4:latest"

    publishDir "${params.output_dir}/module-4", mode: "copy"

    input:
    path(orfs_db)           // all orfs_filt_db* files staged in the work dir
    path(meancov_table)
    path(readscov_table)

    output:
    path("*")

    script:
    // pass only the DB prefix (without file extensions) to the script
    """
    mg-clust-module-4.py \
        --orfs_db        orfs_filt_db \
        --meancov_table  ${meancov_table} \
        --readscov_table ${readscov_table} \
        --output_dir     . \
        --nslots         ${params.nslots} \
        --clust_thres    ${params.clust_thres} \
        --clust_cov_len  ${params.clust_cov_len} \
        --overwrite
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW
// ─────────────────────────────────────────────────────────────────────────────



workflow {

    // Read paired-end samples from reads_dir using the pattern defined in nextflow.config
    reads_ch = Channel.fromFilePairs(
        "${params.reads_dir}/${params.reads_pattern}",
        checkIfExists: true
    )

    // Modules 1 and 2 run per sample in parallel
    mod1_out = MODULE1(reads_ch)
    mod2_out = MODULE2(mod1_out)

    // Collect per-sample outputs into flat lists before running module 3
    orf_files_ch      = mod2_out.map { _s, faa, _m, _r -> faa }.collect()
    meancov_files_ch  = mod2_out.map { _s, _f, mea, _r -> mea }.collect()
    readscov_files_ch = mod2_out.map { _s, _f, _m, rds -> rds }.collect()

    // Module 3 runs once across all samples
    mod3_out = MODULE3(orf_files_ch, meancov_files_ch, readscov_files_ch)

    // Module 4 runs once using module 3 outputs
    // The MMseqs2 DB consists of multiple files — collect() stages them all together
    MODULE4(
        mod3_out.orfs_db.collect(),
        mod3_out.meancov,
        mod3_out.readscov
    )
}
