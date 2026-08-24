// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2: ORF prediction + coverage estimation
// Input:  assembly FASTA + sorted BAM (per sample)
// Output: ORF protein FASTA + mean coverage + reads coverage tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE2 {

    container "ghcr.io/pereiramemo/mg-clust/module-2:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toBoolean() || params.stop_at_module == 2            
            

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(assembly), path(bam)

    output:
    tuple val(sample_name), path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_orfs.faa"),          emit: faa
    tuple val(sample_name), path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_orfs.bed"),          emit: bed
    tuple val(sample_name), path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_orfs_meancov.tsv"),  emit: meancov
    tuple val(sample_name), path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_orfs_readscov.tsv"), emit: readscov

    script:
    """
    mg-clust-module-2.py \
        --assembly_file  ${assembly} \
        --bam_file       ${bam} \
        --sample_name    ${sample_name} \
        --output_dir     ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --nslots         ${task.cpus} \
        --train_file_name ${params.train_file_name} \
        --overwrite
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2 (precomputed input): stage precomputed ORF protein FASTA + BED files
// instead of running FragGeneScanRs. Reuses the same container/script as MODULE2
// and produces the same tuple(sample_name, faa/bed/meancov/readscov) emit shape,
// so downstream modules cannot tell which entry path ran.
// Input:  sorted BAM + precomputed ORF FASTA + precomputed ORF BED (per sample)
// Output: staged ORF protein FASTA + BED + mean coverage + reads coverage tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE2_PRECOMPUTED {

    container "ghcr.io/pereiramemo/mg-clust/module-2:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toBoolean() || params.stop_at_module == 2

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(orfs_faa), path(orfs_bed)

    output:
    tuple val(sample_name), path("module-2/${sample_name}/${sample_name}_orfs.faa"),          emit: faa
    tuple val(sample_name), path("module-2/${sample_name}/${sample_name}_orfs.bed"),          emit: bed
    tuple val(sample_name), path("module-2/${sample_name}/${sample_name}_orfs_meancov.tsv"),  emit: meancov
    tuple val(sample_name), path("module-2/${sample_name}/${sample_name}_orfs_readscov.tsv"), emit: readscov

    script:
    """
    mg-clust-module-2.py \
        --bam_file             ${bam} \
        --precomputed_orfs_faa ${orfs_faa} \
        --precomputed_orfs_bed ${orfs_bed} \
        --sample_name          ${sample_name} \
        --output_dir           module-2/${sample_name} \
        --nslots               ${task.cpus} \
        --overwrite
    """
}
