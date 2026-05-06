// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2: ORF prediction + coverage estimation
// Input:  assembly FASTA + sorted BAM (per sample)
// Output: ORF protein FASTA + mean coverage + reads coverage tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE2 {

    container "ghcr.io/epereira/mg-clust/module-2:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: "copy",
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
        --nslots         ${params.nslots} \
        --train_file_name ${params.train_file_name} \
        --overwrite
    """
}
