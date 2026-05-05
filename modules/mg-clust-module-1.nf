// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1: de novo assembly + read mapping
// Input:  paired-end reads (R1, R2)
// Output: assembly FASTA + sorted BAM
// ─────────────────────────────────────────────────────────────────────────────

process MODULE1 {

    container "ghcr.io/epereira/mg-clust/module-1:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: "copy",
           enabled: params.full_output || params.stop_at_module == 1            

    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name),
          path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/assembly/${sample_name}.contigs.fa"),
          path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_sorted.bam")
    
    script:
    """
    mg-clust-module-1.py \
        --reads1             ${reads[0]} \
        --reads2             ${reads[1]} \
        --sample_name        ${sample_name} \
        --output_dir         ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --nslots             ${params.nslots} \
        --assem_preset       ${params.assem_preset} \
        --min_contig_length  ${params.min_contig_len} \
        --min_seq            ${params.min_seq} \
        --overwrite
    """

}
