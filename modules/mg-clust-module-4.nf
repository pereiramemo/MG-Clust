// ─────────────────────────────────────────────────────────────────────────────
// MODULE 4: Taxonomic annotation of contigs against GTDB using MMseqs2
// Input:  contig sequences from module 1 (per sample)
// Output: per-contig taxonomy annotations
// ─────────────────────────────────────────────────────────────────────────────

process MODULE4 {
    
    container "ghcr.io/pereiramemo/mg-clust/module-4:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toBoolean() || params.stop_at_module == 4           


    tag "${sample_name}"

    input:
    tuple val(sample_name), path(assembly), path(orf_bed)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_contig_tax_annot.tsv"),    emit: tax_annot

    script:
    """
    mg-clust-module-4.py \
        --contigs       ${assembly} \
        --bed_file      ${orf_bed} \
        --gtdb          ${params.gtdb} \
        --sample_name   ${sample_name} \
        --lca_mode      ${params.lca_mode} \
        --sensitivity   ${params.sensitivity} \
        --tax_lineage   ${params.tax_lineage} \
        --nslots        ${task.cpus} \
        --output_dir    ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}
