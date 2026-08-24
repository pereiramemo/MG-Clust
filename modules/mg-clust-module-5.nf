// ─────────────────────────────────────────────────────────────────────────────
// MODULE 5: Functional annotation of ORFs against KO HMMs using pyHMMER
// Input:  ORF protein sequences from module 2 (per sample)
// Output: per-ORF functional annotations
// ─────────────────────────────────────────────────────────────────────────────

process MODULE5 {

    container "ghcr.io/pereiramemo/mg-clust/module-5:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: params.publish_mode,
           enabled: params.full_output.toBoolean() || params.stop_at_module == 5         

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(orfs_faa)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_orfs-minlen${params.min_orf_len}aa-fun_annot.tsv"),       emit: fun_annot
   //  path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name}/${sample_name}_orfs-minlen${params.min_orf_len}aa-fun_annot_all.tsv"), emit: fun_annot_all

    script:
    """
    mg-clust-module-5.py \
        --orfs_faa      ${orfs_faa} \
        --hmm_db        ${params.hmm_db} \
        --db_mode       ${params.db_mode} \
        --ko_list       ${params.ko_list} \
        --sample_name   ${sample_name} \
        --evalue_thres  ${params.evalue_thres} \
        --min_orf_len   ${params.min_orf_len} \
        --edge_tol      ${params.edge_tol} \
        --relax_evalue_factor ${params.relax_evalue_factor} \
        --nslots        ${task.cpus} \
        --output_dir    ${task.process.toLowerCase().replaceFirst('module', 'module-')}/${sample_name} \
        --overwrite
    """
}
