// ─────────────────────────────────────────────────────────────────────────────
// MODULE 6: Merge taxonomy, function, and abundance tables
// Input:  per-sample tax/fun workable TSVs (collected) + module 4 workable tables
// Output: merged orfs_clust_id<N>perc_meancov2taxa2fun_workable.tsv
// ─────────────────────────────────────────────────────────────────────────────

process MODULE6 {

    container "ghcr.io/epereira/mg-clust/module-6:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: "copy",
           enabled: params.full_output || params.stop_at_module == 6                   

    input:
    path(meancov_files)
    path(readscov_files)
    path(clust_tsv)
    path(tax_annot_files)
    path(fun_annot_files)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}")

    script:
    """
    mg-clust-module-6.py \
        --meancov_files  ${meancov_files} \
        --readscov_files ${readscov_files} \
        --clust_tsv      ${clust_tsv} \
        --tax_annot_files ${tax_annot_files} \
        --fun_annot_files ${fun_annot_files} \
        --clust_thres    ${params.clust_thres} \
        --min_orf_len    ${params.min_orf_len} \
        --output_dir     ${task.process.toLowerCase().replaceFirst('module', 'module-')} \
        --overwrite
    """
}