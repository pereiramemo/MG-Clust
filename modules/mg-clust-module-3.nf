// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3: concatenate all samples, filter ORFs, create MMseqs2 DB,
//             and cluster ORFs
// Input:  all per-sample ORF FASTAs (collected)
// Output: filtered ORF MMseqs2 DB + cluster TSV
// ─────────────────────────────────────────────────────────────────────────────

process MODULE3 {

    container "ghcr.io/epereira/mg-clust/${task.process.toLowerCase().replaceFirst('module', 'module-')}:latest"
    publishDir "${params.output_dir}/",
           mode: "copy",
           enabled: params.full_output || params.stop_at_module == 3            


    input:
    path(orf_files)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs_filt_db-minlen${params.min_orf_length}aa/"),         emit: orfs_filt_db
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs_clust-id*/orfs_clust-id*.tsv"), emit: clust_tsv
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs.faa"),   emit: orfs_faa
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs_filt-minlen${params.min_orf_length}aa.faa"),   emit: orfs_filt_faa
    script:
    """
    mg-clust-module-3.py \
        --orf_files      ${orf_files} \
        --output_dir     ${task.process.toLowerCase().replaceFirst('module', 'module-')} \
        --nslots         ${params.nslots} \
        --min_orf_length ${params.min_orf_length} \
        --clust_thres    ${params.clust_thres} \
        --clust_cov_len  ${params.clust_cov_len} \
        --overwrite
    """
}
