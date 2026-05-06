// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3: concatenate all samples, filter ORFs, create MMseqs2 DB,
//             and cluster ORFs
// Input:  all per-sample ORF FASTAs (collected)
// Output: filtered ORF MMseqs2 DB + cluster TSV
// ─────────────────────────────────────────────────────────────────────────────

process MODULE3 {

    container "ghcr.io/epereira/mg-clust/module-3:latest"
    publishDir "${params.output_dir}/intermediate/",
           mode: "copy",
           enabled: params.full_output.toBoolean() || params.stop_at_module == 3            


    input:
    path(orf_files)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs_clust-minlen${params.min_orf_len}aa-id*perc/orfs_clust-minlen${params.min_orf_len}aa-id*perc.tsv"),                      emit: clust_tsv
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs.faa"),                                               emit: orfs_faa
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}/orfs_filt-minlen${params.min_orf_len}aa.faa"),         emit: orfs_filt_faa
    script:
    """
    mg-clust-module-3.py \
        --orf_files      ${orf_files} \
        --output_dir     ${task.process.toLowerCase().replaceFirst('module', 'module-')} \
        --nslots         ${params.nslots} \
        --min_orf_len ${params.min_orf_len} \
        --clust_thres    ${params.clust_thres} \
        --clust_cov_len  ${params.clust_cov_len} \
        --overwrite
    """
}
