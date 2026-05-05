#!/usr/bin/env nextflow

// Include modules
include { MODULE1 } from './modules/mg-clust-module-1.nf'
include { MODULE2 } from './modules/mg-clust-module-2.nf'
include { MODULE3 } from './modules/mg-clust-module-3.nf'
include { MODULE4 } from './modules/mg-clust-module-4.nf'
include { MODULE5 } from './modules/mg-clust-module-5.nf'
include { MODULE6 } from './modules/mg-clust-module-6.nf'

workflow {

    main:
    // Read paired-end samples from reads_dir using the pattern defined in nextflow.config
    reads_ch = channel.fromFilePairs(
        "${params.input_dir}/${params.reads_pattern}",
        checkIfExists: true
    )

    stop = params.stop_at_module as Integer
    log.info "Pipeline will run up to module ${stop}"
    if ("${params.skip_tax_annot}" == "true") {
        log.info "Taxonomic annotation will be skipped"
    }
    if ("${params.skip_fun_annot}" == "true") {
        log.info "Functional annotation will be skipped"
    }

    // MODULE1: Assembly + read mapping (always runs)
    module1_out = MODULE1(reads_ch)

    // MODULE2: ORF prediction + coverage estimation (per sample)
    if (stop >= 2) {
        module2_out = MODULE2(module1_out)
    }

    // MODULE3: Concatenate samples, filter ORFs, create MMseqs2 DB, cluster ORFs
    if (stop >= 3) {
        orf_files_ch = module2_out.faa.map { _sn, faa -> faa }.collect()
        module3_out  = MODULE3(orf_files_ch)
    }

    // MODULE4: Taxonomic annotation of contigs against GTDB using MMseqs2
    if (stop >= 4 && !params.skip_tax_annot) {
        contigs_ch = module1_out
                    .map { sample_name, assembly, _bam -> tuple(sample_name, assembly) }
                    .join(module2_out.bed, by : 0) // joins on sample_name, produces tuple(sample_name, assembly, orf_bed)
            
        module4_out = MODULE4(contigs_ch)
    }

    // MODULE5: Functional annotation of ORFs against KO HMMs using pyHMMER
    if (stop >= 5 && !params.skip_fun_annot) {
        faa = module2_out.faa.map { sample_name, faa_path -> tuple(sample_name, faa_path) }
        module5_out = MODULE5(faa)
    }

    // MODULE6: Build unified OPU-ORF coverage table
    if (stop >= 6) {
        meancov_files  = module2_out.meancov.map { _sn, path -> path }.collect()
        readscov_files = module2_out.readscov.map { _sn, path -> path }.collect()
        tax_annot_files = module4_out.tax_annot.collect()
        fun_annot_files = module5_out.fun_annot.collect()
        clust_tsv      = module3_out.clust_tsv
        module6_out    = MODULE6(meancov_files, readscov_files, clust_tsv,
                                 tax_annot_files, fun_annot_files)
    }

}
