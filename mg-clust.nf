#!/usr/bin/env -S nextflow run

// Include modules
include { MODULE1; MODULE1_ASSEMBLY_ONLY; MODULE1_PRECOMPUTED } from './modules/mg-clust-module-1.nf'
include { MODULE2; MODULE2_PRECOMPUTED } from './modules/mg-clust-module-2.nf'
include { MODULE3 } from './modules/mg-clust-module-3.nf'
include { MODULE4 } from './modules/mg-clust-module-4.nf'
include { MODULE5 } from './modules/mg-clust-module-5.nf'
include { MODULE6 } from './modules/mg-clust-module-6.nf'

workflow {

    main:
    if (params.help) {
        log.info """
        MG-Clust: Operational Protein Units from metagenomic paired-end reads

        Usage: nextflow run main.nf [options]

        General:
          --output_dir      DIR   Output directory (default: ${params.output_dir})
          --nslots          INT   CPU threads per tool (default: ${params.nslots})
          --stop_at_module  INT   Stop after module N, 1-6 (default: ${params.stop_at_module})
          --full_output     BOOL  Publish all intermediate outputs (default: ${params.full_output})
          --publish_mode    STR   publishDir mode: copy | symlink | rellink | link | move (default: ${params.publish_mode})
          --skip_tax_annot  BOOL  Skip MODULE4 taxonomic annotation (default: ${params.skip_tax_annot})
          --skip_fun_annot  BOOL  Skip MODULE5 functional annotation (default: ${params.skip_fun_annot})
          --maxForks        INT   Max parallel process instances (default: ${params.maxForks})

        Input mode (selects which TSV samplesheet MODULE1/MODULE2 read):
          --input_mode          STR  from_reads | from_assembly | from_bam | from_orfs (default: ${params.input_mode})
          --from_reads_tsv       TSV  [from_reads]     sample_name, reads1, reads2                      (default: ${params.from_reads_tsv})
          --from_assembly_tsv    TSV  [from_assembly]  sample_name, assembly, reads1, reads2             (default: ${params.from_assembly_tsv})
          --from_bam_tsv         TSV  [from_bam]       sample_name, assembly, bam                        (default: ${params.from_bam_tsv})
          --from_orfs_tsv        TSV  [from_orfs]      sample_name, assembly, bam, orfs_faa, orfs_bed    (default: ${params.from_orfs_tsv})

        MODULE1 — Assembly:
          --assem_preset    STR   MEGAHIT preset (default: ${params.assem_preset})
          --min_contig_len  INT   Minimum contig length in bp (default: ${params.min_contig_len})
          --min_seq         INT   Minimum reads required to assemble (default: ${params.min_seq})

        MODULE2 — ORF prediction:
          --train_file_name STR   Prodigal training file name (default: ${params.train_file_name})

        MODULE3 — ORF clustering:
          --min_orf_len     INT   Minimum ORF length in aa (default: ${params.min_orf_len})
          --clust_thres     NUM   MMseqs2 identity threshold 0-1 (default: ${params.clust_thres})
          --clust_cov_len   NUM   MMseqs2 coverage threshold 0-1 (default: ${params.clust_cov_len})

        MODULE4 — Taxonomic annotation (GTDB):
          --gtdb            PATH  MMseqs2 GTDB database prefix (default: ${params.gtdb})
          --lca_mode        INT   MMseqs2 LCA mode (default: ${params.lca_mode})
          --sensitivity     NUM   MMseqs2 sensitivity (default: ${params.sensitivity})
          --tax_lineage     INT   MMseqs2 --tax-lineage flag (default: ${params.tax_lineage})

        MODULE5 — Functional annotation (KO HMMs):
          --hmm_db          PATH  HMM profiles file (default: ${params.hmm_db})
          --db_mode         STR   ko | ga | evalue scoring strategy (default: ${params.db_mode})
          --ko_list         PATH  KOfam adaptive threshold file, --db_mode ko only (default: ${params.ko_list})
          --evalue_thres    NUM   E-value threshold (default: ${params.evalue_thres})
          --edge_tol        INT   aa tolerance for truncated-ORF rescue, --db_mode ko only (default: ${params.edge_tol})
          --relax_evalue_factor NUM  E-value tightening factor for rescued hits, --db_mode ko only (default: ${params.relax_evalue_factor})
        """.stripIndent()
        exit 0
    }

    stop = params.stop_at_module as Integer
    log.info "Pipeline will run up to module ${stop}"
    if ("${params.skip_tax_annot}" == "true") {
        log.info "Taxonomic annotation will be skipped"
    }
    if ("${params.skip_fun_annot}" == "true") {
        log.info "Functional annotation will be skipped"
    }

    // MODULE1: Assembly + read mapping, precomputed assembly + read mapping, or
    // a fully precomputed assembly + BAM per sample -- selected by --input_mode
    // (always runs)
    if (params.input_mode == "from_reads") {

        reads_ch = channel.fromPath(params.from_reads_tsv, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_name,
                                 file(row.reads1, checkIfExists: true),
                                 file(row.reads2, checkIfExists: true)) }

        n_samples = file(params.from_reads_tsv).readLines().size() - 1
        log.info "Pipeline will run MODULE1 assembly + mapping on ${n_samples} samples from ${params.from_reads_tsv}"
        log.info "Outputs will be written to ${params.output_dir}"
        module1_out = MODULE1(reads_ch)

    } else if (params.input_mode == "from_assembly") {
        if (!params.from_assembly_tsv) {
            error "--from_assembly_tsv is required when --input_mode is 'from_assembly' " +
                  "(TSV columns: sample_name, assembly, reads1, reads2)"
        }
        log.info "Using precomputed assembly + reads input from ${params.from_assembly_tsv} (skipping MODULE1 assembly, still mapping)"
        log.info "Outputs will be written to ${params.output_dir}"

        assembly_reads_ch = channel.fromPath(params.from_assembly_tsv, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_name,
                                 file(row.assembly, checkIfExists: true),
                                 file(row.reads1, checkIfExists: true),
                                 file(row.reads2, checkIfExists: true)) }

        module1_out = MODULE1_ASSEMBLY_ONLY(assembly_reads_ch)
    } else if (params.input_mode == "from_bam") {
        if (!params.from_bam_tsv) {
            error "--from_bam_tsv is required when --input_mode is 'from_bam' " +
                  "(TSV columns: sample_name, assembly, bam)"
        }
        log.info "Using precomputed assembly + bam input from ${params.from_bam_tsv} (skipping MODULE1 assembly/mapping)"
        log.info "Outputs will be written to ${params.output_dir}"

        module1_in_ch = channel.fromPath(params.from_bam_tsv, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_name,
                                 file(row.assembly, checkIfExists: true),
                                 file(row.bam, checkIfExists: true)) }

        module1_out = MODULE1_PRECOMPUTED(module1_in_ch)
    } else if (params.input_mode == "from_orfs") {
        if (!params.from_orfs_tsv) {
            error "--from_orfs_tsv is required when --input_mode is 'from_orfs' " +
                  "(TSV columns: sample_name, assembly, bam, orfs_faa, orfs_bed)"
        }
        log.info "Using precomputed assembly + bam + ORFs input from ${params.from_orfs_tsv} " +
                 "(skipping MODULE1 assembly/mapping and MODULE2 ORF prediction)"
        log.info "Outputs will be written to ${params.output_dir}"
        
        module1_in_ch = channel.fromPath(params.from_orfs_tsv, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_name,
                                 file(row.assembly, checkIfExists: true),
                                 file(row.bam, checkIfExists: true)) }

        module1_out = MODULE1_PRECOMPUTED(module1_in_ch)

        precomputed_orfs_ch = channel.fromPath(params.from_orfs_tsv, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_name,
                                 file(row.orfs_faa, checkIfExists: true),
                                 file(row.orfs_bed, checkIfExists: true)) }
    } else {
        error "Unknown --input_mode '${params.input_mode}'; expected one of: from_reads, from_assembly, from_bam, from_orfs"
    }

    // MODULE2: ORF prediction + coverage estimation (per sample), or -- when
    // --input_mode from_orfs supplied orfs.faa/orfs.bed directly -- staging of
    // those precomputed files plus fresh coverage estimation only
    if (stop >= 2) {
        if (params.input_mode == "from_orfs") {
            module2_in_ch = module1_out
                .map { sample_name, _assembly, bam -> tuple(sample_name, bam) }
                .join(precomputed_orfs_ch, by: 0) // -> tuple(sample_name, bam, orfs_faa, orfs_bed)
            module2_out = MODULE2_PRECOMPUTED(module2_in_ch)
        } else {
            module2_out = MODULE2(module1_out)
        }
    }

    // MODULE3: Concatenate samples, filter ORFs, create MMseqs2 DB, cluster ORFs
    if (stop >= 3) {
        orf_files_ch = module2_out.faa.map { _sn, faa -> faa }.collect()
        module3_out  = MODULE3(orf_files_ch)
    }

    // MODULE4: Taxonomic annotation of contigs against GTDB using MMseqs2
    if (stop >= 4 && !params.skip_tax_annot.toBoolean()) {
        contigs_ch = module1_out
                    .map { sample_name, assembly, _bam -> tuple(sample_name, assembly) }
                    .join(module2_out.bed, by : 0) // joins on sample_name, produces tuple(sample_name, assembly, orf_bed)

        module4_out = MODULE4(contigs_ch)
    }

    // MODULE5: Functional annotation of ORFs against KO HMMs using pyHMMER
    if (stop >= 5 && !params.skip_fun_annot.toBoolean()) {
        faa = module2_out.faa.map { sample_name, faa_path -> tuple(sample_name, faa_path) }
        module5_out = MODULE5(faa)
    }

    // MODULE6: Build unified OPU-ORF coverage table
    if (stop >= 6) {
        meancov_files   = module2_out.meancov.map { _sn, path -> path }.collect()
        readscov_files  = module2_out.readscov.map { _sn, path -> path }.collect()
        tax_annot_files = params.skip_tax_annot.toBoolean() ? Channel.value([]) : module4_out.tax_annot.collect()
        fun_annot_files = params.skip_fun_annot.toBoolean() ? Channel.value([]) : module5_out.fun_annot.collect()
        clust_tsv       = module3_out.clust_tsv
        module6_out     = MODULE6(meancov_files, readscov_files, clust_tsv,
                                tax_annot_files, fun_annot_files)
    }

}
