#!/usr/bin/env -S nextflow run

// Include modules
include { MODULE1; MODULE1_PRECOMPUTED } from './modules/mg-clust-module-1.nf'
include { MODULE2 } from './modules/mg-clust-module-2.nf'
include { MODULE3 } from './modules/mg-clust-module-3.nf'
include { MODULE4 } from './modules/mg-clust-module-4.nf'
include { MODULE5 } from './modules/mg-clust-module-5.nf'
include { MODULE6 } from './modules/mg-clust-module-6.nf'

// Derive a sample name from a filename + a single-'*' glob pattern, e.g.
// filename "P01-A01-1.contigs.fa" + pattern "*.contigs.fa" -> "P01-A01-1".
// Used to pair up independently-globbed precomputed assembly/BAM files by sample.
def sampleNameFromGlob(filePath, String pattern) {
    def stars = pattern.count('*')
    if (stars != 1) {
        error "Pattern '${pattern}' must contain exactly one '*' wildcard to derive a sample name"
    }
    def starIdx = pattern.indexOf('*')
    def prefix  = java.util.regex.Pattern.quote(pattern.substring(0, starIdx))
    def suffix  = java.util.regex.Pattern.quote(pattern.substring(starIdx + 1))
    def regex   = java.util.regex.Pattern.compile('^' + prefix + '(.*)' + suffix + '$')
    def m = regex.matcher(filePath.getName())
    if (!m.matches()) {
        error "File '${filePath.getName()}' does not match pattern '${pattern}'"
    }
    return m.group(1)
}

workflow {

    main:
    if (params.help) {
        log.info """
        MG-Clust: Operational Protein Units from metagenomic paired-end reads

        Usage: nextflow run main.nf [options]

        General:
          --input_dir       DIR   Input directory with paired-end FASTQ files (default: ${params.input_dir})
          --reads_pattern   STR   Glob pattern for fromFilePairs (default: ${params.reads_pattern})
          --output_dir      DIR   Output directory (default: ${params.output_dir})
          --nslots          INT   CPU threads per tool (default: ${params.nslots})
          --stop_at_module  INT   Stop after module N, 1-6 (default: ${params.stop_at_module})
          --full_output     BOOL  Publish all intermediate outputs (default: ${params.full_output})
          --skip_tax_annot  BOOL  Skip MODULE4 taxonomic annotation (default: ${params.skip_tax_annot})
          --skip_fun_annot  BOOL  Skip MODULE5 functional annotation (default: ${params.skip_fun_annot})
          --maxForks        INT   Max parallel process instances (default: ${params.maxForks})

        Alternative input (skip MODULE1 assembly + mapping):
          --precomputed_input             BOOL  Use precomputed assembly+bam instead of raw FASTQs (default: ${params.precomputed_input})
          --precomputed_assembly_pattern  STR   Glob pattern for precomputed assemblies under --input_dir (default: ${params.precomputed_assembly_pattern})
          --precomputed_bam_pattern       STR   Glob pattern for precomputed bams under --input_dir (default: ${params.precomputed_bam_pattern})

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
          --hmm_db          PATH  KO HMM profiles file (default: ${params.hmm_db})
          --evalue_thres    NUM   E-value threshold (default: ${params.evalue_thres})
          --cut_ga          BOOL  Use per-profile GA cutoffs (default: ${params.cut_ga})
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

    // MODULE1: Assembly + read mapping, or a precomputed assembly + BAM per
    // sample (whole-run switch, always runs)
    if (params.precomputed_input.toBoolean()) {
        log.info "Using precomputed assembly + bam input (skipping MODULE1 assembly/mapping)"

        assembly_ch = channel.fromPath("${params.input_dir}/${params.precomputed_assembly_pattern}", checkIfExists: true)
            .map { f -> tuple(sampleNameFromGlob(f, params.precomputed_assembly_pattern), f) }

        bam_ch = channel.fromPath("${params.input_dir}/${params.precomputed_bam_pattern}", checkIfExists: true)
            .map { f -> tuple(sampleNameFromGlob(f, params.precomputed_bam_pattern), f) }

        module1_in_ch = assembly_ch.join(bam_ch, by: 0, remainder: true)
            .map { sample_name, assembly, bam ->
                if (assembly == null || bam == null) {
                    error "Sample '${sample_name}' is missing its assembly or BAM -- check that " +
                          "--precomputed_assembly_pattern and --precomputed_bam_pattern match the same set of samples"
                }
                tuple(sample_name, assembly, bam)
            }

        module1_out = MODULE1_PRECOMPUTED(module1_in_ch)
    } else {
        // Read paired-end samples from reads_dir using the pattern defined in nextflow.config
        reads_ch = channel.fromFilePairs(
            "${params.input_dir}/${params.reads_pattern}",
            checkIfExists: true
        )
        module1_out = MODULE1(reads_ch)
    }

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
