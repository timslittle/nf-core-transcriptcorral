/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTranscriptcorral.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check that rRNA databases exist for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}
//TODO: Not sure if this warning is functioning correctly.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_TRIMMED } from '../modules/local/multiqc_tsv_from_list'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                } from '../subworkflows/local/input_check'
include { FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastqc_umitools_trimgalore'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'   
include { SPADES as SPADES_SC         } from '../modules/nf-core/spades/main'
include { SPADES as SPADES_RNA        } from '../modules/nf-core/spades/main'
include { TRINITY                     } from '../modules/nf-core/trinity/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { BUSCO                       } from '../modules/nf-core/busco/main' 
include { HISAT2_BUILD                } from '../modules/nf-core/hisat2/build/main' 
include { HISAT2_ALIGN                } from '../modules/nf-core/hisat2/align/main'
include { TRANSDECODER_LONGORF        } from '../modules/nf-core/transdecoder/longorf/main'
include { TRANSDECODER_PREDICT        } from '../modules/nf-core/transdecoder/predict/main' 
include { HMMER_HMMSEARCH             } from '../modules/nf-core/hmmer/hmmsearch/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

// EvidentialGene/EviGene process

process EVIGENE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cd-hit=4.8.1 bioconda::exonerate=2.4 bioconda::blast=2.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        "https://depot.galaxyproject.org/singularity/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:0ed9db56fd54cfea67041f80bdd8b8fac575112f-0" : 
        "quay.io/biocontainers/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:0ed9db56fd54cfea67041f80bdd8b8fac575112f-0" }"

    input:
    tuple val(meta), path(multiassembly)

    output:
    tuple val(meta), path("okayset/*.okay.aa") , emit: metaassemblyOrfs
    tuple val(meta), path("okayset/*.okay.cds"), emit: metaassembly
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2022.05.07' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    #Perl script found in x

    ${workflow.projectDir}/bin/evigene/evigene/scripts/prot/tr2aacds4.pl \\
        -NCPU=$task.cpus \\
        -MAXMEM=${task.memory.mega} \\
        -logfile \\
        -cdnaseq $multiassembly \\
        -tidy \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            EvidentialGene: $VERSION
            blast: \$( blastn -version | head -n1 | awk '{print \$2}')
            cd-hit: \$( cd-hit -h | head -n1 | cut -f 1 -d "(" | cut -f 2 -d "n" )
            exonerate: \$( exonerate -v | head -n1 | cut -f 5 -d " " )
        END_VERSIONS
    """
}

workflow TRANSCRIPTCORRAL {

    ch_versions = Channel.empty()

    // Option to only perform meta-assembly with evigene and no de novo assembly.
    if (!params.only_evigene) {

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //

        // Read in files and combine into individual channels the samples needing concatenating (based on sample prefix in the first column)
        INPUT_CHECK (
            ch_input
        )
        .reads
        .map {
            meta, fastq ->
                def meta_clone = meta.clone()
                meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
                [ meta_clone, fastq ]
        }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        //
        // MODULE: Concatenate FastQ files from same sample if required
        //
        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

        //
        // SUBWORKFLOW: Read QC, extract UMI and trim adapters
        //
        FASTQC_UMITOOLS_TRIMGALORE (
            ch_cat_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read
        )
        ch_versions = ch_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.versions)

        //
        // Filter channels to get samples that passed minimum trimmed read count
        //
        ch_fail_trimming_multiqc = Channel.empty()
        ch_filtered_reads = FASTQC_UMITOOLS_TRIMGALORE.out.reads
        if (!params.skip_trimming) {
            ch_filtered_reads
                .join(FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
                .map {
                    meta, reads, trim_log ->
                        if (!meta.single_end) {
                            trim_log = trim_log[-1]
                        }
                        num_reads = WorkflowTranscriptcorral.getTrimGaloreReadsAfterFiltering(trim_log)
                        [ meta, reads, num_reads ]
                }
                .set { ch_num_trimmed_reads  }

            ch_num_trimmed_reads
                .map { meta, reads, num_reads -> if (num_reads > params.min_trimmed_reads) [ meta, reads ] }
                .set { ch_filtered_reads }

            ch_num_trimmed_reads
                .map {
                    meta, reads, num_reads ->
                    if (num_reads <= params.min_trimmed_reads) {
                        return [ "$meta.id\t$num_reads" ]
                    }
                }
                .set { ch_num_trimmed_reads }
    // TODO: Is MULTIQC_TSV_FAIL_TRIMMED running at all?
            MULTIQC_TSV_FAIL_TRIMMED (
                ch_num_trimmed_reads.collect(),
                ["Sample", "Reads after trimming"],
                'fail_trimmed_samples'
            )
            .set { ch_fail_trimming_multiqc }
        }

        //
        // MODULE: HISAT2_BUILD of the genomes to remove
        //

        // TODO: Make it so that a named list of genomes could be provided.
        if (params.filter_genome){
            ch_filter_genome = Channel.fromPath(params.filter_genome, checkIfExists: true)

            if(params.filter_genome_index){
                ch_hisatIndex = Channel.fromPath(params.filter_genome_index)
            } else {
                HISAT2_BUILD(
                    ch_filter_genome,
                    [],
                    []
                ).index
                .set { ch_hisatIndex }
            }

            // Want the unmapped reads from ALIGN, in .fastq
            HISAT2_ALIGN(
                ch_filtered_reads,
                ch_hisatIndex,
                []
            )
            .fastq
            .set { ch_filtered_reads }
        }

        //
        // MODULE: Sortmerna: Remove ribosomal RNA reads
        //
        ch_sortmerna_multiqc = Channel.empty()
        if (params.remove_ribo_rna) {
            ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

            SORTMERNA (
                ch_filtered_reads,
                ch_sortmerna_fastas
            )
            .reads
            .set { ch_filtered_reads }

            ch_sortmerna_multiqc = SORTMERNA.out.log
            ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
        }

        // Create empty assemblies channel which new assemblies will be added to.
        ch_assembly = Channel.empty()

        //
        // MODULE: Spades_SC
        //

        if(params.assemble_spades_sc){
            // TODO: Need to add elements for the 'pacbio' and 'nanopore' inputs in the tuple.
            ch_spades_input=ch_filtered_reads
                .map { [ it[0], it[1], [], [] ] }

            SPADES_SC (
                ch_spades_input,
                [],
                []
            )
            ch_versions = ch_versions.mix(SPADES_SC.out.versions)
            ch_assembly = ch_assembly.mix(SPADES_SC.out.scaffolds)
        }

        //
        // MODULE: Spades_RNA
        //

        if(params.assemble_spades_rna){
            // TODO: Need to add elements for the 'pacbio' and 'nanopore' inputs in the tuple.
            ch_spades_input=ch_filtered_reads
                .map { [ it[0], it[1], [], [] ] }

            SPADES_RNA (
                ch_spades_input,
                [],
                []
            )
            ch_versions = ch_versions.mix(SPADES_RNA.out.versions)
            ch_assembly = ch_assembly.mix(SPADES_RNA.out.scaffolds)
        }

        //
        // MODULE: Trinity
        //

        if(params.assemble_trinity){
            TRINITY(
                ch_filtered_reads
            )
            ch_versions = ch_versions.mix(TRINITY.out.versions)
            ch_assembly = ch_assembly.mix(TRINITY.out.transcript_fasta)
        }

        //
        // Use collect to ensure that downstream processes wait for all assemblies to be done.
        //
    // TODO: Save the combined assembly file.
        ch_assembly
            .map{ it[1] } // To get the assembly files and not meta
            .collectFile(name: "combined_assemblies.fa.gz", 
                newLine: false, skip: 0,
                storeDir: params.outdir)
                
    } else {
        // Need to provide assembly for meta-assembly as a parameter
        ch_assembly = Channel.fromPath(params.input, checkIfExists: true)
            .map { //Defining meta ID as the file name without the last element in '_' TODO: This only makes sense for paired end files not transcript files.
                def meta = [:]
                meta.id = it.getFileName().toString().split('_')[0..-2].join('_')
                [meta, it] 
                }
    }

    // Process the assemblies. Can use EvidentialGene to filter for best transcripts, find ORFs and translate. Or use Transdecoder to find and translate ORFs.

    if(params.use_evigene || params.only_evigene){

        //
        // PROCESS: EVIGENE
        //

        EVIGENE (
            ch_assembly
        )

        ch_assemblyOrfs = EVIGENE.out.metaassemblyOrfs

    } else {

        //
        // MODULE: Transdecoder - ORF detection
        //

        TRANSDECODER_LONGORF(
            ch_assembly
        )

        ch_assemblyOrfs = TRANSDECODER_PREDICT(
            ch_assembly,
            TRANSDECODER_LONGORF.out.folder
        )
        .pep

        ch_versions = ch_versions.mix(TRANSDECODER_LONGORF.out.versions)
    }

    //
    // MODULE: BUSCO
    // 
    BUSCO (
        ch_assemblyOrfs,
        params.busco_lineage,
        params.busco_lineages_path,
        params.busco_config_file
    )
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    //
    // MODULE: HMMer/hmmsearch
    //

    if(params.hmmsearch_hmmfile){

        ch_hmmsearchInput = ch_assemblyOrfs
            .map{ [it[0], 
                params.hmmsearch_hmmfile,
                it[1],
                [],
                [],
                [] ] }

        HMMER_HMMSEARCH (
            ch_hmmsearchInput
        )

        ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    }

    //
    // MODULE: Pipeline reporting
    //
// TODO: Find out why CUSTOM_DUMP is going forever and then uncomment its sections
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowTranscriptcorral.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowTranscriptcorral.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    if (!params.only_evigene) {
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
