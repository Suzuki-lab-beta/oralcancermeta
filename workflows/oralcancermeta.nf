/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQTK_SEQ as SEQTK_SEQ_CONTIG_FILTER       } from '../modules/nf-core/seqtk/seq/main'
include { BARRNAP                                    } from '../modules/nf-core/barrnap/main'
include { MINIMAP2_INDEX                             } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                             } from '../modules/nf-core/minimap2/align/main'
include { SEMIBIN_SINGLEEASYBIN                      } from '../modules/nf-core/semibin/singleeasybin/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                           } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_oralcancermeta_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ORALCANCERMETA {

    take:
    fastq_gz // [ [ meta ], reads.fastq.gz ]     , reads (mandatory)
    fasta_gz // [ [ meta ], assembly.fasta.gz ]  , assemblies/genomes (mandatory)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    //
    // MODULE: Run Seqtk
    //

    ch_contigs_filtered = SEQTK_SEQ_CONTIG_FILTER ( fasta_gz ).fastx
    ch_versions = ch_versions.mix(SEQTK_SEQ_CONTIG_FILTER.out.versions)


    //
    // MODULE: Filter rRNA with Barrnap
    //

    BARRNAP ( ch_contigs_filtered )
    ch_versions = ch_versions.mix(BARRNAP.out.versions)

    //
    // Module minimap2 index
    //

    MINIMAP2_INDEX ( ch_contigs_filtered )
    ch_index = MINIMAP2_INDEX.out.index

    ch_fastq_contigs_minimap2 = ch_index.combine(fastq_gz, by: 0)
    ch_fastq_contigs_minimap2.unique()

    ch_reads_minimap2 = ch_fastq_contigs_minimap2.map { meta, index, reads -> return [ meta, reads ] }
    ch_index_minimap2 = ch_fastq_contigs_minimap2.map { meta, index, reads -> return [ meta, index ] }

    //
    // MODULE: minimap2 - Map raw reads to contigs and sort them
    //


    MINIMAP2_ALIGN (
            ch_reads_minimap2,
            ch_index_minimap2,
            true,
            false,
            false,
            false
        )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // MODULE: SemiBin2
    //

    //Unzip the 
    ch_semibin_reads = ch_contigs_filtered.join(MINIMAP2_ALIGN.out.bam).view()

    SEMIBIN_SINGLEEASYBIN(ch_semibin_reads)
    ch_versions = ch_versions.mix(SEMIBIN_SINGLEEASYBIN.out.versions)



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
