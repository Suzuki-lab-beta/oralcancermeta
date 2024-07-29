/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQTK_SEQ as SEQTK_SEQ_CONTIG_FILTER                 } from '../modules/nf-core/seqtk/seq/main'
include { BARRNAP                                              } from '../modules/nf-core/barrnap/main'
include { MINIMAP2_INDEX                                       } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                                       } from '../modules/nf-core/minimap2/align/main'
include { SEMIBIN_SINGLEEASYBIN                                } from '../modules/nf-core/semibin/singleeasybin/main'
include { SEQKIT_GREP as SEQKIT_GREP_NON_BACTERIAL             } from '../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_NON_PLASMID               } from '../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_PLASMID                   } from '../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_BACTERIAL                 } from '../modules/nf-core/seqkit/grep/main'
include { GENOMAD_DOWNLOAD                                     } from '../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND as GENOMAD_CONSERVATIVE             } from '../modules/nf-core/genomad/endtoend/main'
include { GENOMAD_ENDTOEND as GENOMAD_DEFAULT                  } from '../modules/nf-core/genomad/endtoend/main'
include { MINIMAP2_INDEX   as MINIMAP2_FUNGI_EUKARYOTIC_VIRUS  } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN   as MINIMAP2_ALIGN_EUKARYOTIC_VIRUS  } from '../modules/nf-core/minimap2/align/main'
include { MSAMTOOLS_FILTER                                     } from '../modules/local/msamtools_filter'
include { MULTIQC                                              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                     } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                               } from '../subworkflows/local/utils_nfcore_oralcancermeta_pipeline'

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

    //Double check here join or combine??
    ch_semibin_reads = ch_contigs_filtered.join(MINIMAP2_ALIGN.out.bam)


    SEMIBIN_SINGLEEASYBIN(ch_semibin_reads)
    ch_versions = ch_versions.mix(SEMIBIN_SINGLEEASYBIN.out.versions)


    //
    // MODULE: Seqkit grep to filter the bacterial contigs
    //

    ch_barrnap = BARRNAP.out.gff
    ch_semibin_fasta = SEMIBIN_SINGLEEASYBIN.out.output_fasta

    ch_gff = BARRNAP.out.gff
                    .splitCsv (skip: 1, sep:"\t" )
                    .map { row -> [row[0].id, row[1][0]] }
                    .unique()


    ch_gff.collectFile(newLine: true) { meta, contigs -> [ "${meta}_contigs.txt", contigs ] }
    .map { contigs -> [ ["id": contigs.simpleName.split("_contigs")[0]], contigs ] }
    .set {ch_contigs}


    ch_barrnap_semibin = ch_contigs.combine(ch_semibin_fasta, by: 0 )
    ch_barrnap_semibin


    ch_filter_gff = ch_barrnap_semibin.map { meta, contigs, semibin_reads -> return [ meta, contigs ] }
    ch_filter_semibin = ch_barrnap_semibin.map { meta, contigs, semibin_reads -> return [ meta, semibin_reads ] }


    SEQKIT_GREP_NON_BACTERIAL ( ch_filter_semibin, ch_filter_gff)
    ch_bacterial_contigs = SEQKIT_GREP_BACTERIAL ( ch_filter_semibin, ch_filter_gff)

    ch_versions = ch_versions.mix( SEQKIT_GREP_BACTERIAL.out.versions )

    //
    // Module: geNomad download
    //

    ch_genomad_db = GENOMAD_DOWNLOAD( ).genomad_db
    ch_versions = ch_versions.mix( GENOMAD_DOWNLOAD.out.versions )

    //
    // Module: geNomad predict-conservative and default
    //

    ch_non_bacterial_conservative = GENOMAD_CONSERVATIVE ( SEQKIT_GREP_NON_BACTERIAL.out.filter, ch_genomad_db )
    //SEQKIT_GREP_NON_BACTERIAL.out.filter
    //ch_non_bacterial_default = GENOMAD_DEFAULT  ( SEQKIT_GREP_NON_BACTERIAL.out.filter, ch_genomad_db )
    ch_virus_summaries_tsv_conservative = GENOMAD_CONSERVATIVE.out.virus_summary.map { meta,virus_summary -> [["conservative_virus":meta], virus_summary.countLines()-1]}
    //ch_virus_summaries_tsv_default = GENOMAD_DEFAULT.out.virus_summary.map { meta,virus_summary -> [ ["default_virus":meta], virus_summary.countLines()-1]}
    ch_plasmid_summaries_tsv_conservative = GENOMAD_CONSERVATIVE.out.plasmid_summary.map { meta,virus_summary -> [["conservative_plasmid":meta], virus_summary.countLines()-1]}
    //ch_plasmid_summaries_tsv_default = GENOMAD_DEFAULT.out.plasmid_summary.map { meta,virus_summary -> [ ["default_plasmid":meta], virus_summary.countLines()-1]}
    ch_versions = ch_versions.mix(GENOMAD_CONSERVATIVE.out.versions)

    // Module: seqkit grep, filter out the viral and plasmid contigs from the non-bacterial

    ch_viral_plasmid = GENOMAD_CONSERVATIVE.out.virus_summary.concat(GENOMAD_CONSERVATIVE.out.plasmid_summary)
                    .splitCsv (skip: 1, sep:"\t" )
                    .map { row -> [row[0].id, row[1][0]] }
                    .unique()

   ch_viral_plasmid.collectFile(newLine: true) { meta, contigs -> [ "${meta}_viral_plasmid_contigs.txt", contigs ] }
   .map { contigs -> [ ["id": contigs.simpleName.split("_viral_plasmid_contigs")[0]], contigs ] }
   .set {ch_viral_plasmid_contigs}

   ch_non_viral_plasmid = ch_viral_plasmid_contigs.combine(SEQKIT_GREP_NON_BACTERIAL.out.filter, by: 0 )

   ch_filter_plasmid = ch_non_viral_plasmid.map { meta, contigs, non_bacterial_reads -> return [ meta, contigs ] }
   ch_filter_non_bacterial = ch_non_viral_plasmid.map { meta, contigs, non_bacterial_reads -> return [ meta, non_bacterial_reads ] }

   SEQKIT_GREP_NON_PLASMID ( ch_filter_non_bacterial, ch_filter_plasmid )
   ch_plasmid_contigs = SEQKIT_GREP_PLASMID ( ch_filter_non_bacterial, ch_filter_plasmid)

   ch_non_plasmid_contigs = SEQKIT_GREP_NON_PLASMID.out.filter

   ch_versions = ch_versions.mix( SEQKIT_GREP_NON_PLASMID.out.versions )


    //
    // Module minimap2 index for plasmids and fungi
    //

    Channel.fromPath(params.eukaryotic_virus_genome)
        .map { file -> tuple(file.baseName, file) }
        .set { ch_eukaryotic_virus }



    ch_non_plasmid_index = MINIMAP2_FUNGI_EUKARYOTIC_VIRUS(ch_eukaryotic_virus ).index

    ch_non_plasmid_contigs_index = ch_non_plasmid_contigs.combine(ch_eukaryotic_virus)


    non_plasmid_index = ch_non_plasmid_contigs_index.map { meta, reads, db, index -> return [ meta, reads ] }
    non_plasmid_contigs = ch_non_plasmid_contigs_index.map { meta, reads, db, index -> return [ meta, index ] }




    MINIMAP2_ALIGN_EUKARYOTIC_VIRUS (
            non_plasmid_contigs,
            non_plasmid_index,
            true,
            false,
            false,
            false
        )


    /*ch_minimap2_mapped_eukaryotic = MINIMAP2_ALIGN_EUKARYOTIC_VIRUS.out.bam
        .map {
            meta, reads ->
                [ meta, reads, [] ]
        }.view()*/

    MSAMTOOLS_FILTER(MINIMAP2_ALIGN_EUKARYOTIC_VIRUS.out.bam).view()

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
