process SEMIBIN_SINGLEEASYBIN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/semibin:2.1.0--pyhdfd78af_0':
        'biocontainers/semibin:2.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(bam)

    output:
    tuple val(meta), path("*.csv")                        , emit: csv
    tuple val(meta), path("*.h5")                         , emit: model            ,optional: true
    tuple val(meta), path("output_prerecluster_bins/*.fa"), emit: output_fasta     ,optional: true
    tuple val(meta), path("output_recluster_bins/*.fa")   , emit: recluster_fasta  ,optional: true
    tuple val(meta), path("*.tsv")                        , emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args3 = task.ext.args3 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    SemiBin2 \\
        single_easy_bin \\
        $args \\
        -i $fasta \\
        -b $bam \\
        -o $prefix \\
        -t $task.cpus \\
        $args3

    mv $prefix/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$( SemiBin2 --version )
    END_VERSIONS
"""
}
