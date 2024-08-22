process MSAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::msamtools=1.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msamtools:1.1.3--he4a0461_0':
        'biocontainers/msamtools:1.1.3--he4a0461_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*filtered.bam")               , emit: filtered_bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    msamtools filter \\
        $args \\
        -t $task.cpus \\
        -b \\
        -l 1000 \\
        -p 95 \\
        -z 50 \\
        $reads \\
        > ${prefix}_filtered.bam
    """


}
