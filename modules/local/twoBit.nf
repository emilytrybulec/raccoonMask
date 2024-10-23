process TWO_BIT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-fatotwobit:469--he8037a5_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-fatotwobit:469--he8037a5_2 ' :
        'quay.io/biocontainers/ucsc-fatotwobit:469--he8037a5_2 ' }"

    input:
    tuple val(meta), path(repeats)

    output:
    tuple val(meta), path("*.2bit"), emit: out

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    faToTwoBit $repeats ${prefix}.2bit
    """
}
