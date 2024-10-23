process REPEAT_MASKER_2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.7p1--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.7p1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(te_curation_fasta)
    tuple val(meta), path(genome_fasta)
    val species
    val soft_mask

    output:
    tuple val(meta), path("*.out"), emit: out
    tuple val(meta), path("*.align") , emit: align
    tuple val(meta), path("*.tbl") , emit: table
    tuple val(meta), path("*.fasta.masked") , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "te_repeatmask_${meta.id}"
    def soft_mask = soft_mask ? "-xsmall" : ''
    def species = species ? "-species ${species}" : ''
    """
    RepeatMasker -s \\
          -lib $te_curation_fasta \\
          -pa $task.cpus \\
          $soft_mask \\
          $species \\
          $genome_fasta -a
    """
}
