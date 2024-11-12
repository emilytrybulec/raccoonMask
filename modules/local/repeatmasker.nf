process REPEAT_MASKER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.7p1--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.7p1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(curation_fasta)
    tuple val(meta), path(genome_fasta)
    val species
    val soft_mask

    output:
    tuple val(meta), path("*.out"), emit: out
    tuple val(meta), path("*.align") , emit: align
    tuple val(meta), path("*.tbl") , emit: table
    tuple val(meta), path("*.masked") , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "repeatmask_${meta.id}"
    def soft_mask = soft_mask ? "-xsmall" : ''
    def species = species ? "-species ${species}" : ''
    """
    RepeatMasker -s \\
          -e ncbi \\
          -lib $curation_fasta \\
          $soft_mask \\
          $species \\
          -pa $task.cpus \\
	  -a \\
          $genome_fasta
    """
}
