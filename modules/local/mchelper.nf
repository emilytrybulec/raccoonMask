process MC_HELPER {
    tag "$meta.id"
    label 'process_low'

    container 'docker://emilytrybulec/mchelper6'

    input:
    tuple val(meta), path(lib)
    tuple val(meta), path(genome)
    path(ref_genes)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genes = ref_genes ? "-b $ref_genes" : ""
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate MCHelper

    python3 /opt/MCHelper/MCHelper.py \\
        -l $lib \\
        -o ${prefix} \\
        -g $genome \\
        --input_type fasta \\
        -b $genes \\
        -a F -z 10 -c 3 -v Y
    """
}
