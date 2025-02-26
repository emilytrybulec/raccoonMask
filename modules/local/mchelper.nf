process MC_HELPER {
    tag "$meta.id"
    label 'process_high'

    container 'docker://plantgenomics/mchelper:1.6.6.0'

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
    set +e  # Disable 'exit on error'
    set +u  # Disable 'treat unset variables as errors'
    set +o pipefail  # Allow pipelines to continue
    set +C  # Disable 'no clobber'

    source /opt/conda/etc/profile.d/conda.sh
    conda activate MCHelper

    python3 /opt/MCHelper/MCHelper.py \\
        -l $lib \\
        -o ${prefix} \\
        -g $genome \\
        --input_type fasta \\
        $genes \\
        -a F -z 10 -c 3 -v Y
    """
}
