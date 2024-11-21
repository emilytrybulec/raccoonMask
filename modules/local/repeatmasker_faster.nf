process genSample {
  tag "$meta.id"
  label 'process_low'

  conda "bioconda::seqkit=2.4.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0':
        'biocontainers/seqkit:2.4.0--h9ee0642_0' }"

  input:
  tuple val(meta), path(genome_fasta)

  output:
  tuple val(meta), path("*.fasta"), emit: out

  script:
  """
  awk 'BEGIN {seq_id=""; seq=""; srand()} {if (\$0 ~ /^>/) {if (seq_id != "" && length(seq) >= 25000) {start = int(rand() * (length(seq) - 25000 + 1)) + 1; print seq_id; print substr(seq, start, 25000); exit} seq_id = \$0; seq = ""} else {seq = seq \$0}} END {if (seq_id != "" && length(seq) >= 25000) {start = int(rand() * (length(seq) - 25000 + 1)) + 1; print seq_id; print substr(seq, start, 25000)}}' $genome_fasta > sample.fasta

  """
}

process warmupRepeatMasker {
  tag "$meta.id"
  label 'process_mid'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.7p1--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.7p1--pl5321hdfd78af_1' }"

  input:
  tuple val(meta), path(small_seq)
  val species

  output:
  tuple val(meta), path("*.rmlog"), emit: out

  script:
  def species = species ? "-species ${species}" : ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  #
  # Run RepeatMasker with "-species" option on a small sequence in order to
  # force it to initialize the cached libraries.  Do not want to do this on the
  # cluster ( in parallel ) as it may cause each job to attempt the build at once.
  #
  RepeatMasker ${species} $small_seq >& ${prefix}.rmlog
  """
}

process twoBit {
    label 'process_low'

    conda "bioconda::ucsc-fatotwobit:469--he8037a5_2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-fatotwobit:469--he8037a5_2 ' :
        'quay.io/biocontainers/ucsc-fatotwobit:469--he8037a5_2 ' }"

    input:
    tuple val(meta), path(genomes)

    output:
    path("*.2bit"), emit: out

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${genomes}"
    """
    faToTwoBit -long $genomes ${genomes.baseName}.2bit
    """
}
process genBatches {
  tag "$meta.id"
  label 'process_mid'

  input:
  tuple val(meta), path(warmuplog)
  val batchSize
  file(inSeqFile)

  output:
  path("*bed") , emit: bed

  script:
  def prefix = task.ext.prefix ?: "${inSeqFile}"
  """
  perl ${projectDir}/assets/genBEDBatches.pl ${inSeqFile.baseName}.2bit $batchSize
  """
}
process twoBittoFa {
  label 'process_mid'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-twobittofa:472--h9b8f530_0' :
        'biocontainers/ucsc-twobittofa:472--h9b8f530_0' }"

  input:
  tuple file(batch_bed), file(inSeqFile)

  output:
  path("*.fa") , emit: out

  script:
  def prefix = task.ext.prefix ?: "${inSeqFile}"
  """
  
  twoBitToFa -bed=$batch_bed $inSeqFile ${batch_bed.baseName}.fa
  """
}

process RepeatMasker {
  tag "$meta.id"
  label 'process_mid'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.7p1--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.7p1--pl5321hdfd78af_1' }"

  input:
  tuple val(meta), path(curation_fasta), path(batch_file)
  val species
  val soft_mask

  output:
  tuple val(meta), path("*.fa.out") , emit: out
  tuple val(meta), path("*.fa.align") , emit: align

  script:
  def species = species ? "-species ${species}" : ''
  def soft_mask = soft_mask ? "-xsmall" : ''
  """
  #
  # Run RepeatMasker and readjust coordinates
  #
  RepeatMasker -s -e ncbi -lib $curation_fasta -pa $task.cpus -a $soft_mask ${batch_file.baseName}.fa >& ${batch_file.baseName}.rmlog
  ${projectDir}/assets/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.out 
  ${projectDir}/assets/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.align
  cp ${batch_file.baseName}.fa.out ${batch_file.baseName}.fa.out.unadjusted
  mv ${batch_file.baseName}.fa.out.adjusted ${batch_file.baseName}.fa.out
  mv ${batch_file.baseName}.fa.align.adjusted ${batch_file.baseName}.fa.align
  """
}

process combineRMOUTOutput {

  input:
  tuple file(twoBitFile), file(outfiles) 

  output:
  file("*.rmout.gz")
  file("*.summary")
  file("combOutSorted-translation.tsv") 
  
  script:
  """
  for f in ${outfiles}; do cat \$f >> combOut; done
  echo "   SW   perc perc perc  query     position in query    matching          repeat       position in repeat" > combOutSorted
  echo "score   div. del. ins.  sequence  begin end   (left)   repeat            class/family begin  end    (left)  ID" >> combOutSorted
  grep -v -e "^\$" combOut | sort -k5,5 -k6,6n -T ${workflow.workDir} >> combOutSorted
  ${projectDir}/assets/renumberIDs.pl combOutSorted > combOutSortedRenumbered
  mv translation-out.tsv combOutSorted-translation.tsv
  export PATH=${ucscToolsDir}/\$PATH
  ${repeatMaskerDir}/util/buildSummary.pl -genome ${twoBitFile} -useAbsoluteGenomeSize combOutSortedRenumbered > ${twoBitFile.baseName}.summary
  gzip -c combOutSortedRenumbered > ${twoBitFile.baseName}.rmout.gz
  """
}

process combineRMAlignOutput {

  input:
  tuple file(twoBitFile), file(alignfiles) 
  file transFile 
  
  output:
  file("*.rmalign.gz")

  script:
  """
  for f in ${alignfiles}; do cat \$f >> combAlign; done
  ####${projectDir}/assets/alignToBed.pl -fullAlign combAlign | ${ucscToolsDir}/bedSort stdin stdout | ${workflow.projectDir}/bedToAlign.pl > combAlign-sorted
  ${projectDir}/assets/alignToBed.pl -fullAlign combAlign > tmp.bed
  # Be mindful of this buffer size...should probably make this a parameter
  sort -k1,1V -k2,2n -k3,3nr -S 3G -T ${workflow.workDir} tmp.bed > tmp.bed.sorted
  ${projectDir}/assets/bedToAlign.pl tmp.bed.sorted > combAlign-sorted
  ${projectDir}/assets/renumberIDs.pl -translation ${transFile} combAlign-sorted > combAlign-sorted-renumbered
  gzip -c combAlign-sorted-renumbered > ${twoBitFile.baseName}.rmalign.gz
  """
}
