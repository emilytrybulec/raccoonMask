process genSample {
  tag "$meta.id"
  label 'process_low'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--hed695b0_5':
        'biocontainers/bioawk:1.0--hed695b0_5' }"

  input:
  tuple val(meta), path(genome_fasta)

  output:
  tuple val(meta), path("*.fasta"), emit: out

  script:
  """
  bioawk -c fastx '{if (length(\$seq) >= 25000) {start = int(rand() * (length(\$seq) - 25000 + 1)) + 1; print ">" \$name; print substr(\$seq, start, 25000); exit}}' $genome_fasta > sample.fasta

  """
}

process warmupRepeatMasker {
  tag "$meta.id"
  label 'process_medium'

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
  label 'process_medium'

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
  label 'process_low'

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
  tag "$meta"
  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.7p1--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.7p1--pl5321hdfd78af_1' }"

  input:
  tuple val(meta), path(batch_file), path(curation_fasta)
  val species
  val soft_mask

  output:
  tuple val(meta), path("*.out") , emit: out
  tuple val(meta), path("*.align") , emit: align
  tuple val(meta), path("*.masked") , emit: masked

  script:
  def species = species ? "-species ${species}" : ''
  def soft_mask = soft_mask ? "-xsmall" : ''
  def lib = "${curation_fasta}"
  def libOpt = lib.contains('.fa') ? "-lib $curation_fasta" : "-species ${species}"
  """
  #
  # Run RepeatMasker
  #

  RepeatMasker -s -e ncbi $libOpt -pa $task.cpus -a $soft_mask ${batch_file.baseName}.fa >| ${batch_file.baseName}.rmlog 2>&1
  """
}
process adjCoordinates {
  tag "$meta"
  label 'process_low'

  input:
  tuple val(meta), path(batch_file), path(out), path(align)

  output:
  tuple val(meta), path("*.out.adjusted") , emit: out
  tuple val(meta), path("*.align.adjusted") , emit: align

  script:

  """
  #
  # Adjust Output
  #

  ${projectDir}/assets/adjCoordinates.pl ${batch_file} ${out}
  ${projectDir}/assets/adjCoordinates.pl ${batch_file} ${align}
  cp ${out} ${batch_file.baseName}.fa.out.unadjusted
  """
}

process combineRMOUTOutput {
  label 'process_medium'

  input:
  tuple val(meta), file(twoBitFile)
  path(outfiles)

  output:
  tuple val(meta), path("*.rmout.gz"), emit: out
  tuple val(meta), path("*.summary"), emit: summary
  path("*.bed"), emit: bed
  path("combOutSorted-translation.tsv"), emit: trans 
  
  script:
  """
  cp ${twoBitFile} ./local.2bit
  for f in ${outfiles}; do cat \$f >> combOut; done
  echo "   SW   perc perc perc  query     position in query    matching          repeat       position in repeat" > combOutSorted
  echo "score   div. del. ins.  sequence  begin end   (left)   repeat            class/family begin  end    (left)  ID" >> combOutSorted
  grep -v -e "^\$" combOut | sort -k5,5 -k6,6n -T ${workflow.workDir} >> combOutSorted
  ${projectDir}/assets/renumberIDs.pl combOutSorted > combOutSortedRenumbered
  mv translation-out.tsv combOutSorted-translation.tsv
  /core/labs/Oneill/jstorer/RepeatMasker/util/buildSummary.pl -genome local.2bit -useAbsoluteGenomeSize combOutSortedRenumbered > ${twoBitFile.baseName}.summary
  gzip -c combOutSortedRenumbered > ${twoBitFile.baseName}.rmout.gz
  ${projectDir}/assets/rm_to_bed.py combOutSortedRenumbered ${twoBitFile.baseName}.bed
  """
}

process combineRMAlignOutput {
  label 'process_medium'

  input:
  tuple val(meta), file(twoBitFile)
  path(alignfiles)
  file transFile 
  
  output:
  tuple val(meta), path("*.rmalign.gz"), emit: align

  script:
  """
  for f in ${alignfiles}; do cat \$f >> combAlign; done
  ${projectDir}/assets/alignToBed.pl -fullAlign combAlign > tmp.bed
  # Be mindful of this buffer size...should probably make this a parameter
  sort -k1,1V -k2,2n -k3,3nr -S 3G -T ${workflow.workDir} tmp.bed > tmp.bed.sorted
  ${projectDir}/assets/bedToAlign.pl tmp.bed.sorted > combAlign-sorted
  ${projectDir}/assets/renumberIDs.pl -translation ${transFile} combAlign-sorted > combAlign-sorted-renumbered
  gzip -c combAlign-sorted-renumbered > ${twoBitFile.baseName}.rmalign.gz
  """
}

process makeMaskedFasta {
  label 'process_low'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
     'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
     'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

  input:
  tuple val(meta), path(genome), path(bed)
  val soft_mask
  
  output:
  tuple val(meta), path("*.masked"), emit: masked

  script:
  def soft_mask_opt = soft_mask ? "-soft" : ''
  """
  bedtools maskfasta -fi $genome -bed $bed $soft_mask_opt -fo ${bed.baseName}.fa.masked
  """
}
