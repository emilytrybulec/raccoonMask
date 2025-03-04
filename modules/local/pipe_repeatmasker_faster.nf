process PIPERepeatMasker {
  tag "$meta"
  label 'process_medium'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.7p1--pl5321hdfd78af_1' :
        'biocontainers/repeatmasker:4.1.7p1--pl5321hdfd78af_1' }"

  input:
  tuple val(meta), path(batch_file), path(curation_fasta)
  val species
  val soft_mask
  path(libDir)

  output:
  tuple val(meta), path("*.out") , emit: out
  tuple val(meta), path("*.align") , emit: align
  tuple val(meta), path("*.masked") , emit: masked

  script:
  def species = species ? "-species ${species}" : ''
  def soft_mask = soft_mask ? "-xsmall" : ''
  def libdir = "${libDir}"
  def libOpt = libdir ? "-libdir ${libdir} -species ${species}" : lib.contains('.fa') ? "-lib ${curation_fasta}" :  "-species ${species}"  """
  #
  # Run RepeatMasker
  #

  RepeatMasker -s -e ncbi $libOpt -pa $task.cpus -a $soft_mask ${batch_file.baseName}.fa >| ${batch_file.baseName}.rmlog 2>&1
  """
}
process PIPEadjCoordinates {
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

process PIPEcombineRMOUTOutput {
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

process PIPEcombineRMAlignOutput {
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

process PIPEmakeMaskedFasta {
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
