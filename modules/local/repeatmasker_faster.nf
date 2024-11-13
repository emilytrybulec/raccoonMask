process warmupRepeatMasker {

  input:
  path(small_seq) from warmup_chan

  output:
  path("*.rmlog") into done_warmup_chan

  script:
  """
  #
  # Run RepeatMasker with "-species" option on a small sequence in order to
  # force it to initialize the cached libraries.  Do not want to do this on the
  # cluster ( in parallel ) as it may cause each job to attempt the build at once.
  #
  hostname > node
  ${repeatMaskerDir}/RepeatMasker ${otherOptions} ${small_seq.baseName}.fa >& ${small_seq.baseName}.rmlog
  """
}

process genBatches {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  path(warmuplog) from done_warmup_chan
  each file(inSeqFile) from inSeqFiles

  output:
  tuple file("${inSeqFile.baseName}.2bit"), file("batch*.bed") into batchChan 

  script:
  """
  # Generate 2bit files if necessary
  if [ ${inSeqFile.extension} == "gz" ]; then
    gunzip -c ${inSeqFile} | ${ucscToolsDir}/faToTwoBit -long stdin ${inSeqFile.baseName}.2bit
  elif [ ${inSeqFile.extension} == "2bit" ]; then
    # Ah....the luxury of 2bit
    sleep 0
  else
    ${ucscToolsDir}/faToTwoBit -long ${inSeqFile} ${inSeqFile.baseName}.2bit
  fi  
 
  # This can magically accept FASTA, Gzip'd FASTA, or 2BIT...but 2Bit is prob. fastest
  export UCSCTOOLSDIR=${ucscToolsDir}
  ${workflow.projectDir}/genBEDBatches.pl ${inSeqFile.baseName}.2bit ${batchSize}
  """
}

process RepeatMasker {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  file inLibFile from opt_libFile
  set path(inSeqTwoBitFile), path(batch_file) from batchChan.transpose()

  output:
  tuple path(inSeqTwoBitFile), path("${batch_file.baseName}.fa.out") into rmoutChan
  tuple path(inSeqTwoBitFile), path("${batch_file.baseName}.fa.align") into rmalignChan

  script:
  def libOpt = inLibFile.name != 'NO_FILE' ? "-lib $inLibFile" : "-species '" + species + "'"
  """
  #
  # Run RepeatMasker and readjust coordinates
  #
  ${ucscToolsDir}/twoBitToFa -bed=${batch_file} ${inSeqTwoBitFile} ${batch_file.baseName}.fa
  ${repeatMaskerDir}/RepeatMasker -pa ${pa_param} -a ${otherOptions} ${libOpt} ${batch_file.baseName}.fa >& ${batch_file.baseName}.rmlog
  export REPEATMASKER_DIR=${repeatMaskerDir}
  ${workflow.projectDir}/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.out 
  ${workflow.projectDir}/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.align
  cp ${batch_file.baseName}.fa.out ${batch_file.baseName}.fa.out.unadjusted
  mv ${batch_file.baseName}.fa.out.adjusted ${batch_file.baseName}.fa.out
  mv ${batch_file.baseName}.fa.align.adjusted ${batch_file.baseName}.fa.align
  """
}

process combineRMOUTOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisAdjOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy'

  input:
  tuple file(twoBitFile), file(outfiles) from rmoutChan.map { tb, outf -> [ tb.toRealPath(), outf ]}.groupTuple()

  output:
  file("*.rmout.gz")
  file("*.summary")
  file("combOutSorted-translation.tsv") into rmAlignTransChan
  
  script:
  """
  for f in ${outfiles}; do cat \$f >> combOut; done
  echo "   SW   perc perc perc  query     position in query    matching          repeat       position in repeat" > combOutSorted
  echo "score   div. del. ins.  sequence  begin end   (left)   repeat            class/family begin  end    (left)  ID" >> combOutSorted
  grep -v -e "^\$" combOut | sort -k5,5 -k6,6n -T ${workflow.workDir} >> combOutSorted
  ${workflow.projectDir}/renumberIDs.pl combOutSorted > combOutSortedRenumbered
  mv translation-out.tsv combOutSorted-translation.tsv
  export PATH=${ucscToolsDir}/\$PATH
  ${repeatMaskerDir}/util/buildSummary.pl -genome ${twoBitFile} -useAbsoluteGenomeSize combOutSortedRenumbered > ${twoBitFile.baseName}.summary
  gzip -c combOutSortedRenumbered > ${twoBitFile.baseName}.rmout.gz
  """
}

process combineRMAlignOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisAdjOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy'

  input:
  tuple file(twoBitFile), file(alignfiles) from rmalignChan.map { tb, alignf -> [ tb.toRealPath(), alignf ]}.groupTuple()
  file transFile from rmAlignTransChan
  
  output:
  file("*.rmalign.gz")

  script:
  """
  for f in ${alignfiles}; do cat \$f >> combAlign; done
  ####${workflow.projectDir}/alignToBed.pl -fullAlign combAlign | ${ucscToolsDir}/bedSort stdin stdout | ${workflow.projectDir}/bedToAlign.pl > combAlign-sorted
  ${workflow.projectDir}/alignToBed.pl -fullAlign combAlign > tmp.bed
  # Be mindful of this buffer size...should probably make this a parameter
  sort -k1,1V -k2,2n -k3,3nr -S 3G -T ${workflow.workDir} tmp.bed > tmp.bed.sorted
  ${workflow.projectDir}/bedToAlign.pl tmp.bed.sorted > combAlign-sorted
  ${workflow.projectDir}/renumberIDs.pl -translation ${transFile} combAlign-sorted > combAlign-sorted-renumbered
  gzip -c combAlign-sorted-renumbered > ${twoBitFile.baseName}.rmalign.gz
  """
}
