/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap       } from 'plugin/nf-validation'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { REPEATMODELER_BUILDDATABASE } from '../modules/nf-core/repeatmodeler/builddatabase/main' 
include { REPEATMODELER_REPEATMODELER } from '../modules/nf-core/repeatmodeler/repeatmodeler/main' 
include { REPEAT_MASKER } from '../modules/local/repeatmasker' 
include { REPEAT_MASKER_2 } from '../modules/local/repeatmasker2' 
include { TE_TRIMMER } from '../modules/local/tetrimmer' 
include { TWO_BIT } from '../modules/local/twoBit' 
include { REPEAT_VIEW } from '../modules/local/repeat_visualization' 
include { MC_HELPER } from '../modules/local/mchelper' 
include { genSample; warmupRepeatMasker; twoBit; genBatches; twoBittoFa; RepeatMasker; combineRMOUTOutput; combineRMAlignOutput } from '../modules/local/repeatmasker_faster' 


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REPEAT_CURATION {

    ch_db_fasta = Channel.fromPath(params.genome_fasta)
    ch_db_fasta
        .map { file -> tuple(id: file.baseName, file)  }
        .set { ch_genome_fasta }

    if (params.consensus_fasta == null) {
        REPEATMODELER_BUILDDATABASE(ch_genome_fasta)
        REPEATMODELER_REPEATMODELER(REPEATMODELER_BUILDDATABASE.out.db)
        ch_consensus_fasta = REPEATMODELER_REPEATMODELER.out.fasta
    } else { 
        ch_consensus = Channel.fromPath(params.consensus_fasta) 
        ch_consensus
            .map { file -> tuple(id: file.baseName, file)  }
            .set { ch_consensus_fasta }
    }

    if (params.te_trimmer == true){
        TE_TRIMMER(ch_consensus_fasta, ch_genome_fasta, params.cons_thr)
        
        if (params.repeat_masker == true){
            if(params.species == null){
                REPEAT_MASKER_2(TE_TRIMMER.out.fasta, ch_genome_fasta, [], params.soft_mask)
            } else {
                REPEAT_MASKER_2(TE_TRIMMER.out.fasta, ch_genome_fasta, params.species, params.soft_mask)
            }
        }
    } else if (params.MC_helper == true){
        MC_HELPER(ch_consensus_fasta, ch_genome_fasta, params.gene_ref)

        if (params.repeat_masker == true){
            if(params.species == null){
                REPEAT_MASKER_2(MC_HELPER.out.fasta, ch_genome_fasta, [], params.soft_mask)
            } else {
                REPEAT_MASKER_2(MC_HELPER.out.fasta, ch_genome_fasta, params.species, params.soft_mask)
            }
        }
    }

    if (params.repeat_masker == true){
        if(params.species == null){
            ch_species = Channel.empty()
        } else {ch_species = params.species}

        if (params.cluster == "xanadu"){
            genSample(ch_genome_fasta)
            warmupRepeatMasker(genSample.out.out, ch_species)
            twoBit(ch_genome_fasta)
            genBatches(warmupRepeatMasker.out.out, params.batchSize, twoBit.out.out)

            genBatches.out.bed
                .flatten()
                .set{ch_batches}

            ch_batches
                .combine(twoBit.out.out)
                .set{ch_batches_2bit}

            twoBittoFa(ch_batches_2bit)

            ch_consensus_fasta
                .combine(twoBittoFa.out.out)
                .set{ch_rm_batches}

            RepeatMasker(ch_rm_batches, ch_species, params.soft_mask)

            RepeatMasker.out.out
                .combine(twoBit.out.out)
                .set{ch_output_combine}

            RepeatMasker.out.align
                .combine(twoBit.out.out)
                .set{ch_align_combine}

            combineRMOUTOutput(ch_output_combine)
            combineRMAlignOutput(ch_align_combine, combineRMOUTOutput.out.trans)
            
            combineRMAlignOutput.out.align
                .set{repeatMasker_align}
            RepeatMasker.out.masked
                .set{repeatMasker_fasta}

        } else {
            if(params.species == null){
                REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, [], params.soft_mask)
            } else {
                REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, params.species, params.soft_mask)
            }
            repeatMasker_fasta = REPEAT_MASKER.out.fasta
            repeatMasker_align = REPEAT_MASKER.out.align

        TWO_BIT(repeatMasker_fasta)

        REPEAT_VIEW(repeatMasker_align, TWO_BIT.out.out)
       }


    }


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
