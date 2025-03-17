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
include { REPEAT_MASKER_2 as PIPE_REPEAT_MASKER} from '../modules/local/repeatmasker2' 
include { TE_TRIMMER } from '../modules/local/tetrimmer' 
include { TWO_BIT } from '../modules/local/twoBit' 
include { REPEAT_VIEW } from '../modules/local/repeat_visualization' 
include { MC_HELPER } from '../modules/local/mchelper' 
include { genSample; warmupRepeatMasker; twoBit; genBatches; twoBittoFa; RepeatMasker; adjCoordinates; combineRMOUTOutput; combineRMAlignOutput; makeMaskedFasta } from '../modules/local/repeatmasker_faster' 
include { PIPERepeatMasker; PIPEadjCoordinates; PIPEcombineRMOUTOutput; PIPEcombineRMAlignOutput; PIPEmakeMaskedFasta } from '../modules/local/pipe_repeatmasker_faster' 


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

    if (params.consensus_fasta == null && params.libdir == null) {
        REPEATMODELER_BUILDDATABASE(ch_genome_fasta)
        REPEATMODELER_REPEATMODELER(REPEATMODELER_BUILDDATABASE.out.db)
        ch_consensus_fasta = REPEATMODELER_REPEATMODELER.out.fasta
    } else if (params.consensus_fasta != null){ 
        ch_consensus = Channel.fromPath(params.consensus_fasta) 
        ch_consensus
            .map { file -> tuple(id: file.baseName, file)  }
            .set { ch_consensus_fasta }
    } else {
        ch_consensus_fasta = Channel.fromPath(params.libdir) 
    }

    if (params.libdir == null){
        ch_libdir = Channel.empty()
    } else {
        ch_libdir = Channel.fromPath(params.libdir) 
    }

    if (params.te_trimmer == true){
        TE_TRIMMER(ch_consensus_fasta, ch_genome_fasta, params.cons_thr)
        
        if (params.repeat_masker == true){
            if(params.species == null){
                PIPE_REPEAT_MASKER(TE_TRIMMER.out.fasta, ch_genome_fasta, [], params.soft_mask) 
            } else {
                PIPE_REPEAT_MASKER(TE_TRIMMER.out.fasta, ch_genome_fasta, params.species, params.soft_mask)      
            }
        pipe_repeatMasker_fasta = PIPE_REPEAT_MASKER.out.fasta
        pipe_repeatMasker_align = PIPE_REPEAT_MASKER.out.align
        } else{ 
        pipe_repeatMasker_fasta = Channel.empty()
        pipe_repeatMasker_align = Channel.empty()
        }
    } else if (params.MC_helper == true){
        MC_HELPER(ch_consensus_fasta, ch_genome_fasta, params.gene_ref)

        if (params.repeat_masker == true){
            if(params.species == null){
                ch_species = Channel.empty()
            } else {ch_species = params.species}

            if (params.cluster == "xanadu" || params.cluster == "mantis"){

            genSample(ch_genome_fasta)
            warmupRepeatMasker(PIPEgenSample.out.out, ch_species)
            twoBit(ch_genome_fasta)
            genBatches(warmupRepeatMasker.out.out, params.batchSize, twoBit.out.out)

            genBatches.out.bed
                .flatten()
                .set{ch_batches}

            ch_batches
                .combine(twoBit.out.out)
                .set{ch_batches_2bit}

            twoBittoFa(ch_batches_2bit)

            twoBittoFa.out.out
                .flatten()
                .map{ file -> tuple(file.baseName, file) }
                .set{batches_meta}

            MC_HELPER.out.fasta
                .map{it[1]}
                .set{pipe_consensus_nometa}

            batches_meta
                .combine(pipe_consensus_nometa)
                .set{pipe_ch_rm_batches}

            PIPERepeatMasker(pipe_ch_rm_batches, ch_species, params.soft_mask)

            ch_batches
                .flatten()
                .map{ file -> tuple(file.baseName, file) }
                .set{pipe_batches_bed_meta}
                
            pipe_batches_bed_meta
                .join(PIPERepeatMasker.out.out)
                .join(PIPERepeatMasker.out.align)
                .set{pipe_ch_rmout}

            PIPEadjCoordinates(pipe_ch_rmout)

            PIPEadjCoordinates.out.out
                .map {it[1]}
                .set{pipe_out_nometa}
            PIPEadjCoordinates.out.align
                .map {it[1]}
                .set{pipe_align_nometa}

            pipe_out_nometa
                .collect()
                .set{pipe_ch_out}
            pipe_align_nometa
                .collect()
                .set{pipe_ch_align}

            PIPEtwoBit.out.out
                .map { file -> tuple(id: file.baseName, file)  }
                .set{pipe_twoBit_meta}

            PIPEcombineRMOUTOutput(pipe_twoBit_meta, pipe_ch_out)
            PIPEcombineRMAlignOutput(pipe_twoBit_meta, pipe_ch_align, PIPEcombineRMOUTOutput.out.trans)

            PIPEcombineRMAlignOutput.out.align
                .set{pipe_repeatMasker_align}

            ch_genome_fasta
                .combine(PIPEcombineRMOUTOutput.out.bed)
                .set{pipe_masked_ch}

            PIPEmakeMaskedFasta(pipe_masked_ch, ch_species)

            PIPEmakeMaskedFasta.out.masked
                .set{pipe_repeatMasker_fasta}
            } else {
                if(params.species == null){
                    PIPE_REPEAT_MASKER(MC_HELPER.out.fasta, ch_genome_fasta, [], params.soft_mask)
                } else {
                    PIPE_REPEAT_MASKER(MC_HELPER.out.fasta, ch_genome_fasta, params.species, params.soft_mask)
                }
            }
            pipe_repeatMasker_fasta = PIPE_REPEAT_MASKER.out.fasta
            pipe_repeatMasker_align = PIPE_REPEAT_MASKER.out.align
        }
        } else {
        pipe_repeatMasker_fasta = Channel.empty()
        pipe_repeatMasker_align = Channel.empty()
        }
    

    if (params.repeat_masker == true){
        if(params.species == null){
            ch_species = Channel.empty()
        } else {ch_species = params.species}

        if (params.cluster == "xanadu" || params.cluster == "mantis"){
            if (params.MC_helper == false){
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

            twoBittoFa.out.out
                .flatten()
                .map{ file -> tuple(file.baseName, file) }
                .set{batches_meta}
            } 
            ch_consensus_fasta
                .map{it[1]}
                .set{consensus_nometa}

            if (params.consensus_fasta != null){
            batches_meta
                .combine(consensus_nometa)
                .set{ch_rm_batches} }
            else {
            batches_meta
                .combine(ch_libdir)
                .set{ch_rm_batches}
            }

            if (params.libdir == null){
                RepeatMasker(ch_rm_batches, ch_species, params.soft_mask)
            } else {
                RepeatMasker(ch_rm_batches, ch_species, params.soft_mask)
            }

            ch_batches
                .flatten()
                .map{ file -> tuple(file.baseName, file) }
                .set{batches_bed_meta}
                
            batches_bed_meta
                .join(RepeatMasker.out.out)
                .join(RepeatMasker.out.align)
                .set{ch_rmout}

            adjCoordinates(ch_rmout)

            adjCoordinates.out.out
                .map {it[1]}
                .set{out_nometa}
            adjCoordinates.out.align
                .map {it[1]}
                .set{align_nometa}

            out_nometa
                .collect()
                .set{ch_out}
            align_nometa
                .collect()
                .set{ch_align}

            twoBit.out.out
                .map { file -> tuple(id: file.baseName, file)  }
                .set{twoBit_meta}

            combineRMOUTOutput(twoBit_meta, ch_out)
            combineRMAlignOutput(twoBit_meta, ch_align, combineRMOUTOutput.out.trans)

            combineRMAlignOutput.out.align
                .set{repeatMasker_align}

            ch_genome_fasta
                .combine(combineRMOUTOutput.out.bed)
                .set{masked_ch}

            makeMaskedFasta(masked_ch, ch_species)

            makeMaskedFasta.out.masked
                .set{repeatMasker_fasta}
            
        } else {
            if(params.species == null){
                if (params.libdir == null){
                    REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, [], params.soft_mask, [])
                } else {
                    REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, [], params.soft_mask, params.libdir)
                } 
            } else {
                if (params.libdir == null){
                    REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, params.species, params.soft_mask, [])
                } else {
                    REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, params.species, params.soft_mask, params.libdir)
                }
            }
            repeatMasker_fasta = REPEAT_MASKER.out.fasta
            repeatMasker_align = REPEAT_MASKER.out.align
       }
        repeatMasker_fasta
            .concat(pipe_repeatMasker_fasta)
            .set{ch_Masker_fasta}

        repeatMasker_align
            .concat(pipe_repeatMasker_align)
            .set{ch_Masker_align}
        ch_Masker_align.view()
        
        TWO_BIT(ch_Masker_fasta)
        TWO_BIT.out.out.view()

        ch_Masker_align
            .join(TWO_BIT.out.out)
            .set{ch_repeat_view}

        REPEAT_VIEW(ch_repeat_view)
    }


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
