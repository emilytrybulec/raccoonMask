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
            REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, [], params.soft_mask)
        } else {
            REPEAT_MASKER(ch_consensus_fasta, ch_genome_fasta, params.species, params.soft_mask)
        }

        TWO_BIT(REPEAT_MASKER.out.fasta)

        REPEAT_VIEW(REPEAT_MASKER.out.align, TWO_BIT.out.out)
    }


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
