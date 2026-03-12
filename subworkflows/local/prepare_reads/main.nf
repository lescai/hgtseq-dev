
// this subworkflow prepares the inputs from fastq reads to bam files
// and performs QC of both reads and resulting bam files

// modules to include in this subworkflow

include { BWAMEM2_INDEX                } from '../../../modules/nf-core/bwamem2/index/main.nf'
include { BWAMEM2_MEM                  } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_INDEX   as BWAMEM1_INDEX } from '../../../modules/nf-core/bwa/index/main.nf'
include { BWA_MEM     as BWAMEM1_MEM   } from '../../../modules/nf-core/bwa/mem/main.nf'
include { TRIMGALORE                   } from '../../../modules/nf-core/trimgalore/main.nf'


workflow PREPARE_READS {

    take:
    reads      // channel: [mandatory] [ val(meta), [ reads ] ]
    fasta      // channel: [mandatory] /path/to/reference/fasta
    aligner    // string:  [mandatory] "bwa-mem" or "bwa-mem2"

    main:
    ch_versions = channel.empty()
    aligned_bam = channel.empty()

    fasta_meta = channel.value(file(fasta)).map { f -> [[id: f.baseName], f] }

    TRIMGALORE ( reads )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions_trimgalore)

    if (aligner == "bwa-mem") {
        // reference is indexed if index not available in iGenomes
        BWAMEM1_INDEX ( fasta_meta )
        ch_versions = ch_versions.mix(BWAMEM1_INDEX.out.versions_bwa)

        // sets bwaindex to correct input
        bwaindex      = params.fasta ? params.bwaindex      ? channel.fromPath(params.bwaindex).collect().map{ it -> [[id:it[0].baseName], it] } : BWAMEM1_INDEX.out.index : []

        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        BWAMEM1_MEM ( TRIMGALORE.out.reads, bwaindex, [[],[]], false )
        ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions_bwa)
        ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions_samtools)
        aligned_bam = BWAMEM1_MEM.out.bam
    } else {
        // reference is indexed if index not available in iGenomes
        BWAMEM2_INDEX ( fasta_meta )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions_bwamem2)

        // sets bwamem2index to correct input
        bwamem2index  = params.fasta ? params.bwamem2index  ? channel.fromPath(params.bwamem2index).collect().map{ it -> [[id: it[0].baseName], it] } : BWAMEM2_INDEX.out.index : []

        // appropriately tagged interleaved FASTQ reads are mapped to the reference
        BWAMEM2_MEM ( TRIMGALORE.out.reads, bwamem2index, [[],[]], false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions_bwamem2)
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions_samtools)
        aligned_bam = BWAMEM2_MEM.out.bam
    }



    emit:
    trimmed_reads = TRIMGALORE.out.reads  // channel: [mandatory] [ val(meta), [ reads ] ]
    bam           = aligned_bam           // channel [mandatory] [ val(meta), [ bam ] ]
    versions      = ch_versions           // channel: [ versions.yml ]


}
