
include { GAWK                     as  GAWK_SINGLE                   } from '../../../modules/local/gawk/main.nf'
include { GAWK                     as  GAWK_BOTH                     } from '../../../modules/local/gawk/main.nf'
include { KRONA_KTIMPORTTAXONOMY   as  KRONA_KTIMPORTTAXONOMY_SINGLE } from '../../../modules/nf-core/krona/ktimporttaxonomy/main.nf'
include { KRONA_KTIMPORTTAXONOMY   as  KRONA_KTIMPORTTAXONOMY_BOTH   } from '../../../modules/nf-core/krona/ktimporttaxonomy/main.nf'
include { RANALYSIS                                                  } from '../../../modules/local/ranalysis/main.nf'


workflow REPORTING {

    take:
    classified_reads_single
    classified_reads_both
    integration_sites
    taxonomy
    sampleids

    main:

    ch_versions = channel.empty()

    GAWK_SINGLE ( classified_reads_single )
    ch_versions = ch_versions.mix(GAWK_SINGLE.out.versions)

    GAWK_BOTH ( classified_reads_both )
    ch_versions = ch_versions.mix(GAWK_BOTH.out.versions)

    // GAWK now emits [ val(meta), path(collated) ] — extract file and re-wrap with group meta for Krona
    fakemeta = channel.value([id: "group"])
    single_input = fakemeta.combine(GAWK_SINGLE.out.collated_reads.map { _meta, f -> f }.collect())

    KRONA_KTIMPORTTAXONOMY_SINGLE ( single_input, taxonomy )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY_SINGLE.out.versions)

    both_input = fakemeta.combine(GAWK_BOTH.out.collated_reads.map { _meta, f -> f }.collect())

    KRONA_KTIMPORTTAXONOMY_BOTH ( both_input, taxonomy )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY_BOTH.out.versions)

    ch_rmarkdown = channel.value(file("$projectDir/assets/analysis_report.Rmd"))

    // Strip meta from tuples and collect all files before passing to RANALYSIS
    ch_single_files = classified_reads_single.map { meta, f -> f }.collect()
    ch_both_files   = classified_reads_both.map { meta, f -> f }.collect()

    RANALYSIS ( ch_single_files, ch_both_files, integration_sites, sampleids, ch_rmarkdown, params.istest, params.taxonomy_id)

    emit:
    single_html = KRONA_KTIMPORTTAXONOMY_SINGLE.out.html
    versions    = ch_versions
}
