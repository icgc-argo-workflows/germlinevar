// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { DEEPVARIANT                                                   } from '../../../modules/nf-core/deepvariant/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF                      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_INDEL                      } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_SNV                        } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF                       } from '../../../modules/nf-core/tabix/tabix/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_INDEL_GERMLINEVARIANT  } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_SNV_GERMLINEVARIANT    } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_INDEL                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_SNV                    } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_DEEPVARIANT {

    take:
        cram                     // channel: [mandatory] [meta, cram, crai, interval]
        dict                     // channel: [optional]
        fasta                    // channel: [mandatory]
        fasta_fai                // channel: [mandatory]
        analysis_json
        versions

    main:

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(versions)

    //Run caller
    DEEPVARIANT(
        cram,     // [meta, cram , crai , intervals]
        fasta,    // [meta , fasta]
        fasta_fai, // [meta , fasta]
        [[id:null],file("NO_FILE")] // [meta , gzi]
        )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())


    //Separate based on intervals
    DEEPVARIANT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }.set{deepvariant_vcf_out}

    // Only when no intervals, compress and index
    TABIX_VC_DEEPVARIANT_VCF(deepvariant_vcf_out.no_intervals)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    // When intervals, merge, compress and index
    MERGE_DEEPVARIANT_VCF(
        deepvariant_vcf_out.intervals
            .map{ meta, vcf ->

                new_meta = [
                            id: meta.id,
                            experimentalStrategy : meta.experimentalStrategy,
                            genomeBuild : meta.genomeBuild,
                            tumourNormalDesignation : meta.tumourNormalDesignation,
                            sampleType : meta.sampleType ,
                            gender : meta.gender,
                            study_id : meta.study_id,
                            num_intervals:  meta.num_intervals,
                        ]
                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})
    ch_versions = ch_versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)

    //Combine single and multi interval VCF channels
    combine_deepvariant_vcf = Channel.empty().mix(
                        MERGE_DEEPVARIANT_VCF.out.vcf,
                        deepvariant_vcf_out.no_intervals,
                        )
                    .map{ meta, vcf->
                        [[
                            id: meta.id,
                            experimentalStrategy : meta.experimentalStrategy,
                            genomeBuild : meta.genomeBuild,
                            tumourNormalDesignation : meta.tumourNormalDesignation,
                            sampleType : meta.sampleType ,
                            gender : meta.gender,
                            study_id : meta.study_id,
                            num_intervals:  meta.num_intervals,
                            tool:  "deepvariant"
                        ], vcf]
                    }

    combine_deepvariant_vcf_tbi = Channel.empty().mix(
                        MERGE_DEEPVARIANT_VCF.out.tbi,
                        TABIX_VC_DEEPVARIANT_VCF.out.tbi,
                        )
                    .map{ meta, tbi ->
                        [[
                            id: meta.id,
                            experimentalStrategy : meta.experimentalStrategy,
                            genomeBuild : meta.genomeBuild,
                            tumourNormalDesignation : meta.tumourNormalDesignation,
                            sampleType : meta.sampleType ,
                            gender : meta.gender,
                            study_id : meta.study_id,
                            num_intervals:  meta.num_intervals,
                            tool:  "deepvariant"
                        ], tbi]
                    }

    // Recombine VCF and TBI to make a single payload
    select_vcf=combine_deepvariant_vcf
    .combine(combine_deepvariant_vcf_tbi)
    .map{ meta, vcf, metaB, tbi ->
        [[
            id: meta.id,
            experimentalStrategy : meta.experimentalStrategy,
            genomeBuild : meta.genomeBuild,
            tumourNormalDesignation : meta.tumourNormalDesignation,
            sampleType : meta.sampleType ,
            gender : meta.gender,
            study_id : meta.study_id,
            num_intervals:  meta.num_intervals,
            tool : meta.tool
        ],vcf,tbi,[]]
    }

    //Filter for Indels
    SELECT_VCF_INDEL(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_INDEL.out.versions)


    //Filter for SNVs
    SELECT_VCF_SNV(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_SNV.out.versions)

    //Manipulate for payload ingestion
    ch_payload_vcf_indel=SELECT_VCF_INDEL.out.vcf.combine(SELECT_VCF_INDEL.out.tbi).combine(analysis_json)
            .map { metaA,vcf,metaB,tbi,analysis_json->
            [
                [
                    id : metaA.id,
                    experimentalStrategy : metaA.experimentalStrategy,
                    genomeBuild : metaA.genomeBuild,
                    tumourNormalDesignation : metaA.tumourNormalDesignation,
                    sampleType : metaA.sampleType ,
                    gender : metaA.gender,
                    study_id : metaA.study_id,
                    dataType : "InDel",
                    tool : metaA.tool
                ]
                ,[vcf, tbi],analysis_json]
            }

    ch_payload_vcf_snv=SELECT_VCF_SNV.out.vcf.combine(SELECT_VCF_SNV.out.tbi).combine(analysis_json)
            .map { metaA,vcf,metaB,tbi,analysis_json->
            [
                [
                    id : metaA.id,
                    experimentalStrategy : metaA.experimentalStrategy,
                    genomeBuild : metaA.genomeBuild,
                    tumourNormalDesignation : metaA.tumourNormalDesignation,
                    sampleType : metaA.sampleType ,
                    gender : metaA.gender,
                    study_id : metaA.study_id,
                    dataType : "SNV",
                    tool : metaA.tool
                ]
                ,[vcf, tbi],analysis_json]
            }

    //Generate payload
    PAYLOAD_VCF_INDEL_GERMLINEVARIANT(
        ch_payload_vcf_indel,
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        false
    )

    PAYLOAD_VCF_SNV_GERMLINEVARIANT(
        ch_payload_vcf_snv,
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        false
    )

    //Gather temporary files
    ch_cleanup=DEEPVARIANT.out.vcf.map{meta,vcf -> [vcf]}
    .mix(MERGE_DEEPVARIANT_VCF.out.vcf.map{meta,vcf -> [vcf]})
    .mix(TABIX_VC_DEEPVARIANT_VCF.out.tbi.map{meta,tbi -> [tbi]})
    .mix(SELECT_VCF_INDEL.out.vcf.map{meta,vcf -> [vcf]})
    .mix(SELECT_VCF_SNV.out.vcf.map{meta,vcf -> [vcf]})
    .mix(PAYLOAD_VCF_INDEL_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
    .mix(PAYLOAD_VCF_SNV_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
    .collect()

    //If Local is true, will be published into "output_dir" directory
    if (params.local==false){
        //Upload variants
        SONG_SCORE_UP_VCF_INDEL(PAYLOAD_VCF_INDEL_GERMLINEVARIANT)
        SONG_SCORE_UP_VCF_SNV(PAYLOAD_VCF_SNV_GERMLINEVARIANT)
        if (params.cleanup){
            CLEANUP(
                ch_cleanup.collect(),
                SONG_SCORE_UP_VCF_INDEL.out.analysis_id
                .combine(SONG_SCORE_UP_VCF_SNV.out.analysis_id).collect()
            )
        }
    } else {
     if (params.cleanup){
        CLEANUP(ch_cleanup.collect(),
        PAYLOAD_VCF_INDEL_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]}
        .combine(PAYLOAD_VCF_SNV_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
        )
        }       
    }


    emit:
    versions = ch_versions
}

