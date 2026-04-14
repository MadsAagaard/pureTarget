#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


//////////// DEFAULT INPUT ///////////////////////

def inputError() {
    log.info"""
    USER INPUT ERROR: The user should point to a samplesheet (--samplesheet parameter) or input folder containing all data to be used as input (--input parameter).
    """.stripIndent()
}

if (!params.samplesheet && !params.input) exit 0, inputError() 

if (params.samplesheet) {


}

if (params.input) {
    if (params.hifiReads){
        inputBam="${params.input}/**/*.hifi_reads.*.bam"
    }
    if (!params.hifiReads){
        inputBam="${params.input}/*_reads/*.bam"
    }
    }

if (!params.input && (params.samplesheet)) {
    if (params.hifiReads){
        inputBam="${params.dataArchive}/**/*.hifi_reads.*.bam"
    }
    if (!params.hifiReads){
        inputBam="${params.dataArchive}/**/*.bam"
    }
    }

if (params.samplesheet) {

    // Input data (ubam) channel
    Channel.fromPath(inputBam, followLinks: true)
    |map { tuple(it.baseName,it) }

    |map {id,bam -> 
        def (samplenameFull,pacbioID,readset,barcode)   =id.tokenize(".")
        def (samplename,material,testlist,gender)       =samplenameFull.tokenize("_")
            meta=[npn:samplename]
            tuple(meta,bam)        
        }
    |groupTuple(sort:true)
    |branch  {meta,bam -> 
        UNASSIGNED: (meta.npn=~/UNASSIGNED/)
                    return [meta,bam]
        samples: true
                    return [meta,bam]
    }
    | set {ubam_input }


    ubam_input.samples
        | map { meta, bam -> tuple(meta.npn,meta,bam) }
        | set {ubam_input_samples_join}    


    // Samplesheet channel
    Channel.fromPath(params.samplesheet, checkIfExists: true)
    .flatMap { file ->

    def metaRunID = file.baseName.replaceFirst(/_\d+.*/, '')

    def lines = file.text.readLines()
      .collect { it.trim() }
      .findAll { it }

    def header = 'Bio Sample Name,Plate Well,Adapter,Adapter2'
    def idx = lines.findIndexOf { it == header }
    if( idx < 0 )
      return []

    def sampleLines = lines[(idx+1)..<lines.size()]

    sampleLines.collect { L ->
      def cols = L.split(/\s*,\s*/, -1)
      tuple(metaRunID, cols[0], cols[1], cols[2], cols[3])
    }
    }
    .map { metaRunID, sampleInfo, flowcell, barcode1, barcode2 ->

    def (npn, material, testlist, gender) = sampleInfo.tokenize("_")
    def sex = (gender == "K") ? "female" : "male"

    def meta = [
        npn             : npn,
        testlist        : testlist,
        sex             : sex,
        id              : "${npn}.${testlist}",
        metaRunID       : metaRunID,
        analysisDate    : date,
        outKey          : metaRunID
    ]
    meta
    }
    |set { samplesheet_full }

    samplesheet_full
    |map {row -> meta2=[row.npn,row]}
    |set {samplesheet_join}


    //join samplesheet and bam channels:
    samplesheet_join.join(ubam_input_samples_join)
    |map { samplename, metaSS, metaData, bam ->
        def merged = [
            // Keys from metaSS are authoritative
            npn          : metaSS.npn,
            testlist     : metaSS.testlist,
            sex          : metaSS.sex,
            id           : metaSS.id,
            metaRunID    : metaSS.metaRunID,
            analysisDate : metaSS.analysisDate,
            outKey       : metaSS.outKey,
            // Keys from metaData that add useful info not in metaSS
            rundate      : metaData.rundate,
            gender       : metaData.gender,
        ]
        tuple(merged, bam)
    }
    |set {finalUbamInput}
}

if (!params.samplesheet && params.input) {
    def inputRunID = file(params.input).name
    Channel.fromPath(inputBam, followLinks: true)
    |map { tuple(it.baseName, it) }
    |map { id, bam ->
        def (samplenameFull, pacbioID, readset, barcode) = id.tokenize(".")
        def (instrument, date2, time)                    = pacbioID.tokenize("_")
        def (samplename, material, testlist, gender)     = samplenameFull.tokenize("_")
        def sex = (gender == "K") ? "female" : "male"   // same logic as SS branch
        def meta = [
            npn          : samplename,                   
            id           : "${samplename}.${testlist}", 
            testlist     : testlist,
            sex          : sex,                         
            gender       : gender,                      
            rundate      : date2,
            metaRunID    : inputRunID,
            analysisDate : date,
            outKey       : inputRunID     
        ]
        tuple(meta, bam)
    }
    |groupTuple(sort: true)
    |branch { meta, bam ->
        UNASSIGNED: (meta.npn =~ /UNASSIGNED/)
            return [meta, bam]
        samples: true
            return [meta, bam]
    }
    |set { ubam_input }
    
    ubam_input.samples
    |set {finalUbamInput}
}


/////////////////// MODULES ///////////////////////
include {create_fofn;
        pbmm2_align_mergedData;
        inputFiles_symlinks_ubam;
        trgt4_pureTarget;
        trgt5_pureTarget;
        methylationBW;
        methylationSegm;
        trgt4_pureTarget_plots;
        trgt4_pureTarget_plots_meth;
        trgt5_pureTarget_plots;
        } from "./modules/dnaModules.nf" 

////////////////// WORKFLOWS AND PROCESSES ///////////////////////

workflow PREPROCESS {

    take:
    finalUbamInput     
   
    main:

    inputFiles_symlinks_ubam(finalUbamInput)

    create_fofn(finalUbamInput)
    pbmm2_align_mergedData(create_fofn.out)

    emit:
    aligned=pbmm2_align_mergedData.out.bam
    
}

workflow {



    PREPROCESS(finalUbamInput)

    PREPROCESS.out.aligned
    | map {meta,bam,bai -> tuple(meta,[bam,bai])}
    |set {alignedFinal}

    trgt4_pureTarget(alignedFinal)
    trgt5_pureTarget(alignedFinal)

    trgt4_pureTarget.out.trgt_full.combine(params.puretargetPlotGenes)
    |map {meta,bam,bai,vcf,tbi,genes -> 
    tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi,strID:genes])}
    //tuple(meta,bam,genes)}
    |set {trgt4_plot_ch}
    trgt4_pureTarget_plots(trgt4_plot_ch)


    trgt4_pureTarget.out.trgt_full
    |map {meta,bam,bai,vcf,tbi -> 
    tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi])}
    |set {trgt4_plot_ch_meth}

    trgt4_pureTarget_plots_meth(trgt4_plot_ch_meth)

    trgt5_pureTarget.out.trgt_full.combine(params.puretargetPlotGenes)
    |map {meta,bam,bai,vcf,tbi,genes -> 
    tuple(meta,[bam:bam,bai:bai,vcf:vcf,tbi:tbi,strID:genes])}
    |set {trgt5_plot_ch}
    trgt5_pureTarget_plots(trgt5_plot_ch)

 methylationBW(alignedFinal)

}
      