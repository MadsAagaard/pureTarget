#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
*/
date=new Date().format( 'yyMMdd' )
date2=new Date().format( 'yyMMdd HH:mm:ss' )
user="$USER"
runID="${date}.${user}"

//////////////////////////// SWITCHES ///////////////////////////////// 


log.info """\
======================================================
Clinical Genetics Vejle: PacBio PureTarget DEV
======================================================
Genome       : $params.genome
Genome FASTA : $params.genome_fasta
Genome MMI   : $params.genome_mmi
Genome ver.  : $params.genome_version
RunID        : $runID
Script start : $date2
ReadSet      : $params.readSet
outputDir    : $params.outputDirTMP 
"""


////////////////////////////////////////////
/////// ------- PREPROCESS + ALN ------- ///
////////////////////////////////////////////


process create_fofn {
    publishDir {"${params.outBase(meta)}/documents/"}, mode: 'copy',pattern: '*.fofn'

    input:
    tuple val(meta), path(data) //ubam

    output:
    tuple val(meta), path("${meta.id}.fofn")
    script:
    """
    `realpath ${data} > ${meta.id}.fofn`
    """
} 

process inputFiles_symlinks_ubam{
    label "low"
    publishDir {"${params.outBase(meta)}/inputSymlinks/"}, mode: 'symlink', pattern: '*.{bam,pbi}'

    input:
    tuple val(meta), path(data)   
    //data (default): 0:ubam, 1:ubam pbi

    output:
    tuple val(meta), path(data)

    script:
    """
    """
}

process pbmm2_align {
    errorStrategy 'ignore'
    tag "$meta.id"
    conda "${params.pbmm2}"

    publishDir {"${params.outBase(meta)}/alignments/"}, mode: 'copy', pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.bam"), path("${meta.id}.${params.genome_version}.${params.readSet}*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genome_mmi} \
    ${data[0]} \
    ${meta.id}.${params.genome_version}.${params.readSet}.bam
    """
}

process pbmm2_align_mergedData {
    label "medium"
    tag "$meta.id"
    conda "${params.pbmm2}"

    publishDir {"${params.outBase(meta)}/alignments/"}, mode: 'copy', pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(fofn)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.bam"), path("${meta.id}.${params.genome_version}.${params.readSet}*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genome_mmi} \
    ${fofn} \
    ${meta.id}.${params.genome_version}.${params.readSet}.bam
    """
}



///////////////////////////////////////////////////
/////// ------- PSEUDO, VNTR AND REPEATS ------- //
///////////////////////////////////////////////////


process trgt4_pureTarget{
    tag "$meta.id"
    label "medium"
    conda "${params.trgt4}"

    publishDir {"${params.outBase(meta)}/TRGT/bam/"}, mode: 'copy', pattern: "*.trgt4.ba*"
    publishDir {"${params.outBase(meta)}/TRGT/vcf/"}, mode: 'copy', pattern: "*.trgt4.vcf.*"

/*
    publishDir {params.groupedOutput ? "${outputDir}/TRGT/pureTarget/bam" : "${outputDir}/${meta.id}/TRGT/pureTarget/bam/"}, mode: 'copy', pattern: "*.sorted.ba*"

    publishDir {params.groupedOutput ? "${outputDir}/TRGT/pureTarget/perSampleVcf/" : "${outputDir}/${meta.id}/TRGT/pureTarget/perSampleVcf/"}, mode: 'copy', pattern: "*.sorted.vcf.*"
*/


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt4.bam"), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt4.bam.bai"), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt4.vcf.gz"), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt4.vcf.gz.tbi"),emit: trgt_full
    tuple val(meta),path ("*.trgt4.*")
    path("${meta.id}.${params.genome_version}.${params.readSet}.trgt4.vcf.gz"),emit: trgt_vcf
    script:

    //def sex=meta.sex=="male" ? "--karyotype XY" : "--karyotype XX"
    def karyotype=(meta.sex=="male"||meta.gender=="M") ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --preset targeted \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    $karyotype \
    --reads ${data[0]} \
    --output-prefix ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.STRchive

    bcftools sort -Oz -o ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.vcf.gz ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.vcf.gz

    samtools sort -o ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.bam ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.STRchive.spanning.bam
    samtools index ${meta.id}.${params.genome_version}.${params.readSet}.trgt4.bam
    """
}

process trgt4_pureTarget_plots{
    tag "$meta.id"
    label "medium"
    conda "${params.trgt4}"
    publishDir {"${params.outBase(meta)}/TRGT/plots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    """
    trgt plot \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id ${data.strID} \
    --squished \
    -o ${data.strID}.${meta.id}.${params.genome_version}.${params.readSet}.alleleSquished.pdf

    trgt plot \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id ${data.strID} \
    --plot-type waterfall \
    -o ${data.strID}.${meta.id}.${params.genome_version}.${params.readSet}.waterfall.pdf

    """
}

process trgt4_pureTarget_plots_meth{
    tag "$meta.id"
    label "medium"
    conda "${params.trgt4}"

    publishDir {"${params.outBase(meta)}/TRGT/methylationPlots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    """
    trgt plot \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --show meth \
    --squished \
    --max-allele-reads 75 \
    -o FXS_FMR1.${meta.id}.${params.genome_version}.${params.readSet}.METH.alleleSquished.pdf

    trgt plot \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --plot-type waterfall \
    --show meth \
    --max-allele-reads 75 \
    -o FXS_FMR1.${meta.id}.${params.genome_version}.${params.readSet}.METH.waterfall.pdf

    """
}



process trgt5_pureTarget{
    tag "$meta.id"
    label "medium"
    conda "${params.trgt5}"

    publishDir {"${params.outBase(meta)}/TRGT5/bam/"}, mode: 'copy', pattern: "*.trgt5.ba*"
    publishDir {"${params.outBase(meta)}/TRGT5/vcf/"}, mode: 'copy', pattern: "*.trgt5.vcf.*"


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt5.bam"), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt5.bam.bai"), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt5.vcf.gz"), path("${meta.id}.${params.genome_version}.${params.readSet}.trgt5.vcf.gz.tbi"),emit: trgt_full
    tuple val(meta),path ("*.trgt5.*")
    path("${meta.id}.${params.genome_version}.${params.readSet}.trgt5.vcf.gz"),emit: trgt_vcf
    script:

    //def sex=meta.sex=="male" ? "--karyotype XY" : "--karyotype XX"
    def karyotype=(meta.sex=="male"||meta.gender=="M") ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --preset targeted \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    $karyotype \
    --reads ${data[0]} \
    --output-prefix ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.STRchive

    bcftools sort -Oz -o ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.vcf.gz ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.vcf.gz

    samtools sort -o ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.bam ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.STRchive.spanning.bam
    samtools index ${meta.id}.${params.genome_version}.${params.readSet}.trgt5.bam
    """
}

process trgt5_pureTarget_plots{
    tag "$meta.id"
    label "medium"
    conda "${params.trgt5}"

    publishDir {"${params.outBase(meta)}/TRGT5/plots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"


    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    """
    trgt plot \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id ${data.strID} \
    --squished \
    -o ${data.strID}.${meta.id}.${params.genome_version}.${params.readSet}.alleleSquished.pdf

    trgt plot \
    --genome ${params.genome_fasta} \
    --repeats ${params.tr_pathogenic_v2} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id ${data.strID} \
    --plot-type waterfall \
    -o ${data.strID}.${meta.id}.${params.genome_version}.${params.readSet}.waterfall.pdf

    """
}


///////////////////////////////////////////////////
/////// ------- METHYLATION ------- ///////////////
///////////////////////////////////////////////////
process methylationBW{
    tag "$meta.id"
    errorStrategy 'ignore'
    cpus 20
    conda "${params.pbCPGtools}"

    publishDir {"${params.outBase(meta)}/methylation/BigWigBed/"}, mode: 'copy',  pattern: "*.methylation.{hap1,hap2,combined}.*"


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.methylation*")
    
    script:
    """
    aligned_bam_to_cpg_scores \
    --bam ${data[0]} \
    --output-prefix ${meta.id}.${params.genome_version}.${params.readSet}.methylation
    """
}

process methylationSegm{
    tag "$meta.id"
    errorStrategy 'ignore'
    cpus 20
    conda "${params.methbat}"

    publishDir {"${params.outBase(meta)}/methylation/"}, mode: 'copy'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.met.*")
    
    script:
    """
    methbat profile \
    --input-prefix ${meta.id}.methylation \
    --input-regions ${methylationBackground} \
    --output-region-profile ${meta.id}.met.profile
    """
}


///////////////////////////////////////////////////
/////// ------- QUALITY CONTROL ------- ///////////
///////////////////////////////////////////////////


process qualimap {
    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 10
    maxForks 6

    //publishDir {params.groupedOutput ? "${outputDir}/QC/qualimap/" : "${outputDir}/${meta.id}/QC/qualimap/"}, mode: 'copy'


    conda '/lnx01_data3/shared/programmer/miniconda3/envs/qualimapSamtools/' 
    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}_qualimap"), emit: qualimap_out
    path("${meta.id}_qualimap"), emit: multiqc

    //path ("versions.yml"), emit: versions

    script:
    use_bed = ''//qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    qualimap --java-mem-size=40G bamqc \
    -nt ${task.cpus} \
    -outdir ${meta.id}_qualimap \
    -bam ${data[0]} 

    """
}

process mosdepthROI {
    tag "$meta.id"
    errorStrategy 'ignore'

    //publishDir {params.groupedOutput ? "${outputDir}/QC/mosdepth/" : "${outputDir}/${meta.id}/QC/mosdepth/"}, mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/mosdepth/' 
    cpus 8

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}_roi.*"),emit: mosdepth_roi
    tuple val(meta), path("*.region.dist.txt"), emit:multiqc
    script:
    def callable=params.genome=="hg38" ? "--by ${CALLABLE_ROI}" : "--by 1000"
    """
    mosdepth \
    -t ${task.cpus} \
    $callable \
    ${meta.id}.${params.genome_version}_roi \
    ${data[0]}

    """
}


process whatsHap_stats {
    tag "$meta.id"
    errorStrategy 'ignore'

    //publishDir {params.groupedOutput ? "${outputDir}/QC/whatsHap/" : "${outputDir}/${meta.id}/QC/whatsHap/"}, mode: 'copy'


    conda '/lnx01_data3/shared/programmer/miniconda3/envs/whatshap/' 
    cpus 8

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.whatshap.stats.tsv"),emit:multiqc

    script:
    """
    whatshap stats \
    ${data[0]} \
    --tsv=${meta.id}.whatshap.stats.tsv
    """
}

process cramino {
    tag "$meta.id"
    errorStrategy 'ignore'

    //publishDir {params.groupedOutput ? "${outputDir}/QC/cramino/" : "${outputDir}/${meta.id}/QC/cramino/"}, mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/cramino/' 
    cpus 8

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.craminoQC.txt")

    script:
    """
    cramino \
    -t ${task.cpus} \
    ${data[0]} > ${meta.id}.${params.genome_version}.${params.readSet}.craminoQC.txt
    """
}

process nanoStat {
    tag "$meta.id"
    errorStrategy 'ignore'
   
    //publishDir {params.groupedOutput ? "${outputDir}/QC/nanoStat/" : "${outputDir}/${meta.id}/QC/nanoStat/"}, mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/nanostats/' 
    cpus 8

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.id}.${params.genome_version}.${params.readSet}.nanostat.txt"),emit: multiqc
    path("${meta.id}.${params.genome_version}.${params.readSet}.nanostat.txt")
    script:
    """
    NanoStat \
    -t ${task.cpus} \
    -n ${meta.id}.${params.genome_version}.${params.readSet}.nanostat.txt \
    --bam ${data[0]}
    """
}


process multiQC {
    tag "$meta.id"
    errorStrategy 'ignore'

    //publishDir {params.groupedOutput ? "${outputDir}/QC/" : "${outputDir}/${meta.id}/QC/"}, mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/multiqc/' 

    when:
    !params.groupedOutput

    input:
    tuple val(meta),  path(data)  

    output:
    path ("*MultiQC*.html")

    script:
    """
    multiqc \
    -c ${multiqc_config} \
    -f -q ${launchDir}/${outputDir}/${meta.id}/QC/ \
    -n ${meta.id}.MultiQC.DNA.html
    """
}

process multiQC_ALL {
    
    //errorStrategy 'ignore'
    //publishDir "${outputDir}/", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/multiqc/' 

    input:
    tuple val(meta),path(data)  
    //   path("_fastqc.*").collect().ifEmpty([])
    // path("${meta.id}.samtools.sample.stats.txt").collect().ifEmpty([])
    // path("bamQC/*").collect().ifEmpty([]) 
    //path("${meta.id}.picardWGSmetrics.txt").collect().ifEmpty([]) 

    output:
    path ("${params.rundir}.MultiQC.ALL.html")

    script:
    """
    multiqc \
    -c ${multiqc_config} \
    -f -q ${launchDir}/${outputDir}/*/QC/ \
    -n ${params.rundir}.MultiQC.ALL.html
    """
}

/////////////// TO DO /////////////////////

process vntyper2 {
    errorStrategy 'ignore'
    //publishDir "${outputDir}/MUC1-VNTR_kestrel/", mode: 'copy'
    cpus 16

    input:
    tuple val(meta), path(reads)

    output:
    //tuple val(meta), path("vntyper${meta.id}.vntyper/*")
    tuple val(meta), path("*/*.{tsv,vcf}")
    script:
    
    def reads_command = "--fastq1 ${reads[0]} --fastq2 ${reads[1]}"
    
    """
    singularity run -B ${s_bind} ${simgpath}/vntyper20.sif \
    -ref ${vntyperREF}/chr1.fa \
    --fastq1 ${r1} --fastq2 ${r2} \
    -t ${task.cpus} \
    -w vntyper \
    -m ${vntyperREF}/hg19_genic_VNTRs.db \
    -o ${meta.id} \
    -ref_VNTR ${vntyperREF}/MUC1-VNTR_NEW.fa \
    --fastq \
    --ignore_advntr \
    -p /data/shared/programmer/vntyper/VNtyper/
    """
}


process advntr {

    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 8
    //publishDir "${outputDir}/${meta.id}/advntr/", mode: 'copy', pattern: "*.advntr.*"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/advntr15/'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.advntr.*")
    
    script:
    """
    advntr genotype \
    -f ${data[0]} \
    --pacbio \
    -m ${vntr_defaultModel} \
    -o ${meta.id}.advntrDefault
    """
}



