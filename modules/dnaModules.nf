#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
*/
// The run header and the date/user/runID values moved to main.nf (entry workflow):
// strict syntax allows no statements at the top level of a module file.


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
    conda "${params.condaEnvs.pbmm2}"

    publishDir {"${params.outBase(meta)}/alignments/"}, mode: 'copy', pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genomeMmi} \
    ${data[0]} \
    ${meta.id}.${params.genomeVersion}.${params.readSet}.bam
    """
}

process pbmm2_align_mergedData {
    label "medium"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    publishDir {"${params.outBase(meta)}/alignments/"}, mode: 'copy', pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(fofn)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genomeMmi} \
    ${fofn} \
    ${meta.id}.${params.genomeVersion}.${params.readSet}.bam
    """
}



///////////////////////////////////////////////////
/////// ------- PSEUDO, VNTR AND REPEATS ------- //
///////////////////////////////////////////////////


process trgt4_pureTarget{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"

    publishDir {"${params.outBase(meta)}/TRGT/bam/"}, mode: 'copy', pattern: "*.trgt4.ba*"
    publishDir {"${params.outBase(meta)}/TRGT/vcf/"}, mode: 'copy', pattern: "*.trgt4.vcf.*"

/*
    publishDir {params.groupedOutput ? "${outputDir}/TRGT/pureTarget/bam" : "${outputDir}/${meta.id}/TRGT/pureTarget/bam/"}, mode: 'copy', pattern: "*.sorted.ba*"

    publishDir {params.groupedOutput ? "${outputDir}/TRGT/pureTarget/perSampleVcf/" : "${outputDir}/${meta.id}/TRGT/pureTarget/perSampleVcf/"}, mode: 'copy', pattern: "*.sorted.vcf.*"
*/


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam.bai"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz.tbi"),emit: trgt_full
    tuple val(meta),path ("*.trgt4.*")
    path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz"),emit: trgt_vcf
    script:

    //def sex=meta.sex=="male" ? "--karyotype XY" : "--karyotype XX"
    def karyotype=(meta.sex=="male"||meta.gender=="M") ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --preset targeted \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    $karyotype \
    --reads ${data[0]} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.STRchive

    bcftools sort -Oz -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz

    samtools sort -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.STRchive.spanning.bam
    samtools index ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam
    """
}

process trgt4_pureTarget_plots{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"
    publishDir {"${params.outBase(meta)}/TRGT/plots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    def geneList = params.puretargetPlotGenes.join(' ')
    """
    for gene in ${geneList}; do
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --squished \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.alleleSquished.pdf

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --plot-type waterfall \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.waterfall.pdf
    done
    """
}

process trgt4_pureTarget_plots_meth{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"

    publishDir {"${params.outBase(meta)}/TRGT/methylationPlots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    """
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --show meth \
    --squished \
    --max-allele-reads 75 \
    -o FXS_FMR1.${meta.id}.${params.genomeVersion}.${params.readSet}.METH.alleleSquished.pdf

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --plot-type waterfall \
    --show meth \
    --max-allele-reads 75 \
    -o FXS_FMR1.${meta.id}.${params.genomeVersion}.${params.readSet}.METH.waterfall.pdf

    """
}



process trgt5_pureTarget{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt5}"

    publishDir {"${params.outBase(meta)}/TRGT5/bam/"}, mode: 'copy', pattern: "*.trgt5.ba*"
    publishDir {"${params.outBase(meta)}/TRGT5/vcf/"}, mode: 'copy', pattern: "*.trgt5.vcf.*"


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam.bai"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz.tbi"),emit: trgt_full
    tuple val(meta),path ("*.trgt5.*")
    path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz"),emit: trgt_vcf
    script:

    //def sex=meta.sex=="male" ? "--karyotype XY" : "--karyotype XX"
    def karyotype=(meta.sex=="male"||meta.gender=="M") ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --preset targeted \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    $karyotype \
    --reads ${data[0]} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.STRchive

    bcftools sort -Oz -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz

    samtools sort -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.STRchive.spanning.bam
    samtools index ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam
    """
}

process trgt5_pureTarget_plots{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt5}"

    publishDir {"${params.outBase(meta)}/TRGT5/plots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"


    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    
    script:
    def geneList = params.puretargetPlotGenes.join(' ')
    """
    for gene in ${geneList}; do
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --squished \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.alleleSquished.pdf
     || echo "WARN: squished plot failed for \$gene" >&2

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --plot-type waterfall \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.waterfall.pdf
    || echo "WARN: waterfall plot failed for \$gene" >&2
    done
    """
}


///////////////////////////////////////////////////
/////// ------- METHYLATION ------- ///////////////
///////////////////////////////////////////////////
process methylationBW{
    tag "$meta.id"
    errorStrategy 'ignore'
    cpus 20
    conda "${params.condaEnvs.pbCPGtools}"

    publishDir {"${params.outBase(meta)}/methylation/BigWigBed/"}, mode: 'copy',  pattern: "*.methylation.{hap1,hap2,combined}.*"


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.methylation*")
    
    script:
    """
    aligned_bam_to_cpg_scores \
    --bam ${data[0]} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.readSet}.methylation
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
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.craminoQC.txt")

    script:
    """
    cramino \
    -t ${task.cpus} \
    ${data[0]} > ${meta.id}.${params.genomeVersion}.${params.readSet}.craminoQC.txt
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
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.nanostat.txt"),emit: multiqc
    path("${meta.id}.${params.genomeVersion}.${params.readSet}.nanostat.txt")
    script:
    """
    NanoStat \
    -t ${task.cpus} \
    -n ${meta.id}.${params.genomeVersion}.${params.readSet}.nanostat.txt \
    --bam ${data[0]}
    """
}
