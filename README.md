# KG Vejle PacBio PureTarget pipeline


## General info:
This pipeline is used for PacBio PureTarget at Clinical Genetics, Vejle

## Default analysis steps and tools used

- Alignment (pbmm2)
- Repeat expansions incl. per-locus plots (TRGT) 
- Methylation profiles (pb-cpg-tools and methBat)

# Usage

## Default samplesheet format (use samplesheet used for sequencing)
The default samplesheet is simply the same as used to sequence the samples (.csv file)
The script will extract relevant information from the samplesheet for each sample, including:
- Samplename
- Biological material code
- LabWare testlist
- Gender


## Options and parameters:
    --help                  Show this help menu with available options
    
    --samplesheet   [path]: Path to samplesheet to use. Required
    
    --input         [path]: Path to data to use as input. 
                                Default: Not set. Instead, Search KG Vejle archive for input unmapped bams (search across all previous PacBio runs)

    ### Slurm Execution parameters:
    -profile slurm:         Run pipeline using KGVejle SLURM cluster
                                Default: Run pipeline on local server (where script is started)
    --slurmA        [bool]: Use secondary fast tmp storage (nfs_fast_a)
                                Default: Use primary fast tmp storage location at KGVejle


## Usage examples

#### Default: Analyze all samples in default samplesheet. Use all unmapped bam files available (across multiple SMRTcells) for each sample. Run all default analysis steps:
   
    nextflow run MadsAagaard/pureTarget -r main --samplesheet /path/to/samplesheet.csv
