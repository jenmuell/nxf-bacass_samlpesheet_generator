#!/usr/bin/env nextflow

/*
========================================================================================
    nxf-bacass_samplesheet_generator
========================================================================================
    Github : 
    Website: 
    Slack  : 
----------------------------------------------------------------------------------------
*/

/*
========================================================================================
    GENERAL SETTINGS 
========================================================================================
*/

nextflow.enable.dsl = 2
include { helpMessage; parameter_Message } from './commandline_messages'

/*
* ANSI escape codes to color output messages
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

// exit early if help
params.help = ""
if (params.help) {
    helpMessage()
    exit(0)
}

parameter_Message()

/*
========================================================================================
    CHECK PARAMETER INPUT
========================================================================================
*/

// check pipeline options
def samplesheet_header = ""
def genomesize = 0f
def reads = params.reads
def pipeline = params.pipeline

if(pipeline == 'bacass'){
    samplesheet_header = params.samplesheet_header
    genomesize = params.genomesize
} else if (pipeline == 'unicycler'){
    samplesheet_header = null
} else if (pipeline == 'plasmident'){
    samplesheet_header = null
} else if (pipeline == 'viralrecon'){
    samplesheet_header = params.samplesheet_header
}else{
    error "${ANSI_RED}Parameter pipeline: ${params.pipeline} is not one of the valid options: bacass, viralrecon, unicylcer or plasmident${ANSI_RESET}"
}


// check reads options & mapping file
// Read inputs
def int_reads_input = null
def ont_reads_input = null
def ont_reads_outdir = null
def ont_reads_ids = null

// Mapping file parameters
def mapping_file = null
def illumina_ids_col = []
def nanopore_ids_col = []
def illumina_ids = null
def nanopore_ids = null

if(reads == 'short'){
    int_reads_input = "${params.int_reads_input}".replaceFirst(/$/, "/")
    int_reads_input_ch = Channel.fromPath(int_reads_input, checkIfExists: true, type: 'dir')
        .ifEmpty { error "Can not find folder ${int_reads_input}" }

    // Mapping file & Read IDs
    mapping_file = file(params.mapping_file)
                    .splitEachLine(/\t/){list -> 
                        illumina_ids_col.add(list[0])
                    }


    // Value channel with Illumina read IDs
    if (illumina_ids_col.size() > 0){
        illumina_ids = Channel.from(illumina_ids_col)
    } else {
        error "${ANSI_RED}No Illumina IDs were defined. Check your mapping file ${params.mapping_file}${ANSI_RESET}"
    }

} else if (reads == 'long'){
    ont_reads_input = "${params.ont_reads_input}".replaceFirst(/$/, "/")
    ont_reads_input_ch = Channel.fromPath(ont_reads_input, checkIfExists: true, type: 'dir')
        .ifEmpty { error "Can not find folder ${ont_reads_input}" }

    ont_reads_outdir = "${params.ont_reads_outdir}".replaceFirst(/$/, "/")
    
    // Mapping file & Read IDs 
    mapping_file = file(params.mapping_file)
                    .splitEachLine(/\t/){list -> 
                        nanopore_ids_col.add(list[0])
                    }

    // Value channel with Nanopore read IDs
    if (nanopore_ids_col.size() > 0){
        nanopore_ids = Channel.from(nanopore_ids_col)
    } else {
        error "${ANSI_RED}No Nanpore IDs were defined. Check your mapping file ${params.mapping_file}${ANSI_RESET}"
    }
    
   
} else if (reads == 'hybrid'){
    // Read inputs
    int_reads_input = "${params.int_reads_input}"
    int_reads_input_ch = Channel.fromPath(int_reads_input, checkIfExists: true, type: 'dir')
        .ifEmpty { error "Can not find folder ${int_reads_input}" }

    ont_reads_input = "${params.ont_reads_input}".replaceFirst(/$/, "/")
    ont_reads_input_ch = Channel.fromPath(ont_reads_input, checkIfExists: true, type: 'dir')
        .ifEmpty { error "Can not find folder ${ont_reads_input}" }

    ont_reads_outdir = "${params.ont_reads_outdir}".replaceFirst(/$/, "/")

    ont_reads_ids = params.ont_reads_ids
    ont_reads_ids_ch = Channel.from(ont_reads_ids)

    // Mapping file & Read IDs 
    mapping_file = file(params.mapping_file)
                    .splitEachLine(/\t/){list -> 
                        illumina_ids_col.add(list[0])
                        nanopore_ids_col.add(list[1])
                    }

    // Value channel with Illumina read IDs
    if (illumina_ids_col.size() > 0){
        illumina_ids = Channel.from(illumina_ids_col)
    } else {
        error "${ANSI_RED}No Illumina IDs were defined. Check your mapping file ${params.mapping_file}${ANSI_RESET}"
    }

    // Value channel with Nanopore read IDs
    if (nanopore_ids_col.size() > 0){
        nanopore_ids = Channel.from(nanopore_ids_col)
    } else {
        error "${ANSI_RED}No Nanpore IDs were defined. Check your mapping file ${params.mapping_file}${ANSI_RESET}"
    }

} else{
    error "${ANSI_RED}Parameter reads: ${params.reads} is  not one of the options: short, long, hybrid${ANSI_RESET}"
}


mqc_config = file(params.multiqc_config) // this is needed, otherwise the multiqc config file is not available in the docker image
outdir = file(params.outdir)


/*
========================================================================================
    WORKFLOWS FOR PREPROCESSING AND SAMPLESHEET GENERATION
========================================================================================
*/
process MERGE_NANOPORE_DATA {
    tag "merge nanopore data from ${ont_reads_input_ch}"
    publishDir "${params.ont_reads_outdir}", mode: 'copy'
    
    input:
        path ont_reads_input_ch
        each nanopore
    
    output:
        file 'barcode*.fastq'
        

    script:
    """
    echo cat ${ont_reads_input_ch}/${nanopore}/PAE*.fastq > ${nanopore}.fastq
    cat ${ont_reads_input_ch}/${nanopore}/PAE*.fastq > ${nanopore}.fastq
    """

}

process GENERATE_SAMPLESHEET_WITH_LONG {
    tag "generate samplesheet for long read assembly"
    publishDir params.outdir, mode: 'copy'


    output:
        file "samplesheet_${pipeline}.csv"
    
    script:
    """
    python3 ${baseDir}/samplesheet_generation.py null $ont_reads_outdir null ${launchDir}/${params.mapping_file} $reads $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.csv
    """
}

process GENERATE_SAMPLESHEET_WITH_HYBRID {
    tag "generate samplesheet for hybrid assembly"
    publishDir params.outdir, mode: 'copy'

    input:  
        path int_reads_input_ch                 // all Illumina reads

    output:
        file "samplesheet_${pipeline}.csv"
    
    script:
    """
    python3 ${baseDir}/samplesheet_generation.py $int_reads_input_ch $ont_reads_outdir null ${launchDir}/${params.mapping_file} $reads $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.csv 
    """
}

process GENERATE_SAMPLESHEET_WITH_SHORT {
    tag "generate samplesheet for short read assembly"
    publishDir params.outdir, mode: 'copy'

    input:
        path int_reads_input_ch
    
    output:
        file "samplesheet_${pipeline}.csv"

    script:
    int_dir = int_reads_input_ch.toString()
    """
    python3 ${baseDir}/samplesheet_generation.py ${launchDir}/${int_reads_input} null null ${launchDir}/${params.mapping_file} $reads $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.csv  
    """
}

process GENERATE_SAMPLESHEET_PLASMIDENT {
    tag "generate samplesheet for plasmIDent pipeline"
    publishDir params.outdir, mode: 'copy'

    output:
        file "samplesheet_${pipeline}.tsv"
    
    
    script:
    """
    python3 ${baseDir}/samplesheet_generation.py null $ont_reads_outdir ${params.plas_assemblies} ${launchDir}/${params.mapping_file} null $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.tsv
    """
}



/*
========================================================================================
    WORKFLOWS FOR ANALYSIS WITH PLASMIDENT AND UNICYCLER
========================================================================================
*/
// TODO: does not work, problem with command
process RUN_PLASMIDENT {
    tag "run plasmIDent pipeline"

    input:
        file samplesheet_plasmident
    
    script:
    """
    sudo nextflow run imgag/plasmIDent --input $samplesheet_plasmident --outDir ${params.plas_outdir} -profile docker
    """
}

// TODO: does not work, problem with conda enviroments
process RUN_UNICYCLER {
    tag "run unicycler"
    publishDir "${params.outdir}/results_unicycler", mode: 'copy'
    conda "unicycler"

    input: 
        file samplesheet 
        val reads_process
        file ont_files
    
    output:
        path "results_unicycler/" 
    
    script:
    if (reads_process == 'short')
    """
    conda activate unicycler
    while IFS="," read sample read1 read2; \
     do unicycler -1 \${read1} -2 \${read2} -o \${sample}; \
     done < ${samplesheet}
    conda deactivate
    """
    else if (reads_process == 'long')
    """
    conda activate unicycler
    while IFS="," read sample longread; \
     do unicycler -l \${longread} -o \${sample}; \
     done < ${samplesheet}
    conda deactivate
    """
    else if (reads_process == 'hybrid')
    """
    conda activate unicycler
    while IFS="," read sample read1 read2 longread; \
     do unicycler -1 \${read1} -2 \${read2} -l \${longread} -o \${sample}; \
     done < ${samplesheet}
    conda deactivate
    """
    else 
        error "${ANSI_RED}Unicycler run error: ${reads_process} is not one of the valid options: short, long, hybrid${ANSI_RESET}"
}


// Need an additional process for plasmident and nf-core/viralrecon pipeline
workflow {
    
    if (pipeline == 'plasmident'){
        GENERATE_SAMPLESHEET_PLASMIDENT()
        //RUN_PLASMIDENT(GENERATE_SAMPLESHEET_PLASMIDENT.out)
    } else if (pipeline == 'viralrecon'){
        GENERATE_SAMPLESHEET_WITH_SHORT(int_reads_input_ch)
    }else if (reads == 'short'){
        GENERATE_SAMPLESHEET_WITH_SHORT(int_reads_input_ch)  
    } else if (reads == 'long'){
        MERGE_NANOPORE_DATA(ont_reads_input_ch, nanopore_ids)
        GENERATE_SAMPLESHEET_WITH_LONG()
    } else if (reads == 'hybrid'){
        MERGE_NANOPORE_DATA(ont_reads_input_ch, nanopore_ids)
        GENERATE_SAMPLESHEET_WITH_HYBRID(int_reads_input_ch)  

        //if (pipeline == 'unicycler'){
          //  RUN_UNICYCLER(GENERATE_SAMPLESHEET_WITH_HYBRID.out, reads, MERGE_NANOPORE_DATA.out)
        //}

    }else {
        error "${ANSI_RED}No process available for reads parameter: ${params.reads}${ANSI_RESET}"
    }

    

}


