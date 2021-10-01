#!/usr/bin/env nextflow

/*
========================================================================================
    bacass_samplesheet_generator
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

/*
========================================================================================
    CHECK PARAMETER INPUT
========================================================================================
*/

// check pipeline options
def samplesheet_header = ""
def genomesize = 0
def reads = params.reads
def pipeline = params.pipeline

if(pipeline == 'bacass'){
    samplesheet_header = params.samplesheet_header
    genomesize == params.genomesize
} else if (pipeline == 'unicycler'){
    samplesheet_header = null
} else{
    error "${ANSI_RED}Parameter pipeline: ${params.pipeline} is not one of the valid options: bacass or unicylcer${ANSI_RESET}"
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
    int_reads_input = "${params.int_reads_input}".replaceFirst(/$/, "/")
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
    PRINT PARAMETER SUMMARY
========================================================================================
*/
log.info """
        ===========================================
         ${ANSI_GREEN}Bacterial hybrid assembly samplesheet generator and PlasmIDent${ANSI_RESET}

         Used parameters:
        -------------------------------------------
         --pipeline            : ${params.pipeline}
         --reads               : ${params.reads}
         --mapping_file        : ${params.mapping_file}
         --outdir              : ${params.outdir}
         --multiqc_config      : ${params.multiqc_config}

         Illumina options:
        -------------------------------------------
         --int_reads_input     : ${params.int_reads_input}
        
        Nanopore options:
        -------------------------------------------
         --ont_reads_input     : ${params.ont_reads_input}
         --ont_reads_outdir    : ${params.ont_reads_outdir}
         --ont_dir_ids         : ${params.ont_reads_ids}

        Nf-core/bacass options:
          --samplesheet_header : ${params.samplesheet_header}
          --genomesize         : ${params.genomesize}

         Runtime data:
        -------------------------------------------
         Running with profile   : ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Used container         : ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
         Running as user        : ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir             : ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir               : ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         Nextflow version       : ${ANSI_GREEN}${nextflow.version}${ANSI_RESET}
        ===========================================
         """
         .stripIndent()

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/
def helpMessage() {
log.info """
        ===========================================
         ${ANSI_GREEN}I N T E R O P   and   B C L 2 F A S T Q   P I P E L I N E${ANSI_RESET}

         This pipeline takes an Illumina run output folder and runs the Illumina executables
         InterOp summary (sequencing run metrics) and bcl2fastq (converts bcl to fastq).
         By default, the file SampleSheet.csv from the runfolder is used in bcl2fastq, but another sample sheet file can be provided.
         The resulting fastq files are saved under ${ANSI_GREEN}results-bcl/fastq${ANSI_RESET}.
         The number of threads used by bcl2fastq are set to 4. If you know what you are doing,
         you can set them using the respective parameters.

         Usage:
        -------------------------------------------
         --runfolder            : Illumina run folder
         --outdir               : where results will be saved, default is "results-bcl"
         --samplesheet          : sample sheet file, default is runfolder/SampleSheet.csv
         --multiqc_config       : config file for MultiQC, default is "multiqc_config.yml"
         --title                : MultiQC report title, default is "InterOp and bcl2fastq summary"
         --load_threads         : Number of threads used for loading BCL data. 4 by default.
         --proc_threads         : Number of threads used for processing demultiplexed data. 4 by default.
         --write_threads        : number of threads used for writing FASTQ data. ${ANSI_RED}Must not be higher than number of samples!${ANSI_RESET} 4 by default.
         --barcode_mismatches:  : number of allowed barcode mismatches per index, 1 by default. Accepted values: 0,1,2.
        ===========================================
         """
         .stripIndent()

}


/*
========================================================================================
    WORKFLOW FOR PIPELINE
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
    python3 ${baseDir}/samplesheet_generation.py null $ont_reads_outdir ${launchDir}/${params.mapping_file} $reads $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.csv
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
    python3 ${baseDir}/samplesheet_generation.py $int_reads_input_ch $ont_reads_outdir ${launchDir}/${params.mapping_file} $reads $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.csv 
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
    """
    python3 ${baseDir}/samplesheet_generation.py $int_reads_input_ch null ${launchDir}/${params.mapping_file} $reads $pipeline $samplesheet_header $genomesize > samplesheet_${pipeline}.csv  
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
    
    if (reads == 'short'){
        GENERATE_SAMPLESHEET_WITH_SHORT(int_reads_input_ch)  
    } else if (reads == 'long'){
        MERGE_NANOPORE_DATA(ont_reads_input_ch, nanopore_ids)
        GENERATE_SAMPLESHEET_WITH_LONG()
    } else if (reads == 'hybrid'){
        MERGE_NANOPORE_DATA(ont_reads_input_ch, nanopore_ids)
        GENERATE_SAMPLESHEET_WITH_HYBRID(int_reads_input_ch)  

        if (pipeline == 'unicycler'){
            RUN_UNICYCLER(GENERATE_SAMPLESHEET_WITH_HYBRID.out, reads, MERGE_NANOPORE_DATA.out)
        }

    } else {
        error "${ANSI_RED}No process available for reads parameter: ${params.reads}${ANSI_RESET}"
    }

    

}

/*
 * Process finish line  ---------------------------------------------------------------------------------------------------------------------
 */
//=============================
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            ${ANSI_GREEN}Finished in ${workflow.duration}
            The samplesheet is here ==> ${ANSI_RESET}$params.outdir/samplesheet_${params.pipeline}.csv
            ===========================================
            """
            .stripIndent()
    }
    else {
        log.info """
            ===========================================
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}
