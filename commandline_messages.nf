#!/usr/bin/env nextflow

ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"


/*
========================================================================================
    PRINT PARAMETER SUMMARY
========================================================================================
*/
def parameter_Message(){
    if(params.pipeline === 'bacass'){
        params.plas_assemblies = 'not used'
        params.plas_outdir = 'not used'
    } else if ((params.pipeline === 'plasmident') ||(params.pipeline === 'unicycler')){
        params.samplesheet_header = 'not used'
        params.genomesize = 'not used'
    }

    if (params.reads === 'short'){
        params.ont_reads_input = 'not used'
        params.ont_reads_outdir = 'not used'
    }else if (params.reads === 'long'){
        params.int_reads_input = 'not used'
    }
    log.info """
            ===========================================
            ${ANSI_GREEN}Bacterial hybrid assembly, Unicycler and PlasmIDent samplesheet generator${ANSI_RESET}

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

            Nf-core/bacass options:
            -------------------------------------------
            --samplesheet_header : ${params.samplesheet_header}
            --genomesize         : ${params.genomesize}
            
            plasmIDent options:
            -------------------------------------------
            --plas_assemblies    : ${params.plas_assemblies}
            --plas_outdir        : ${params.plas_outdir}

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
}

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/
def helpMessage() {
    log.info """
            ===========================================
            ${ANSI_GREEN}Bacterial hybrid assembly, Unicycler and PlasmIDent samplesheet generator${ANSI_RESET}

            This pipeline generates the samplesheets for the various options to perform a bacterial assembly.
            You can choose between the three different assembly types short, long and hybrid as well as between
            the two different assembly pipelines ${ANSI_GREEN}nf-core/bacass and unicycler${ANSI_RESET}.
            If long reads are specified a concatination of the individual FastQ files in the barcodeX directories will be performed.
            As extension the samplesheet for the ${ANSI_GREEN}imgag/plasmIDent${ANSI_RESET} pipeline can be generated.

            General parameters:
            -------------------------------------------
            --pipeline            : Pipeline for which a samplesheet should be generated | nf-core/bacass (default), unicycler or imgag/plasmIDent
            --reads               : Choice between short (default), long or hybrid assembly
            --mapping_file        : Path to tab-separated mapping file, has to contain the sample IDs for short reads and sample barcodes for long reads
                                    If you want to perform a hybrid assembly the first column has to contain the sample IDs and the second the long read barcodes
            --outdir              : The output directory were the results will be saved
            --multiqc_config      : Path to multiqc configuration file, ${ANSI_RED}currently not used${ANSI_RESET}

            Illumina options:
            -------------------------------------------
            --int_reads_input     : Path to the Illumina read directory, must contain the same name as the sample IDs in the mapping file
            
            Nanopore options:
            -------------------------------------------
            --ont_reads_input     : Path to the nanopore directories for the individual barcodes, these will be concatinated and used for the assembly.
                                    must contain the same name as the barcodes in the mapping file
            --ont_reads_outdir    : The output directory for the concatinated and renamed FastQ files for each individual barcode

            Nf-core/bacass options:
            -------------------------------------------
            --samplesheet_header : String with the header information of the samplesheets. Default is the bacass header since all other do not require a header
            --genomesize         : Expected genomesize in Megabases. Only used by the canu assembler in the nf-core/bacass pipeline. Default: 0
            
            plasmIDent options:
            -------------------------------------------
            --plas_assemblies    : The assembly directory. If no directory is given the pipeline will search for directories with the specified sample ID
            --plas_outdir        : The output directory for the plasmIDent pipeline, ${ANSI_RED}currently not used${ANSI_RESET}
            ===========================================
            """
            .stripIndent()

}

/*
========================================================================================
    PROCESS FINISHED MESSAGE
========================================================================================
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