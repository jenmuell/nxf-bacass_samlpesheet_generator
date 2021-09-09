// Interop - bcl2fastq - MultiQC pipeline

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
 * pipeline input parameters  ---------------------------------------------------------------------------------------------------------------------
 */
// use bcl test data from Illumina basespace
 
// get docker ncores, not nextflow env ncores!!				
 int ncores       = Runtime.getRuntime().availableProcessors(); 
 int usable_cores = Math.floor(ncores*0.8).toInteger()				// JM: Why is this code used? --> garantee a 80 percent allocation of the cpu

 params.runfolder   = ""							// JM: parameters which can be changed
 params.samplesheet = "${params.runfolder}/SampleSheet.csv"
 
 // find number of samples in order to dynamically set write_threads, must NOT be higher than number of samples
 def sample_index = file(params.samplesheet)
        .readLines()
        .findLastIndexOf {it =~ "Sample_ID"}                                    // JM: get number of lines until header of Samples 
 
 def total_lines = file(params.samplesheet)					// JM: all lines in the files
        .countLines()
 
 def nsamples = total_lines - (sample_index + 1)				
 //println "total_lines: $total_lines, sample_index: $sample_index, nsamples: $nsamples"
 
 params.outdir = "$workflow.launchDir/results-bcl"				// JM: rest of the parameters 
 params.title = "InterOp and bcl2fastq summary"
 params.multiqc_config = "$baseDir/multiqc_config.yml" //in case ncct multiqc config needed   // JM: can be replaced against multiqc_config_ncct.yml example in nxf-fastqc --> MS in generell not used
 params.load_threads = usable_cores
 params.proc_threads = usable_cores
 params.barcode_mismatches = 1 // default of bcl2fastq. Allows 0,1,2, so set limits also here  // JM: allow (2n+1) mismatches 
 params.scratch = false // used in special cases, stages the bcl process in a local dir

 if (nsamples >= usable_cores) {						// JM: number of threads can not exceed number of samples
     params.write_threads = usable_cores
     } else {
     params.write_threads = nsamples
         } //must not be higher than number of samples

 mqc_config = file(params.multiqc_config) // this is needed, otherwise the multiqc config file is not available in the docker image


/*
 * Standard output per run  ---------------------------------------------------------------------------------------------------------------------
 */
log.info """
        ===========================================
         ${ANSI_GREEN}I N T E R O P   and   B C L 2 F A S T Q   P I P E L I N E${ANSI_RESET}

         Used parameters:
        -------------------------------------------
         --runfolder            : ${params.runfolder}
         --samplesheet          : ${params.samplesheet}
         --outdir               : ${params.outdir}
         --multiqc_config       : ${params.multiqc_config}
         --title                : ${params.title}
         --load_threads         : ${params.load_threads}
         --proc_threads         : ${params.proc_threads}
         --write_threads        : ${params.write_threads}
         --barcode_mismatches   : ${params.barcode_mismatches}

         Runtime data:
        -------------------------------------------
         Running with profile   : ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Used container         : ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
         Running as user        : ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir             : ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir               : ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         Number of host cores   : ${ANSI_GREEN}${ncores}${ANSI_RESET}
         Estimated # of samples : ${ANSI_GREEN}${nsamples}${ANSI_RESET}
         Nextflow version       : ${ANSI_GREEN}${nextflow.version}${ANSI_RESET}
        ===========================================
         """
         .stripIndent()

/*
 * Help message  ---------------------------------------------------------------------------------------------------------------------
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
 * Start processes and parameter testing  ---------------------------------------------------------------------------------------------------------------------
 */
// in case trailing slash in runfolder not provided:
runfolder_repaired = "${params.runfolder}".replaceFirst(/$/, "/")

Channel										// JM: channel for runfolder connected to input of InterOP and bcl2fastq
    .fromPath(runfolder_repaired, checkIfExists: true, type: 'dir')
    .ifEmpty { error "Can not find folder ${runfolder_repaired}" }
    .into {runfolder_interop; runfolder_bcl}

Channel										// JM: channel for samplesheet connected to input of bcl2fastq
    .fromPath(params.samplesheet, checkIfExists: true, type: 'file')		//     Is there a Error message if the samplesheet is not definied?
    .set {samplesheet_ch}

//runfolder_ch.into {runfolder_interop; runfolder_bcl}


process interop {
    tag "interop on $x"								// JM: label to find the process in the log file

    input:
        path x from runfolder_interop

    // path is preferred over file as an output qualifier
    output:
        path 'interop_file' into interop_ch

    script:
    """
    interop_summary --csv=1 $x > interop_file
    """
    // JM: --csv=1 to be compatible with MultiQC 
}

// I want to print the ncores and nsamples to adjust the bcl2fastq parameters in subsequent runs
process bcl {

    tag "bcl2fastq"
    publishDir params.outdir, mode: 'copy', pattern: 'fastq/**fastq.gz'   // JM: path to the fastq files, copy the files from the pipeline work directory
    //publishDir params.outdir, mode: 'copy', pattern: 'fastq/Stats/*'
    publishDir params.outdir, mode: 'copy', pattern: 'bcl_out.log'        // JM: log file if the run was successful
    scratch params.scratch
    /* I use this scratch to be able to tail -f the bcl_out file, to monitor progress in the shiny app
    echo works but with delay, while this is real time*/
    
    input:
        path x from runfolder_bcl
        path y from samplesheet_ch

    // path is preferred over file as an output qualifier
    
    output:
        path 'fastq/Stats/Stats.json' into bcl_ch
        path 'fastq/**fastq.gz' // ** is for recursive match, directories are omitted (the fastq files might be in fastq/someproject/...)
        path 'bcl_out.log'
    
    // default to --ignore-missing all?
    script:
    
    """
    bcl2fastq -R $x \
    -o fastq \
    --sample-sheet $y \
    --no-lane-splitting \
    --ignore-missing-bcls \
    --barcode-mismatches ${params.barcode_mismatches} \
    -r ${params.load_threads} \
    -p ${params.proc_threads} \
    -w ${params.write_threads} >bcl_out.log 2>&1
    """
    /* JM: parameter description
    * -R: location of the run folder
    * -o: location of the demultiplexed output (fastq files)
    * --sample-sheet: path to the sample sheet
    * --no-lane-splitting: do not split FASTQ files lane
    * --barcode-mismatches: allowed mismatches per index adapter (0,1,2) default=0
    * (--ignore-missing): available for bcls, filter, positions --> ignore the missing and false information
    */
}

process multiqc {
    publishDir params.outdir, mode: 'copy'

    input:
        file interop_file from interop_ch		// JM: logs from InterOP
        file 'Stats.json' from bcl_ch			// logs from bcl2fastq
        file mqc_config					// configuration file for MultiQC

    output:
        path 'multiqc_report.html'
        path 'multiqc_report_data/multiqc_general_stats.txt' into excel_ch
        path 'multiqc_report_data' type 'dir'

    script:
    """
    multiqc --force --interactive \
    --title "${params.title}" \
    --config $mqc_config \
    --filename "multiqc_report.html" \
    $interop_file Stats.json
    """
    /* JM: parameter description
     * --force: overwrite any existing reports
     * (--exclude): exclude this modules
     */
}

// take the multiqc_general_stats.txt and write an excel file from it
process excel {
    publishDir params.outdir, mode: 'copy'

    input:
        path mqc_stats from excel_ch
    
    output:
        path 'multiqc_general_stats.xlsx'
    
    script:
    """
    write_excel.R $mqc_stats
    """

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
            See the report here      ==> ${ANSI_RESET}$params.outdir/multiqc_report.html${ANSI_GREEN}
            The fastq files are here ==> ${ANSI_RESET}$params.outdir/fastq/
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
