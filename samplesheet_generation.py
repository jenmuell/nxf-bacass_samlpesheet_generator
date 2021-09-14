#!/usr/bin/python

from os import listdir
from os.path import isfile, join
import sys

def generate_samplesheet_short(illumina_files, mapping_file, pipeline, genomesize):
    samplesheet_string = ''
    illumina_read_files = [f for f in listdir(illumina_files) if isfile(join(illumina_files,f))]           # list of all illumina reads

    for line in mapping_file:  # Read mapping file line by line
        illumina_reads_id = [file for file in illumina_read_files if file.contain(line)]  # R1 and R1 for sample 
        if pipeline == 'bacass':
            samplesheet_string += str(line) + '\t' + illumina_reads_id[0] + '\t' + illumina_reads_id[1] + '\tNA\tNA\t' + str(genomesize) + '\n'
    
        elif pipeline == 'unicycler':
            samplesheet_string += str(line) + '\t' + illumina_reads_id[0] + '\t' + illumina_reads_id[1] + '\tNA\n'
    
        else:
            print('Samplesheet generation error: --pipeline input is not bacass or unicycler')
            break
    
    return samplesheet_string

def generate_samplesheet_long(nanopore_files, mapping_file, pipeline, genomesize):
    samplesheet_string = ''


    for line in mapping_file:
        if pipeline == 'bacass':
            samplesheet_string += str(id) + '\tNA\tNA\t' + nanopore_files + nanopore_ids + '.fastq\tNA\t' + str(genomesize) + '\n'
    
        elif pipeline == 'unicycler':
            samplesheet_string += str(id) + '\tNA\tNA\t' + nanopore_files + nanopore_ids + '.fastq\n'

        else:
            print('Samplesheet generation error: --pipeline input is not bacass or unicycler')
            break
    
    
    return samplesheet_string

def generate_samplesheet_hybrid(illumina_files, nanopore_files, mapping_file, pipeline, genomesize):
    samplesheet_string = ''
    illumina_read_files = [f for f in listdir(illumina_files) if isfile(join(illumina_files,f))]           # list of all illumina reads

    for line in mapping_file:
        split_line = line.split(',')        # split each line by separater (',' comma) --> expect mapping order sampleID, nanopore ID (barcodeXX)
        illumina_reads_id = [file for file in illumina_read_files if file.contain(split_line[0])] # get Illumina reads matching to the sample ID

        if pipeline == 'bacass':

            samplesheet_string += str(split_line[0]) + '\t' + illumina_reads_id[0] + '\t' + illumina_reads_id[1] + '\t' + nanopore_files + split_line[1] + '.fastq\tNA\t' + str(genomesize) + '\n'

        elif pipeline == 'unicycler':
            samplesheet_string += str(split_line[0]) + '\t' + illumina_reads_id[0] + '\t' + illumina_reads_id[1] + '\t' + nanopore_files + split_line[1] + '.fastq\n'

def wrapper_process(illumina_files, nanopore_files, reads, pipeline, genomesize):
    samplesheet_string = ''
    
    if reads == 'short':
        print('Generate file for short read assembly')
        samplesheet_string =generate_samplesheet_short(illumina_files, mapping_file, pipeline, genomesize)

    elif reads == 'long':
        print('Generate file for long read assembly')
        samplesheet_string = generate_samplesheet_long(nanopore_files, mapping_file, pipeline, genomesize)
    
    elif reads == 'hybrid':
        print('Generate file for hybrid assembly')
        samplesheet_string = generate_samplesheet_hybrid
    else:
        print('Wrong assembly type. Enter one of the following flags: short, long, hybrid')
        return samplesheet_string

    return samplesheet_string
    
if __name__ == '__main__':

    print(sys.argv)
    # Argument order: $illumina_files $nanopore_files $mapping_file_ch $reads $pipeline $samplesheet_header $genomesize
    illumina_files = sys.argv[0]
    nanopore_files = sys.argv[1]
    mapping_file = open(sys.argv[2], 'r')
    reads = sys.argv[3]
    pipeline = sys.argv[4]

    samplesheet_string = ''
    if pipeline == 'bacass':
        samplesheet_string = samplesheet_header.replaceAll(',', '\t') + '\n'   # header for nf-core/bacass pipeline
    
    
