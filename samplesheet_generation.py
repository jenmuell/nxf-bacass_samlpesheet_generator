#!/usr/bin/python

from os import listdir
from os.path import isfile, join
import sys

def generate_samplesheet_short(illumina_files, mapping_file, pipeline, genomesize):
    samplesheet_string = ''
    illumina_read_files = [f for f in listdir(illumina_files) if isfile(join(illumina_files,f))]           # list of all illumina reads
    for line in mapping_file:  # Read mapping file line by line
        illumina_reads_id = [file for file in illumina_read_files if line[:-1] in file]  # R1 and R1 for sample 
        if pipeline == 'bacass':
            samplesheet_string += str(line[:-1]) + '\t' + illumina_files + '/' + illumina_reads_id[0] + '\t' + illumina_files + '/' + illumina_reads_id[1] + '\tNA\tNA\t' + str(genomesize) + 'm\n'
    
        elif pipeline == 'unicycler':
            samplesheet_string += str(line[:-1]) + ',' + illumina_files + '/' + illumina_reads_id[0] + ',' + illumina_files + '/' + illumina_reads_id[1] + '\n'
    
        else:
            print('Samplesheet generation error: --pipeline input is not bacass or unicycler')
            break
    
    return samplesheet_string

def generate_samplesheet_long(nanopore_files, mapping_file, pipeline, genomesize):
    samplesheet_string = ''


    for line in mapping_file:
        if pipeline == 'bacass':
            samplesheet_string += str(line[:-1]) + '\tNA\tNA\t' + nanopore_files + line[:-1] + '.fastq\tNA\t' + str(genomesize) + 'm\n'
    
        elif pipeline == 'unicycler':
            samplesheet_string += str(line[:-1]) + ',' + nanopore_files + line[:-1] + '.fastq\n'

        else:
            print('Samplesheet generation error: --pipeline input is not bacass or unicycler')
            break
    
    
    return samplesheet_string

def generate_samplesheet_hybrid(illumina_files, nanopore_files, mapping_file, pipeline, genomesize):
    samplesheet_string = ''
    illumina_read_files = [join(illumina_files, f) for f in listdir(illumina_files) if isfile(join(illumina_files, f))]           # list of all illumina reads

    for line in mapping_file:
        split_line = line.replace('\n','').split('\t')        # split each line by separater (',' comma) --> expect mapping order sampleID, nanopore ID (barcodeXX)
        illumina_reads_id = [file for file in illumina_read_files if split_line[0] in file] # get Illumina reads matching to the sample ID
        if pipeline == 'bacass':
            samplesheet_string += str(split_line[0]) + '\t'  + illumina_reads_id[0] + '\t' +  illumina_reads_id[1] + '\t' + nanopore_files + split_line[1] + '.fastq\tNA\t' + str(genomesize) + 'm\n'

        elif pipeline == 'unicycler':
            samplesheet_string += str(split_line[0]) + ',' +  illumina_reads_id[0] + ','  + illumina_reads_id[1] + ',' + nanopore_files + split_line[1] + '.fastq\n'
            
    return samplesheet_string

def generate_samplesheet_plasmident(mapping_file, nanopore_files, assembly_files):
    samplesheet_string = ''

    for line in mapping_file:
        split_line = line.replace('\n', '').split('\t')
        if assembly_files == 'null' :
            samplesheet_string += str(split_line[0]) + '\t' + str(split_line[0]) + '/assembly.fasta\t' + nanopore_files + split_line[1] + '.fastq\n'
        else:
            samplesheet_string += str(split_line[0]) + '\t' + assembly_files + str(split_line[0]) + '/assembly.fasta\t' + nanopore_files + split_line[1] + '.fastq\n'

    return samplesheet_string

def wrapper_process(illumina_files, nanopore_files, assembly_files, mapping_file, reads, pipeline, genomesize):
    samplesheet_string = ''
    
    if pipeline == 'plasmident':
        samplesheet_string = generate_samplesheet_plasmident(mapping_file, nanopore_files, assembly_files)
    elif reads == 'short':
        samplesheet_string = generate_samplesheet_short(illumina_files, mapping_file, pipeline, genomesize)

    elif reads == 'long':
        samplesheet_string = generate_samplesheet_long(nanopore_files, mapping_file, pipeline, genomesize)
    
    elif reads == 'hybrid':
        samplesheet_string = generate_samplesheet_hybrid(illumina_files, nanopore_files, mapping_file, pipeline, genomesize)
    else:
        print('Wrong assembly type. Enter one of the following flags: short, long, hybrid')
        return samplesheet_string

    return samplesheet_string
    
if __name__ == '__main__':

    #print(sys.argv)
    # Argument order: $illumina_files $nanopore_files $mapping_file_ch $reads $pipeline $samplesheet_header $genomesize
    # argv[0] = path to this script
    illumina_files = sys.argv[1]
    nanopore_files = sys.argv[2]
    assembly_files = sys.argv[3]
    mapping_file = open(sys.argv[4], 'r')
    reads = sys.argv[5]
    pipeline = sys.argv[6]
    samplesheet_header = sys.argv[7]
    genomesize = sys.argv[8]

    samplesheet_string = ''
    if pipeline == 'bacass':
        samplesheet_string = samplesheet_header.replace(',', '\t') + '\n'   # header for nf-core/bacass pipeline
    elif pipeline == 'viralrecon':
        samplesheet_string = samplesheet_header + '\n' 
        pipeline = 'unicycler'
    samplesheet_string += wrapper_process(illumina_files, nanopore_files, assembly_files, mapping_file, reads, pipeline, genomesize)
    print(samplesheet_string)
