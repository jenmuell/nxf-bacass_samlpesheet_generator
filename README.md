# nxf-bacass_samlpesheet_generator

A simple [nextflow](https://www.nextflow.io/) pipeline to generate the sample sheets for the [unicycler](https://github.com/rrwick/Unicycler), [nf-core/bacass](https://nf-co.re/bacass/2.0.0/parameters) and the [plasmIDent](https://github.com/imgag/plasmIDent) pipeline.
It can handle all types of assemblies and was tested with Illumina and Nanopore Reads.
Current the pipeline run on the local maschine and require an installation of python (v3.6). 

## Default parameters
--pipeline          to specify one of the above mentioned pipelines <br>
--reads             to specify the assembly type by the corresponding reads. You can selcet ***short, long or hybrid*** <br>
--mapping_file      specfy the file names for the corresponding reads. Default for Nanopore reads is ***barcode***

### The mapping file
For short- and long-read assemblies it has to contain a single column with the respectiv name for the short and long read samples. <br>
For the hybrid assembly it contain two tab-separated columns, were the first include all short-read and the second all long-read sample names.

## Running commands
For more information about the pipeline parameter run:
```bash
    nextflow run jenmuell/nxf-bacass_samplesheet_generator --help
```
