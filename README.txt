This is a list of instructions for running the pipeline in its current form as seperate scripts
This also includes a run through of the steps of this process as is

1. First script will take all of the genomes in a given directory, create copies with all contigs
that are <250bp removed. This is as EBI does not accept contigs below this size. The script will then
generate BAM files from these.

2.the next script will generate new copies of these files that are seperated in >1kbp and <1kbp. Tiara
will then be run to generate a Tiara report.

3. A .yml file will be needed to run FCS, so this will have to be manually created at present.

4. FCS can then be run on the >1kbp files to generate a FCS report.

5. The comparison script will then be run on the tiara and fcs-gx report. the result of this will be 
a taxonomic assignments tsv file, which will be fed into blobtools.

6.the >250bp version of the fasta files will then be fed into blobtools, along with the .bam files and the
taxonomic assignment tsv using the provided script.

7. contigs flagged as contaminants by this process will then be removed from the >250bp fasta files, which
will need another assessment step.



 
