This is a list of instructions for running the pipeline in its current form as seperate scripts.
A script is included for generating .bam files from the raw reads and assembly as blobtoolkit required
as specific style of indexing where a .bam.csi is generated instead of a .bam.bai. 

1.The script "tiara.sh" is used to run tiara on the assembly file. it requires tiara to be installed, and takes
an input directory and output directory as inputs. It is designed to loop through a given directory, as tiara is very
quick. FCS-GX is much much slower so processes 1 file at a time.

2. The script "FCS-GX.sh" requires an input fasta file and  NCBI taxonid, as well as a installation of the FCS-GX database
"/gxdb/all.gxi". If not taxon ID is available your target species, you can provide the taxon id for the genus of your target species.

FCS-GX outputs 2 files. one is a .txt file containing the contigs FCS-GX has classified as "contamination". the second is
a .rpt file containing the classifications of all the contigs. This is what the pipeline uses, and requires some processing 
to be readable by R, which is carried out by FCS_reformat.sh. So once you have sucessfully run the FCS-GX.sh script, you then run FCS_reformat.sh.

3. There is then an R script which processed the results of FCS-GX and Tiara and compares their outputs. This generates several
reports as well as a final TSV which is what blobtools uses to tag contigs with taxonomy. To run the R script you need to open it and 
add your file path for tiara and fcs results and an out directory. then you need to run the "comparison_R.sh" script, as this is a wrapper 
script which runs R for you via slurm.

4. At this point you should have a .bam file, a fasta.fa file, a taxonomy.tsv file. you will also need to create a .yml file. simply replace
the species name and taxon id in the provided example .yml file. 


5. You will then input the file paths for these along with your busco.tsv (if you have it, only the .bam &.fasta files are essential) into the blobtoolkit (btk) script.
This script builds the blobdirectoryand adds each data set to it. you will need to install blobtoolkit into a mamba env (issues when using btk and conda) and then tell the script where
you have installed it. You then should run the blobtoolkit command "blobtools validate ." from inside your new blobdirectory. This will check if all of
your .json files (how blobtools stores your input info) were created correctly. if correct, you will then have to transfer the blobdirectory to your local
machine and run btk on a docker container, as btk has some issues displaying on the cluster. The docker command will be provided here also.


docker command to be run on local machine "docker run -d --rm --name EGP017_25_Com1_blobdir     -v ~/FSP/contamination_detection/btk:/blobtoolkit/datasets     -p 8000:8000 -p 8080:8080     -e VIEWER=true genomehubs/blobtoolkit:latest"


once you have successfully opened your blobdirectory on localhost, you use the user interface to manually trim which contigs you wish to include. once that is done you can export the remaining "clean"
contigs as a fasta file.
