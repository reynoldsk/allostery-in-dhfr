McCormick et. al 2020

This is the Readme for the first step in the Turbidostat Lit/Dark DHFR-LOV2 chimera mutagenesis experiments (Step_1_Usearch). 



The purpose of this was to combine the separate reads of the FASTQ files using the program Usearch11.0.667. Usearch is a sequence analysis tool that is faster than blast and most importantly updates the Qscores when combining an overlapping paired end read. 


Usearch is hosted here: https://www.drive5.com/usearch/ with the publication at: https://doi.org/10.1093/bioinformatics/btq461


The binary (usearch11.0.667_i86linux32) is included in this directory. 


It along with the unedited fastq files from GeneWiz were put on the UTSouthwestern BioHPC which uses the SLURM workload manager. 


runBashScript was submitted which kicked off the bash script: UCOMBINER.bsh 


Files from genewiz were unjoined illumina fastq.gz files separated by index. The forward and reverse reads were then combined using usearch v11.0.667 using the i86linux32 package. Due to 32 bit memory limitations, the large files had to be split. (split -l 40000000 #bash) This is tolerated due to the fact that the forward and reverse pair are in the same position in the file. After they are joined the files are then merged back into a single fastq file (all within UCOMBINER.bsh).  This directory does not contain the fastq files input or output due to size reasons. 


The scripting here is done for automation. The essential command to usearch is: 
$USEARCH -fastq_mergepairs ${name}R1aa -reverse ${name}R2aa -fastqout ${name}Raa_out  -fastq_maxdiffs 30  -fastq_nostagger  -fastq_minmergelen 140 -fastq_maxmergelen 250 -report ${name}.txt >> Logfile_COMBINER.txt
As the script and input files use absolute pathing (actual directory paths used have replaced with /path/ ). This setup will require some modification if the user wants to run it on their own machine. 