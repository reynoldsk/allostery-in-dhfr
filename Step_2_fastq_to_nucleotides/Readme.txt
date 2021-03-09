McCormick et. al. 2020

This is the Readme for the second step in the Turbidostat Lit/Dark DHFR-LOV2 chimera mutagenesis experiments (Step_2_fastq_to_nucleotides). 


The goal of this section was to take the fastq files merged by usearch in step 1 and convert them to a count of mutations and nucleotides sorted by sublibrary in a human readable text document. As this script has to parse billions of lines of fastq data, it was run in parallel for each Index on the UT Southwestern BioHPC which uses the SLURM workload manager. Due to the large size of the fastq files, they are not located in this directory. 


The runPython script was submitted to the slurm manager to start the python (3.4) script DL121_fastq_analysis.py in parallel with one instance per input file (included in this dectory, e.g. SL1T0_input.txt). The input file identifies the absolute path to the merged fastq file as well as the vial and timepoint that the fastq represents (T0 samples were pooled then divided up into vials for lit/dark and is instead divided by sublibrary) 


DL121_fastq_analysis.py is the only python script in this code compilation that is not written for Jupyter notebook and thus has no markdown. The script parses the fastq files and input files to identify the sublibrary, timepoint, and vial (which is lit/dark) from which the read originates. It then filters the coding region of the read (the parts of DHFR for that sub ibrary that were mutated) and discards any read that has a base call in this region with a Q score less than 30 (base call accuracy  > 99.9%. For more information on Q scores see:
https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf


The script then counts the number of wild type reads as well as the number of nucleotide mutations. Multiple mutations (i.e. multiple nucleotide changes not in the same codon) are discarded as multimutants. This is done because the library construction produces only single amino acid mutations and thus multimutants are probably from some sort of sequencing or PCR error (bad base call).  


The script then writes this information and information about the number of reads that failed, passed, and the % that passed the qscore filter threshold of 30. 


Currently the script only produces .txt files of nucleotide counts, but can also produce amino acid counts through uncommenting the second to last line. 


While the script on the BioHPC writes the files to the directory the script is in. Those output files have been moved to the output directory in this directory where the following scripts can access them (step 3). 


All absolute paths have been shortened to /path/ to remove PII. 
