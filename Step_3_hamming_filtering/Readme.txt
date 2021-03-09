McCormick et. al. 2020

This is the Readme for the third step in the Turbidostat Lit/Dark DHFR-LOV2 chimera mutagenesis experiments (Step_3_hamming_filtering). 


The Hamming_analysis.ipynb (Python 3.4) script is written for jupyter notebooks and contains markdown and can be run in this directory on a personal computer producing all the data needed for step 4, calculating allostery. 


The notebook takes the nucleotide counts produced in step 2, FastQ to Nucleotides and converts them to the amino acid counts while accounting for and subtracting expected error in illumina sequencing, termed Hamming Filtering. 
Previously in FastQ to Nucleotides (analysis step 2) we threw out any read where the base calls in the coding region (or coding region of the sublibrary of interest, specifically) had a Q score of less than 30. And recorded the number of each nucleotide mutation as well as wild type counts in .txt files in ../Step_2_fastq_to_nucleotides/output. The python notebook takes these output files as names in ./masterfile.txt and performs hamming analysis by each of the 37 out files (T0_sl1 T0_sl2 T0_sl3 T0_sl4 have been merged into a single file, T0). 


Hamming Analysis, created by James McCormick (in this context) is the process for accounting for the error rate in Illumina sequencing. Because each read is counted as a single mutant, a base call with a quality score Q40 has a 99.99% base call accuracy. But when there are a large number of reads, numbering into the billions for this set of experiments, the contribution of miscalled mutations that are actually wild type (or a mutant that is misread as wild type) becomes an important contribution.


This source of noise is also not uniform. If we sequence only an identical wild type sequence with a very large number of reads, apparent mutations that require two adjacent mutations are much less likely than ones that are a single base miscall away. Such as ATG being called falsely as ATC (met to ile) vs ATG to AAC (met to asn).


From the Hamming distance Wikipedia article: "In information theory, the Hamming distance between two strings of equal length is the number of positions at which the corresponding symbols are different. In other words, it measures the minimum number of substitutions required to change one string into the other, or the minimum number of errors that could have transformed one string into the other. In a more general context, the Hamming distance is one of several string metrics for measuring the edit distance between two sequences. It is named after the American mathematician Richard Hamming."


The probability of two errant base calls happening (a hamming distance of two) is 0.0001 x 0.0001, three errant base calls is 0.0001 x 0.0001 x 0.0001 or P^H where P is the probability of a errant base call and H is the hamming distance.


Once we have the probability of an errant mutant stemming from another sequence, we can calculate the number of errant calls we would expect from a given number of other counts, e.g. wild type. If there are 1,000,000 reads of the wild type sequence and mutation A is one hamming distance from wild type, we can estimate that there are (1,000,000*0.0001^1) = 100 mutant A's that are errantly called wild type. We can repeat this for every observed mutation type at that codon.

This process makes a few simplifying assumptions:

1.	That every other observed mutant/wild type is correct.

2.	There are no mutants that are misread in one location and then misread again to produce this mutation. e.g. ATTCGG --> ATCCGA.

3.	The probability of the sequencer missing a base call is equal to the resultant mutation. i.e. that an erroneous ATC -> ATG is equally as likely as an erroneous ATC -> ATT

4.	That every base call has a Q score equal to the average as reported by genewiz This is not true, and is a source of bias, as the average Q score is better at the start of a read, which is why USEARCH (step 1) and paired end coverage is employed. Scores below Q30 were thrown out in the previous script and USEARCH increases the Q score from that of illumina base calling.

These assumptions result in a "hard filtering" where more reads are thrown away than strictly need to be.

After the counts are adjusted they are converted to the amino acid count (counts of different codons encoding the same residue are combined) saved in text files corresponding to the timepont and vial and written to ./output where the growth rate and allosteric effect can be analyzed in step 4. 