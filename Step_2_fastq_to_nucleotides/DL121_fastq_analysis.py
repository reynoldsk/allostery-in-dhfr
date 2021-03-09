#The goal of this script is to take the fastq files and convert them to a count of mutations and nucleotides sorted by sublibrary
#This script takes a very long time to run, it is set up to run on the UT Southwestern bioHPC
#as it has to loop over billions of lines of fastq data, it was run in parallel for each Index

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import stats
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from itertools import islice


###Establishing  parameters of the experiment:### 


#sequence of the wild type amplicon from 4N to 4N for each of the four sublibraries. 
#due to read length restrictions for high quality base pair calls on the illumina platform,
#the DHFR-LOV2 construct is broken up into 4 sublibraries which in total cover all 159 positions on DHFR that are mutated. 

#The reference wild type for each sublibrary below is what the unmutated DHFR region appears and is the defalt expectation
#due to multiplexing limitations of available indexs, the sublibraries are combined, they are id'ed in this script. 

#on NNNN: first 4 nucleotides are random from primers to generate high base diversity for read calibration. 
wt_ref_sl1 = 'NNNNATCACCATCATCACCACAGCCAGGATCCGATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCGATCGGTNNNN'

wt_ref_sl2 = 'NNNNGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCGATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTGAGCTCACAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCAATTGCGGCGTGTGGTGACGTACCAGAAATCNNNN'

wt_ref_sl3 = 'NNNNGTAACGTGGGTGAAGTCGGTGGATGAAGCAATTGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGCCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAGCTTTATCTGACGCATATCGACGCAGAAGTGGAACTGGCCACCACTCTAGAGCGCATCGAGNNNN'           

wt_ref_sl4 = 'NNNNAAGAAGACCGCCGAGAACATCGACGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAACATCCGTCGACAAGCTTGCGGCCGCATAATGCTTAAGNNNN'

#indices of amplicon where the coding region of interest begins & ends
#this is to ignore the NNNN and Homology region which is outside the sequence of interest (for primer binding)
sl1_coding_ix = [33,33+120]
sl2_coding_ix = [27,27+120]
sl3_coding_ix = [31,31+120] #last two amino acids are from LOV Insertion,
sl4_coding_ix = [28,28+117] 

#pulling out the wt sequence of the coding regions of interest
sl1_coding_region = wt_ref_sl1[sl1_coding_ix[0]:sl1_coding_ix[1]]
sl2_coding_region = wt_ref_sl2[sl2_coding_ix[0]:sl2_coding_ix[1]]
sl3_coding_region = wt_ref_sl3[sl3_coding_ix[0]:sl3_coding_ix[1]]
sl4_coding_region = wt_ref_sl4[sl4_coding_ix[0]:sl4_coding_ix[1]]

#SL1 amino acid positions
#these cover residue positions 1-40 in DHFR
sl1_aa_pos = np.arange(1,40+1) 
sl1_nuc_pos = np.arange(1,120+1) #changed for nucleotides (multiplied by 3)

#SL2 amino acid positions
#these cover positions 41-80 in DHFR 
sl2_aa_pos = np.arange(41,80+1)
sl2_nuc_pos = np.arange(121,240+1) #changed for nucleotides (multiplied by 3)

#SL3 amino acid positions
#these cover positions 81-120 in DHFR 
sl3_aa_pos = np.arange(81,120+1)
sl3_nuc_pos = np.arange(241,360+1) #changed for nucleotides (multiplied by 3)

#SL4 amino acid positions
#these cover positions 121-159 in DHFR
sl4_aa_pos = np.arange(121,159+1)
sl4_nuc_pos = np.arange(361,475+1) #changed for nucleotides (multiplied by 3)


###Setting up definitions###

#function to import names of the fastq files to open
def fileNames(line):
    spLine = line.split('\t')
    #print(spLine)
    fwdFileName = spLine[0]
    slIndexName = spLine[1]
    outName = spLine[2].strip('\n')
    return fwdFileName, slIndexName, outName

#function to check and fail a read based on the quality (Q) score. 
def qscore_filter(fwd):
    qscores = fwd #concatenating both reads
    low_quality = 'False' #default assumes Qscore is unacceptable
    for qscore in qscores: 
        if (ord(qscore) -33) <= qthreshold: #ord() converts ASCII encoding to numerical qscore
            low_quality = 'True'
    return low_quality

#function to take a full read and identify which sublibrary it belongs to. 
def id_sublibrary(read): 
    #Coding Regions of the protein. These are expected to have one mutation, which is two separate regions are checked.
    #These are 'fingerprints' which serve to identify which sublibrary the read belongs to.
    sl1_id_A = 'ATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAA' 
    sl1_id_B = 'ATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAA' 
    
    sl2_id_A = 'ATGGGCCGCCATACCTGGGAATCGATCGGTCGTCCGTTGCCAGGACGCAAAAATA' 
    sl2_id_B = 'CCTGAGCTCACAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGAT' 
    
    sl3_id_A = 'TTGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGCCGCGTTTATG'
    sl3_id_B = 'CTTGCCAAAAGCGCAAAAGCTTTATCTGACGCATATCGACGCAGAAG' 
    
    sl4_id_A = 'GCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAA' 
    sl4_id_B = 'CGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGC' 
        
    sublib_id = 'NEITHER' #default state of the id. 
    coding_ix = [0,0]

    #SL1
    stringer = read.find(sl1_id_A)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl1'
        coding_ix = [stringer-3,stringer+117]
        return sublib_id, coding_ix
        
    else:
        stringer = read.find(sl1_id_B)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl1'
        coding_ix = [stringer-57,stringer+len(sl1_id_B)+6]
        return sublib_id, coding_ix 
    else:
    #SL2
        stringer = read.find(sl2_id_A)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl2'
        coding_ix = [stringer-3,stringer+len(sl2_id_A)+62]
        return sublib_id, coding_ix 
    else:
        stringer = read.find(sl2_id_B)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl2'
        coding_ix = [stringer-62,stringer+len(sl2_id_B)+3]
        return sublib_id, coding_ix 
    else:
    #SL3
        stringer = read.find(sl3_id_A)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl3'
        coding_ix = [stringer-4,stringer+116]
        return sublib_id, coding_ix 
    else:
        stringer = read.find(sl3_id_B)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl3'
        coding_ix = [stringer-68,stringer+len(sl3_id_B)+5]
        return sublib_id, coding_ix 
    else:
    #SL4
        stringer = read.find(sl4_id_A)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl4'
        coding_ix = [stringer-1,stringer+116]
        return sublib_id, coding_ix 
    else:
        stringer = read.find(sl4_id_B)
    if not stringer == -1 and sublib_id == 'NEITHER':
        sublib_id = 'sl4'
        coding_ix = [stringer-62,stringer+len(sl4_id_B)+5]
        return sublib_id, coding_ix     

    return sublib_id, coding_ix 

#data metrics for each file looked at. 
def return_stats():
    print( 'read_total: ' + str(read_total)) 
    print( 'read_fail: ' + str(read_fail))
    print( 'read_pass: '+ str(read_pass))
    print( 'read_unassigned: '+ str(read_unassigned))
    print( 'qscore_used was: '+ str(qthreshold))
    print( str(percent_pass)+'% reads pass qscore filter threshold - '+str(qthreshold))

#this returns only the coding region of the protien
def trim_read(fwd):
    fwd_trim = fwd[coding_region[0]:coding_region[1]]         
    return fwd_trim

#turn the DNA into protein
def translate_seq(seq):
    seq = Seq(seq,generic_dna)
    seq_translate = str(seq.translate().strip())
    return seq_translate

#Identifies the mutation type both as codon and amino acid change
def identifyMutant(seq_fwd,reference_fwd,seq_len):
    mutants = []
    nuc_list = []
    aa_seq_fwd = translate_seq(seq_fwd)
    nuc_seq_fwd = seq_fwd
    aa_reference_fwd = translate_seq(reference_fwd)
    nuc_reference_fwd = reference_fwd
    def fwdMut():
        for ix,aa in enumerate(aa_seq_fwd):        
            ione = (np.multiply(ix,3))
            itwo = (np.multiply(ix,3) +1)
            ithree = (np.multiply(ix,3) +2)
            if (nuc_seq_fwd[ione]+nuc_seq_fwd[itwo]+nuc_seq_fwd[ithree]) != (nuc_reference_fwd[ione]+nuc_reference_fwd[itwo]+nuc_reference_fwd[ithree]):
                if sl_id == 'sl1':
                    mutants.append(aa_reference_fwd[ix]+str(sl1_aa_pos[ix])+aa)
                    nuc_list.append(nuc_reference_fwd[ione]+nuc_reference_fwd[itwo]+nuc_reference_fwd[ithree]                                        +str(sl1_aa_pos[ix])+(nuc_seq_fwd[ione]+nuc_seq_fwd[itwo]+nuc_seq_fwd[ithree]))
                elif sl_id == 'sl2':
                    mutants.append(aa_reference_fwd[ix]+str(sl2_aa_pos[ix])+aa)
                    nuc_list.append(nuc_reference_fwd[ione]+nuc_reference_fwd[itwo]+nuc_reference_fwd[ithree]                                        +str(sl2_aa_pos[ix])+(nuc_seq_fwd[ione]+nuc_seq_fwd[itwo]+nuc_seq_fwd[ithree]))
                elif sl_id == 'sl3':
                    mutants.append(aa_reference_fwd[ix]+str(sl3_aa_pos[ix])+aa)
                    nuc_list.append(nuc_reference_fwd[ione]+nuc_reference_fwd[itwo]+nuc_reference_fwd[ithree]                                        +str(sl3_aa_pos[ix])+(nuc_seq_fwd[ione]+nuc_seq_fwd[itwo]+nuc_seq_fwd[ithree]))
                elif sl_id == 'sl4':
                    mutants.append(aa_reference_fwd[ix]+str(sl4_aa_pos[ix])+aa)
                    nuc_list.append(nuc_reference_fwd[ione]+nuc_reference_fwd[itwo]+nuc_reference_fwd[ithree]                                        +str(sl4_aa_pos[ix])+(nuc_seq_fwd[ione]+nuc_seq_fwd[itwo]+nuc_seq_fwd[ithree]))
   
    #compare reference and sequence read
    if nuc_seq_fwd == nuc_reference_fwd:
        nuc_list.append('WT') 
        mutants.append('WT')
    else:
        fwdMut()    
    return mutants, nuc_list

#for each fastq file, make a fresh mutant_counts dictionary and counting residue mutants by parsing the mutant dictionary 
def record_mut_counts(mutant_dict):
    
    mutant_counts = {'SL1':{},'SL2':{},'SL3':{},'SL4':{}
                     
    def count_mutant(mutant,sl):
        if mutant in mutant_counts[sl].keys():
            mutant_counts[sl][mutant] +=1
        else:
            mutant_counts[sl][mutant] = 1

    for sl in mutant_dict.keys(): #iterate through each sub-library
        for mut_list in mutant_dict[sl]:
            if len(mut_list) == 1:
                count_mutant(mut_list[0],sl)
            elif len(mut_list) >=2: 
                count_mutant('fail_multimutant',sl) #only one mutation inserted in library per plasmid, so multiple is removed. 
                
    return mutant_counts
#for each fastq file, make a fresh mutant_counts dictionary and counting codon mutants by parsing the mutant dictionary 
def record_nuc_counts(nuc_list_dict):
                     
    nuc_counts = {'SL1':{},'SL2':{},'SL3':{},'SL4':{}

    def count_nuc(nuc,sl):
        if nuc in nuc_counts[sl].keys():
            nuc_counts[sl][nuc] +=1
        else:
            nuc_counts[sl][nuc] = 1

    for sl in nuc_dict.keys(): #iterate through each sub-library
        for nuc_list in nuc_dict[sl]:
            if len(nuc_list) == 1:
                count_nuc(nuc_list[0],sl)
            elif len(nuc_list) >=2: 
                count_nuc('fail_multimutant',sl) #only one mutation inserted in library per plasmid, so multiple is removed. 
                
    return nuc_counts

#write statistics & mutant residue counts into one .txt file                    
def writeOutputFile(outName):
    
    output_file = open(outName+'.txt','w')
    #write out statistics
    output_file.write('read_total:\t'+str(read_total)+'\t'+str(sl_id)+'\n')
    output_file.write('read_fail:\t'+str(read_fail)+'\n')
    output_file.write('read_pass:\t'+str(read_pass)+'\n')
    output_file.write('read_unassigned:\t'+str(read_unassigned)+'\n')
    output_file.write('qscore_used was:\t'+str(qthreshold)+'\n')
    output_file.write(str(percent_pass)+'\t% reads pass qscore filter threshold -\t'+str(qthreshold)+'\n')

    #write out mutant counts 
    for sl in mutant_counts.keys():
        for key in mutant_counts[sl].keys():
            output_file.write(sl+'\t'+key+'\t'+str(mutant_counts[sl][key])+'\n')

    output_file.close()
                  
#write statistics & mutant codon counts into one .txt file      
def writeOutputFile_nuc(outName):

    #writes read statistics & mutant counts into one convenient .txt file  
    
    output_file = open(outName+'nuc'+'.txt','w')
    #write out statistics
    output_file.write('read_total:\t'+str(read_total)+'\t'+str(sl_id)+'\n')
    output_file.write('read_fail:\t'+str(read_fail)+'\n')
    output_file.write('read_pass:\t'+str(read_pass)+'\n')
    output_file.write('read_unassigned:\t'+str(read_unassigned)+'\n')
    output_file.write('qscore_used was:\t'+str(qthreshold)+'\n')
    output_file.write(str(percent_pass)+'\t% reads pass qscore filter threshold -\t'+str(qthreshold)+'\n')

    #write out mutant counts 
    for sl in nuc_counts.keys():
        for key in nuc_counts[sl].keys():
            output_file.write(sl+'\t'+key+'\t'+str(nuc_counts[sl][key])+'\n')

    output_file.close()
                  
###Portion of code that loops through each line and executes the identification and categorization defined above###
                  
input_filenames = open(sys.argv[1]).readlines()
mutant_dict = {'SL1':[],'SL2':[],'SL3':[],'SL4':[]
nuc_dict = {'SL1':[],'SL2':[],'SL3':[],'SL4':[]

for l in input_filenames: #this unpacks each sample's reads, sublibary, and output files. 
    line = fileNames(l) #open fastq files for this sample

    slIndexName = line[1]

    #initialize output filename
    outName = line[2].strip('\n') #FILES MUST BE TAB DELINIATED.

    #thse flags are to identify what each line in the Fastq file is. (clusterID, break, coding region, Qscore). 
    seq_flag = 0 #if 0, you are not on a sequence line, if 1, you are on the sequence line
    seq_tracker = 0

    read_unassigned = 0
    read_total = 0
    read_pass = 0
    read_fail = 0
    qthreshold = 30 #Minimum of Phred Quality Score of 30 or must have more than 99.9% confidence in each base call in coding region.
            #for more on Qscore: https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf
    i = 0

    with open(line[0]) as fwdFile: #this bit is to read huge files 4 lines at a time. 
        #There is no line recognition beyond this. Illumina + flash + ucombine outputs use the same order, if different: be aware.
        while True:
            next_n_lines= list(islice(fwdFile, 4))
            if not next_n_lines:
                break
            read_total += 1
            fwd_seq = next_n_lines[1]
            qscore_line_fwd = next_n_lines[3]

            #want to get length of sequence
            seq_length = len(fwd_seq)
            #determine if it is in sublibrary1 or sublibrary2 
            sl_id, coding_region = id_sublibrary(fwd_seq)

            if sl_id == 'sl1':
                #trim to coding regions of interest
                fwd_seq_trim = trim_read(fwd_seq)
                qscore_line_fwd_trim = trim_read(qscore_line_fwd)
                fail = qscore_filter(qscore_line_fwd_trim) #call function for qscore filter 
                if fail == 'True':
                    read_fail+=1
                    continue

                else: 
                    #function for determining mutant
                    reference_fwd = sl1_coding_region
                    mutants, nuc_list = identifyMutant(fwd_seq_trim,reference_fwd,seq_length)
                    mutant_dict['SL1'].append(mutants)
                    nuc_dict['SL1'].append(nuc_list)

            elif sl_id == 'sl2':
                #trim to coding regions of interest
                fwd_seq_trim = trim_read(fwd_seq)
                qscore_line_fwd_trim = trim_read(qscore_line_fwd)

                fail = qscore_filter(qscore_line_fwd_trim) #call function for qscore filter 
                if fail == 'True':
                    read_fail+=1
                    continue

                else:#function for determining mutant
                    reference_fwd = sl2_coding_region
                    mutants, nuc_list = identifyMutant(fwd_seq_trim,reference_fwd,seq_length)
                    mutant_dict['SL2'].append(mutants)
                    nuc_dict['SL2'].append(nuc_list)

            elif sl_id == 'sl3':
                #trim to coding regions of interest
                fwd_seq_trim = trim_read(fwd_seq)
                qscore_line_fwd_trim = trim_read(qscore_line_fwd)
                fail = qscore_filter(qscore_line_fwd_trim) #call function for qscore filter 
                if fail == 'True':
                    read_fail+=1
                    continue

                else:#function for determining mutant
                    reference_fwd = sl3_coding_region
                    mutants, nuc_list = identifyMutant(fwd_seq_trim,reference_fwd,seq_length)
                    mutant_dict['SL3'].append(mutants)
                    nuc_dict['SL3'].append(nuc_list)

            elif sl_id == 'sl4':
                #trim to coding regions of interest
                fwd_seq_trim = trim_read(fwd_seq)
                qscore_line_fwd_trim = trim_read(qscore_line_fwd)
                fail = qscore_filter(qscore_line_fwd_trim) #call function for qscore filter 
                if fail == 'True':
                    read_fail+=1
                    continue

                else:  #function for determining mutant
                    reference_fwd = sl4_coding_region
                    mutants, nuc_list = identifyMutant(fwd_seq_trim,reference_fwd,seq_length)
                    mutant_dict['SL4'].append(mutants)
                    nuc_dict['SL4'].append(nuc_list)

            else: #if it cant figure out what the sublibrary is, it is counted here. e.g. a mutation or bad call in both fingerprints
                read_fail+=1
                read_unassigned+=1
                continue

            if seq_flag == 0:
                if'@M02894' in line: #sequence is right after this id line
                    seq_flag = 1
                else: 
                    continue

    mutant_counts = record_mut_counts(mutant_dict)
    nuc_counts = record_nuc_counts(nuc_dict)

    read_pass = read_total-read_fail
    percent_pass = read_pass*100/read_total
    return_stats()

    #writing statistics and mutant counts to a txt file 
    #writeOutputFile(outName) #this has been commented out. In step 3 we hamming filter than convert to residues. 
            #if you want to get residues without any filtering or to compare the list. Uncommenting this will produce them. 
    writeOutputFile_nuc(outName)

