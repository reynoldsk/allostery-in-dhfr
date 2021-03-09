McCormick et. al. 2020

This is the Readme for the fourth step in the Turbidostat Lit/Dark DHFR-LOV2 chimera mutagenesis experiments (Step_4_calculating_allostery). 


The python (3.4) script Growth_Rate_and_allostery.ipynb is written for jupyter notebooks, contains markdown, and can be run in this directory on a personal computer producing all the data needed for step 5, calculating allostery.
 

The notebook takes the hamming adjusted residue counts produced in step3, hamming_filtering, and computes relative growth rates and allosteric effects. The data is then saved into ./output as both pickle databases (.pkl files for saving python variables containing the data) as well as human readable text files. 


This data is then analyzed and used for figure generation in Step_5_analysis_and_figures. 


