This file contains details about BE623_Project1 - Sequence analysis 
Group members - Riya,Pragya, Swathi
The genes analysed are:
SELENOK(Seelnoprotein)- Riya 
MT-ND5(Mitochondrial gene)- Pragya
TP53(Nuclear gene)- Swathi

Project Pipeline
This pipeline can take in mRNA sequence of the gene as input and produce its corresponding protein sequence as output. Further the translated protein sequence is compared against protein sequence retrievd from NCBI database.

Steps:
1. Takes nucleotide fasta file for a gene as input 
2. Extracts the CDS region based on the coordinates given as start and end position
- Checks if CDS region is extracted correctly by checking its divisibilty by 3(1 codon = 3 nucleotide) 
3. Translates the CDS region using custom made genetic code dictionary
- Standard genetic code for nuclear gene.
- For Selenoprotein, TGA assigned as U(Selenocysteine) instead of stop codon.
- For Mitochondrial gene, AGA/AGG acts as stop codon istead of coding for Arginine and TGA codes for Tryptophan instead of stop codon and ATA codes for methionine instead of Isoleucine.
4. Compares custom translated sequence with retrieved sequence from NCBI protein database using Zip function, if difference present at any position, returns number of differences.
5. Presents the results as a summary table which gets downloaded as .tsv file displaying the following information for each gene:
- Gene Name 
- CDS length(nt)
- Translated protein length 
- Deposited protein length 
- Total differences
- First mismatch 
- Status (match/mismatch)

Tools Used:
NCBI nucleotide and Protein databases
Python
Git and GitHub


Commands used:
Functions in Python:
1. read_fasta()- to read the fasta file
2. extract_cds()- to extract the cds region based on start and end 
3. check_cds_length- to check the length of the extracted CDS region
4. translate_cds- to translate the extracted cds according to our genetic code dictionary 
5. compare_proteins()- compares the translated sequence to the retrieved NCBI sequence
6. process_gene()- main function that runs the analysis for one gene
