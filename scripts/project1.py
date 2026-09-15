from Bio.Seq import Seq

#input files and CDS information as dictionary

genes = {
    "TP53": {
        "category": "nuclear",
        "mrna_file": "data/TP53_mRNA.fasta",
        "protein_file": "data/TP53_protein.fasta",
        "cds_start": 143,
        "cds_end": 1324
    },

    "SELENOK": {
        "category": "selenoprotein",
        "mrna_file": "data/selenok.fasta",
        "protein_file": "data/selenok_prtn.fasta",
        "cds_start": 72,
        "cds_end": 356
    },

    "MT-ND5": {
        "category": "mitochondrial",
        "mrna_file": "data/nd5.fasta",
        "protein_file": "data/nd5_protein.fasta",
	"cds_start": 1,
	"cds_end": 1812
    },
}

#creating function to read the fasta

def read_fasta(filename):
    with open(filename, "r") as fh:
        lines = fh.readlines()

    sequence = ""

    for line in lines:
        if not line.startswith(">"):
            sequence = sequence + line.strip()

    return sequence

#creating function to extract cds

def extract_cds(sequence, start, end):
	return sequence[start - 1:end]


#checking cds length to see if it is divisible by 3

def check_cds_length(cds_sequence):
	length = len(cds_sequence)

	if length % 3 == 0:
		return length, True
	else:
		return length, False

#creating genetic code


genetic_code_nuclear = {
    "TTT" : "F",
    "TTC" : "F",
    "TTA" : "L",
    "TTG" : "L",
    "CTT" : "L",
    "CTC" : "L",
    "CTA" : "L",
    "CTG" : "L",
    "ATT" : "I",
    "ATC" : "I",
    "ATA" : "I",
    "ATG" : "M",
    "GTT" : "V",
    "GTC" : "V",
    "GTA" : "V",
    "GTG" : "V",
    "TCT" : "S",
    "TCC" : "S",
    "TCA" : "S",
    "TCG" : "S",
    "CCT" : "P",
    "CCC" : "P",
    "CCA" : "P",
    "CCG" : "P",
    "ACT" : "T",
    "ACC" : "T",
    "ACA" : "T",
    "ACG" : "T",
    "GCT" : "A",
    "GCC" : "A",
    "GCA" : "A",
    "GCG" : "A",
    "TAT" : "Y",
    "TAC" : "Y",
    "TAA" : "*",
    "TAG" : "*",
    "CAT" : "H",
    "CAC" : "H",
    "CAA" : "Q",
    "CAG" : "Q",
    "AAT" : "N",
    "AAC" : "N",
    "AAA" : "K",
    "AAG" : "K",
    "GAT" : "D",
    "GAC" : "D",
    "GAA" : "E",
    "GAG" : "E",
    "TGT" : "C",
    "TGC" : "C",
    "TGA" : "*",
    "TGG" : "W",
    "CGT" : "R",
    "CGC" : "R",
    "CGA" : "R",
    "CGG" : "R",
    "AGT" : "S",
    "AGC" : "S",
    "AGA" : "R",
    "AGG" : "R",
    "GGT" : "G",
    "GGC" : "G",
    "GGA" : "G",
    "GGG" : "G",
}










genetic_code_selenocys = {
    #Phenylalanine (F)
    "TTT": "F", "TTC": "F",
    #Leucine (L)
    "TTA": "L", "TTG": "L", "CTT": "L",  
    "CTC": "L", "CTA": "L", "CTG": "L",
    #Isoleucine (I)
    "ATT": "I", "ATC": "I", "ATA": "I",
    #Methionine (M) which is a start codon
    "ATG": "M",
    #Valine (V)
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    #Serine (S)
    "TCT": "S", "TCC": "S", "TCA": "S",
    "TCG": "S", "AGT": "S", "AGC": "S",
    #Proline (P)
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    #Threonine (T)
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    #Alanine (A)
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    #Tyrosine (Y)
    "TAT": "Y", "TAC": "Y",
    #Histidine (H)
    "CAT": "H", "CAC": "H",
    #Glutamine (Q)
    "CAA": "Q", "CAG": "Q",
    #Asparagine (N)
    "AAT": "N", "AAC": "N",
    #Lysine (K)
    "AAA": "K", "AAG": "K",
    #Aspartic acid (D)
    "GAT": "D", "GAC": "D",
    #Glutamic acid (E)
    "GAA": "E", "GAG": "E",
    #Cysteine (C)
    "TGT": "C", "TGC": "C",
    #Tryptophan (W)
    "TGG": "W",
    #Arginine (R)
    "CGT": "R", "CGC": "R", "CGA": "R", 
    "CGG": "R", "AGA": "R", "AGG": "R",
    #Glycine (G)
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    #Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*",     #for selenocysteine TGA codes for U
}

















def translate_cds(cds_sequence, genetic_code):
	protein = ""

	for i in range(0, len(cds_sequence), 3):
		codon = cds_sequence[i:i+3]
		amino_acid = genetic_code[codon]
		protein = protein + amino acid

	return protein


def compare_proteins(protein1, protein2):
    differences = 0
    first_mismatch = 0
    position = 1

    for a, b in zip(protein1, protein2):
        if a != b:
            differences = differences + 1

            if first_mismatch == 0:
                first_mismatch = position

        position = position + 1

    return differences, first_mismatch

def process_gene(gene_name, gene_info, genetic_code):
    mrna = read_fasta(gene_info["mrna_file"])

    cds_sequence = extract_cds(
        mrna,
        gene_info["cds_start"],
        gene_info["cds_end"]
    )

    translated_protein = translate_cds(cds_sequence, genetic_code)

    deposited_protein = read_fasta(gene_info["protein_file"])

    differences, first_mismatch = compare_proteins(
        translated_protein,
        deposited_protein
    )

    return cds_sequence, translated_protein, deposited_protein, differences, first_mismatch


for gene_name in genes:
    gene_info = genes[gene_name]
    genetic_code = genetic_codes[gene_name]

    results = process_gene(
        gene_name,
        gene_info,
        genetic_code
    )
