from Bio.Seq import Seq

#input files and CDS information as dictionary

genes = {
    "TP53": {
        "category": "nuclear",
        "mrna_file": "../data/TP53_mRNA.fasta",
        "protein_file": "../data/TP53_protein.fasta",
        "cds_start": 143,
        "cds_end": 1324
    },

    "SELENOK": {
        "category": "selenoprotein",
        "mrna_file": "../data/selenok.fasta",
        "protein_file": "../data/selenok_prtn.fasta",
        "cds_start": 72,
        "cds_end": 356
    },

    "MT-ND5": {
        "category": "mitochondrial",
        "mrna_file": "../data/nd5.fasta",
        "protein_file": "../data/nd5_protein.fasta",
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

#creating genetic code for nuclear gene


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








#creating genetic code for selenoprotein

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











#creating genetic code  mitochondria

genetic_code_mito = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "M",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCG": "A",
    "GCA": "A",
    "GCC": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "W",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "*",
    "AGG": "*",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def translate_cds(cds_sequence, genetic_code):
	protein = ""

	for i in range(0, len(cds_sequence), 3):
		codon = cds_sequence[i:i+3]
		amino_acid = genetic_code[codon]
		protein = protein + amino_acid

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

    return (
        cds_sequence,
        translated_protein,
        deposited_protein,
        differences,
        first_mismatch
    )


#genetic code for each gene
genetic_codes = {
    "SELENOK": genetic_code_nuclear,
    "MT-ND5": genetic_code_mito,
    "TP53": genetic_code_nuclear
}


# Process each gene and store results
summary = []

for gene_name in genes:
    gene_info = genes[gene_name]
    genetic_code = genetic_codes[gene_name]

    (
        cds_sequence,
        translated_protein,
        deposited_protein,
        differences,
        first_mismatch
    ) = process_gene(gene_name, gene_info, genetic_code)

    status = "Match" if differences == 0 else "Mismatch"

    summary.append([
        gene_name,
        len(cds_sequence),
        len(translated_protein),
        len(deposited_protein),
        differences,
        first_mismatch,
        status
    ])

    print("\n" + "=" * 50)
    print("Gene:", gene_name)
    print("CDS length:", len(cds_sequence), "nt")
    print("Translated protein length:", len(translated_protein))
    print("Deposited protein length:", len(deposited_protein))
    print("Total differences:", differences)
    print("First mismatch:", first_mismatch)
    print("Status:", status)


# Final summary table
print("\n" + "=" * 90)
print("SUMMARY TABLE")
print("=" * 90)

print(
    f'{"Gene":<12}'
    f'{"CDS(nt)":<10}'
    f'{"Translated":<12}'
    f'{"Deposited":<11}'
    f'{"Diff":<8}'
    f'{"First":<8}'
    f'{"Status"}'
)

print("-" * 90)

for row in summary:
    print(
        f"{row[0]:<12}"
        f"{row[1]:<10}"
        f"{row[2]:<12}"
        f"{row[3]:<11}"
        f"{row[4]:<8}"
        f"{str(row[5]):<8}"
        f"{row[6]}"
    )

# Save summary table
with open("../results/summary_table.tsv", "w") as fh:
    fh.write("Gene\tCDS(nt)\tTranslated\tDeposited\tDiff\tFirst\tStatus\n")

    for row in summary:
        fh.write(
            f"{row[0]}\t"
            f"{row[1]}\t"
            f"{row[2]}\t"
            f"{row[3]}\t"
            f"{row[4]}\t"
            f"{row[5]}\t"
            f"{row[6]}\n"
        )
