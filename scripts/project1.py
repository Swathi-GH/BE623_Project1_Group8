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


def extract_cds(sequence, start, end):
	return sequence[start - 1:end]



