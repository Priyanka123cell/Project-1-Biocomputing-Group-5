print("PIPELINE ALL GENES")
print("DRD2 PIPELINE")

# 2(a) Extraction of cds from mRNA
with open("data_Priyanka_DRD2/refseq_mrna.fasta","r") as refseq:
    mrna =""
    for seq in refseq:
        if(seq[0]!=">"):
            mrna = mrna+seq.strip()

cds = mrna[354:1686]
print("Extracted coding sequence: DRD2")
print(cds)
cds_len = len(cds)
print("Length of DRD2 extracted coding sequence:",cds_len)
print("Number of codons in extracted coding sequence for DRD2:",cds_len / 3)

# 2(b) Translation of DRD2 CDS
codon_table = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

amino_seq =""
for r in range(0,len(cds),3):
    codon_seq = cds[r:r+3]
    if codon_seq in codon_table:
        amino_seq = amino_seq+codon_table[codon_seq]

print("Translated sequence manually: DRD2")
print(amino_seq)

amino_30 =""
for r in range(0,90,3):
    codon_30 = cds[r:r+3]
    if codon_30 in codon_table:
        amino_30 = amino_30+codon_table[codon_30]

len_amino_seq  = len(amino_seq)
print("Total length of translated sequence (DRD2): ",len_amino_seq)
print("First 30 residues (DRD2):", amino_30)

# 2(c) Deposited and manual translation comparison for DRD2
with open("data_Priyanka_DRD2/protein_sequence.fasta","r") as protseq:
    depo_seq =""
    for seq1 in protseq:
        if(seq1[0]!=">"):
            depo_seq = depo_seq+seq1.strip()

print("Deposited protein sequence: DRD2")
print(depo_seq)
depo_len = len(depo_seq)
print("Length of deposited sequence(DRD2) :",depo_len)
print("My translated sequence length DRD2:",len_amino_seq)

#IDENTICAL OR NOT AND MISMATCH POSITION
diff = 0
for a,d in zip(amino_seq,depo_seq):
    if a!=d:
        diff+=0

print("Difference between my translated sequence and deposited sequence:", diff)

drd2_identical_OR_NOT = "Y"
drd2_mismatch_pos = "None"
if len(amino_seq)==len(depo_seq):
    print(drd2_identical_OR_NOT)
    print(drd2_mismatch_pos)

else:
    drd2_identical_OR_NOT = "N"
    drd2_mismatch_pos = min(len(amino_seq), len(depo_seq)) + 1
    for i in range(min(len(amino_seq), len(depo_seq))):
        if (amino_seq[i] != depo_seq[i]):
            drd2_identical_OR_NOT = "N"
            drd2_mismatch_pos = i + 1
            break 

print("Identical (Y/N):",drd2_identical_OR_NOT)
print("Mismatch Position (due to stop codon at last position):", drd2_mismatch_pos)


#MTCYB
print("MTCYB_pipeline")

# 2(a) Extracting the CDS region for MTCYB from the mitochondrial reference genome
with open("Sparshika_MTCYB_data/Mitochondrial_genome.fasta", "r") as fh:
    genome_data = ""
    for sequence in fh:
        if(sequence[0]!=">"):
            genome_data = genome_data + sequence.strip()
New_sequence = genome_data[14746:15887]
New_sequence = New_sequence + "AA"

with open("Sparshika_MTCYB_data/MTCYB_CDS_corrected.fasta", "w") as fh:
    fh.write(New_sequence)

MTCYB_CDS_length = len(New_sequence)
print(MTCYB_CDS_length)
print(MTCYB_CDS_length / 3)

# 2(b) Mitochondrial genetic codon table
Table = {
    "TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
    "TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
    "TTA":"L", "TCA":"S", "TAA":"",  "TGA":"W",
    "TTG":"L", "TCG":"S", "TAG":"",  "TGG":"W",
    "CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
    "CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
    "CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
    "CTG":"L", "CCG":"P", "CAG":"Q", "CGC":"R",
    "ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
    "ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
    "ATA":"M", "ACA":"T", "AAA":"K", "AGA":"",
    "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"",
    "GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
    "GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
    "GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
    "GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
}

with open("Sparshika_MTCYB_data/MTCYB_CDS_corrected.fasta", "r") as fh:
    seq = fh.read().replace("\n", "")

Translated_protein_seq = ""
for i in range(0, len(seq), 3):
    codon = seq[i:i+3]
    AA_residue = Table[codon]
    Translated_protein_seq = Translated_protein_seq + AA_residue

print(Translated_protein_seq)
print("First 30 residues:", Translated_protein_seq[:30])
MTCYB_translated_length = len(Translated_protein_seq)
print("Total length of the protein sequence:", MTCYB_translated_length)

# 2(c) Comparing my translated sequence with the deposited protein
with open("Sparshika_MTCYB_data/MTCYB_protein.fasta", "r") as fh:
    deposited_protein = ""
    for seq1 in fh:
        if(seq1[0]!=">"):
            deposited_protein = deposited_protein + seq1.strip()

seq_1 = Translated_protein_seq
seq_2 = deposited_protein
print(MTCYB_translated_length)
MTCYB_deposited_length = len(deposited_protein)
print(MTCYB_deposited_length)
MTCYB_is_identical_to_deposited_protein = "Y"
MTCYB_mismatch_position = "N/A"

differences = 0
for a, b in zip(seq_1, seq_2):
    if a != b:
        differences = differences + 1

if MTCYB_translated_length != MTCYB_deposited_length:
	MTCYB_is_identical_to_deposited_protein  = "N"
print(differences)
print(MTCYB_mismatch_position)

#SELENOH
# CDS
print("PIPELINE FOR SELENOH")
start_coordinate=106;
end_coordinate=474;
length= ((end_coordinate-start_coordinate)+1)
print("Length of the coding sequence",length)
# CODONS
Codons=(length/3)
print("Number of codons",Codons)
# mRNA
mRNA= "AGAGCTTCCGGGCTGCGCTCTTCGTTGCCCAGTTTCCGCTCAGTGGTCGCGTCTCCGCCCCCCACCCACCAGTCCCGCTGCATTCTCGGCCGGGCTCTAGGCGCCATGGCTCCCCGCGGGAGGAAGCGTAAGGCTGAGGCCGCGGTGGTCGCCGTAGCCGAGAAGCGAGAGAAGCTGGCGAACGGCGGGGAGGGAATGGAGGAGGCGACCGTTGTTATCGAGCATTGCACTAGCTGACGCGTCTATGGGCGCAACGCCGCGGCCCTGAGCCAGGCGCTGCGCCTGGAGGCCCCAGAGCTTCCAGTAAAGGTGAACCCGACGAAGCCCCGGAGGGGCAGCTTCGAGGTGACGCTGCTGCGCCCGGACGGCAGCAGTGCGGAGCTCTGGACTGGGATTAAGAAGGGGCCCCCACGCAAACTCAAATTCCCTGAGCCTCAAGAGGTGGTGGAAGAGTTGAAGAAGTACCTGTCGTAGGGAGATTTGGGTAGAAGCCCTCATGCTGAGCTTTGTGTCCCTGGTGATGTTGGAACATTAATGATGGAACATGGCCAAACTTCAGTCATGATCCTGAAGCCATGGTTTCTTCCCTGCCAGAAATGAAGGTTCAGTTATGAGGCAACCCTCTAGTAAGGCATTGTAAAAGTTACTGGATTTGGTTTAATAAAAGTTGAAATAAAGTATTTGAGTAATAGAGTAAAACTTAATCTTAAGGGATGGACAAAGCCAGCCTTGTGGAGTTGGTAGTCCTTGTGTTTTTCATAGCCCAAGATAGCCATGTGACCTGGAACTTTTAAAAAAAAATTTTTTTTTAAGAGCCAGAGTCTCACTTTGTCACCCATGCTGGAGTATAGTGGTGCTATCTCGGCTCACTGCAACCATTGTCTCCCAGGTTCAAGCGATACTCCTGCCTCAGTCTCCTGAGTAGCTGGGATCACAGGTGCACACCACCACGCCTGGCTAATTTTTTTATTTTTATTTTTAGTAAAGGCGGGGTTTCACCATGTTGGCCAGGCTGGTCTCAAACTCCTGACCTCAAGTGATATGCCTGCCTTGACCTCCCAAAGTGCTGGGGTTACAGGCATGATGAGCCACCGCACCTGGCCTGACCTGCAACTTTTATGATCTCTGGTAAAACATGAGAAATACTTGTAAAAATTATGGGAGCAAAGTAAATAATTTATTTTTTAAAAAA"

# CODON TABLE
CDS= mRNA[105:474]
Genetic_code= {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",

    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",

    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",

    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"};

# TRANSLATION
protein = ""
for i in range(0, len(CDS), 3):
    codon = CDS[i:i+3]
    protein = protein + Genetic_code[codon]
print("The protein is:", protein)

# Cell 6
print("CDS length:", len(CDS))
print("First thirty residues:", protein[:30])
print("Total protein length:", len(protein))

# DEPOSITED PROTEIN
deposited_protein = "MAPRGRKRKAEAAVVAVAEKREKLANGGEGMEEATVVIEHCTSURVYGRNAAALSQALRLEAPELPVKVNPTKPRRGSFEVTLLRPDGSSAELWTGIKKGPPRKLKFPEPQEVVEELKKYLS"
print("len_protein:", len(protein))
print("len_deposited:", len(deposited_protein))

# DEPOSITED AND TRANSLATION COMPARISON
for i in range(len(protein)):
    if protein[i]!=deposited_protein[i]:
        print("the first difference:", i+1)
        print("The protein:", protein[i])
        print("The deposited protein:", deposited_protein[i])
        break
identical_y_n = "N"
mismatch ="44"
# =====================================================================
# SECTION 3: WRITE CONSOLIDATED SUMMARY TABLE CSV
# =====================================================================
headers = [
    "gene name", "category", "CDS length", "translated length", 
    "deposited length (from database)", "identical (Y/N)", "position of first mismatch"
]

row_drd2 = [
    "DRD2", "Nuclear", str(cds_len), str(len_amino_seq),
    str(depo_len), drd2_identical_OR_NOT, str(drd2_mismatch_pos)
]

row_mtcyb = [
    "MTCYB", "Mitochondrial", str(MTCYB_CDS_length), str(MTCYB_translated_length),
    str(MTCYB_deposited_length), str(MTCYB_is_identical_to_deposited_protein), str(MTCYB_mismatch_position)
]

row_selenoh = [
    "SELENOH", "Selenoprotein", str(len(CDS)), str(len(protein)),
    str(len(deposited_protein)),identical_y_n,mismatch  
]


# THIS CODE WRITES AND SAVES THE CSV FILE
output_filename = "Group5_table.csv"
with open(output_filename, "w") as file:
    file.write(",".join(headers) + "\n")
    file.write(",".join(row_drd2) + "\n")
    file.write(",".join(row_mtcyb) + "\n")
    file.write(",".join(row_selenoh) + "\n")

print("\nPipeline summary output saved to:", output_filename)
