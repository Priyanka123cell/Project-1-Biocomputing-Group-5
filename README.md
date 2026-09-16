# Project-1-Biocomputing-Group-5
Steps performed:
1. Downloading FASTA sequence and genbank data from NCBI and genbank respectively.
2. Extracting CDS from mRNA (using slicing in python) for nuclear and the gene coding for selenoprotein. For MTCYB CDS was extracted from the mitochondrial genome and two AA residues were added to the 3' end of the sequence. Divided the CDS by 3 to check if it gives whole number or not.
3. Made codon table using dictionary (in python)
4. Slicing the codons(3 bases) from the CDS and stored in a variable. When it is called using dictionary as a key, it returns its value which is the amino acid it codes for. This is the translation of the coding sequence.
5. Comparing the translated sequence with the deposited sequence and report their lengths and if mismatch is there then report that. Also compared translated sequence with Biopython trnaslated sequence and report if any difference is present or not.Used zip for comparison.
6. Buiding a pipeline to report all three genes and writing a single summary table containing, per gene: gene name, category, CDS length, translated length, deposited length (from database),identical (Y/N) and position of first mismatch. Script written in nano (project_Q4_G5.py).
## Project Information

| Question | Description | File |
|----------|-------------|------|
| 1(a) | MT-CYB sequence search and FASTA record | MTCYB_protein.fasta |
| 1(b) | GenBank record and CDS information | Mitochondrial_genome.gb |
| 1(c) | Deposited protein sequence | deposited_protein.fasta |
| 1(d) | Header lines | headers.txt |
| 2(a) | CDS extraction and length check | MTCYB_CDS_corrected.fasta |
| 2(b) | Translation of CDS | Sparshika_MTCYB_Script.ipynb |
# Biocomputing Project 1

## 1. Project Overview

This project involved retrieving a selected gene sequence from NCBI, obtaining its GenBank and protein records, extracting the CDS, translating the coding sequence, and comparing the translated protein sequence with the deposited protein sequence.

## 2. Team Members and Gene Selection

| Member | Category | Gene |
|--------|----------|------|
| Sparshika Mishra | Mitochondrial gene | MT-CYB |
| Sampurna Ganguly | Selenoprotein gene | SELENOH |
| Priyanka Deb | Nuclear gene | DRD2 |
