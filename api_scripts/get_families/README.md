

Script: get_families

The script goes through the list of families in 'families.txt' and fetches
for each the HMM, alignment, sequences, and the tree.
It saves it into 'treefam_data' and names them with the family_name, e.g. 

alignment (fasta): TF101001.aln
sequences (protein): TF101001.fa
sequences (cds): TF101001.cds.fa
tree: TF101001.nhx

Input file: families.txt
Output directory: treefam_data

usage: get_families.pl -all -hmm -seq -tree -aln
       -all: gets hmm, alignment, sequences (protein+cds), tree for each family (default option) 
       -hmm:  gets hmm for each family
       -seq: gets sequences (protein+cds) for each family
       -tree: gets tree for each family
