## Description:

Extracts the promotor sequence (default defined as 1000 bp upstream sequence, can be changed via promotorregion option) for the given genes from the genome fasta file.
It uses in the background the bedtools command line tool. An example output can be found under the /example_data/Promotor.fasta.

# Usage:

python BEDT.py bedfile genelist genomefasta intermediatefolder shscript --promotorregion

## with example data:

python ../BEDT.py CH.bed ./genelist_query_x.txt ./CH.fasta ./test42 ../bedtools.sh
