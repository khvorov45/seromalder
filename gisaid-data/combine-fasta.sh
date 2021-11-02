cat gisaid-data/sequences*.fasta > gisaid-data/combined.fasta
seqkit fx2tab gisaid-data/combined.fasta > gisaid-data/combined-fasta.tsv
