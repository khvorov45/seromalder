mafft --auto summary-measures/cdc-fasta.fasta > summary-measures/cdc-fasta-aligned.fasta
seqkit fx2tab summary-measures/cdc-fasta-aligned.fasta > summary-measures/cdc-fasta-aligned.tsv
