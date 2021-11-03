mafft --auto summary-measures/cdc-fasta.fasta > summary-measures/cdc-fasta-aligned.fasta
seqkit fx2tab summary-measures/cdc-fasta-aligned.fasta > summary-measures/cdc-fasta-aligned.tsv
seqkit translate summary-measures/cdc-fasta-aligned.fasta > summary-measures/cdc-fasta-aligned-aa.fasta
seqkit fx2tab summary-measures/cdc-fasta-aligned-aa.fasta > summary-measures/cdc-fasta-aligned-aa.tsv
