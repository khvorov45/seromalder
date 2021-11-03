mafft --auto virus-distances/cdc-fasta.fasta > virus-distances/cdc-fasta-aligned.fasta
seqkit fx2tab virus-distances/cdc-fasta-aligned.fasta > virus-distances/cdc-fasta-aligned.tsv
seqkit translate virus-distances/cdc-fasta-aligned.fasta > virus-distances/cdc-fasta-aligned-aa.fasta
seqkit fx2tab virus-distances/cdc-fasta-aligned-aa.fasta > virus-distances/cdc-fasta-aligned-aa.tsv
