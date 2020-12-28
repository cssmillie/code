#!/bin/bash

f=$1 # Input FASTA file

# If FASTQ, convert to FASTA and blast against UniVec
if echo $f | grep -q -e fsq -e fastq; then
    #python ~/box/fastq_to_fasta.py $f | blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 -db ~/aviv/db/univec/univec.fst > test.blastn
    #cat test.blastn | cut -f 1 | sort -u > test.remove
    #cat test.remove | python ~/box/fasta_remove.py --fsq $f > test.fastq
    python ~/box/fastq_to_fasta.py $f | blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 -db ~/aviv/db/univec/univec.fst | cut -f 1 | sort -u | python ~/box/fasta_remove.py --fsq $f
elif echo $f | grep -q -e fst -e fasta; then
    cat $f | blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 -db ~/aviv/db/univec/univec.fst | cut -f 1 | sort -u | python ~/box/fastq_remove.py --fsq $f
else
    echo "Error: file suffix not fst, fasta, fsq, or fastq"
fi
