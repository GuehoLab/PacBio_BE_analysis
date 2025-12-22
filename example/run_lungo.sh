#!/bin/bash

mkdir -p lungo_output 

CRISPRlungo ref.fasta demulti/DIP2B_sg692_10x.fastq --control demulti/DIP2B_empty_10x.fastq lungo_output/DIP2B_sg692_10x CGGATCCTGTAGCCGGGTGT --using_all_mutations --window 10000 -t 20

python3 ../organize_data.py lungo_output/DIP2B_sg692_10x/ CGGATCCTGTAGCCGGGTGT