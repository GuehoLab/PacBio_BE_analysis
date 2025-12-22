#!/bin/bash

minimap2 -x map-hifi ref.fasta example_hifi_reads.fastq.gz > example_hifi_reads.paf
PacBe_demultiplex example_hifi_reads.fastq.gz example_hifi_reads.paf index_table.txt demulti/


mkdir -p lungo_output 

CRISPRlungo ref.fasta demulti/DIP2B_sg692_10x.fastq --control demulti/DIP2B_empty_10x.fastq lungo_output/DIP2B_sg692_10x CGGATCCTGTAGCCGGGTGT --using_all_mutations --window 10000 -t 20

PacBe_organize_data lungo_output/DIP2B_sg692_10x/ CGGATCCTGTAGCCGGGTGT