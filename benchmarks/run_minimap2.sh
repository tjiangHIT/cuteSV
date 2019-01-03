#!/usr/bin/bash 

Thaliana_ref=/home/tjiang/Trio/Thaliana/GCF_000001735.4_TAIR10.1_genomic.fa
Thaliana_p1=/data2/tjiang/benchmarks/Arabidopsis_trio/col.fa
Thaliana_p2=/data2/tjiang/benchmarks/Arabidopsis_trio/cvi.fa
Thaliana_child=/data2/tjiang/benchmarks/Arabidopsis_trio/f1.fa

ref=${}
# not include .fa
query=${}
temp=${}
# minimap2
minimap2 -ax map-pb ${ref} ${query}.fa -t 16 --MD -Y > ${query}_minimap2.sam
samtools view -Sb ${query}_minimap2.sam | samtools sort -O bam -T ${temp} - > ${query}_minimap2.bam && samtools index ${query}_minimap2.bam

# sniffles
sniffles -m ${query}_minimap2.bam -v ${query}_minimap2_sniffles.vcf -s 10 -l 50 --genotype --report_BND

# pbsv
pbsv discover ${query}_minimap2.bam ${query}_minimap2_pbsv.svsig.gz -s ${query}
pbsv call ${ref} ${query}_minimap2_pbsv.svsig.gz ${query}_minimap2_pbsv.vcf -j 16 -m 50

# SVIM
/home/tjiang/miniconda3/bin/python3 /home/tjiang/Tools/svim/SVIM.py alignment --min_sv_size 50 svim/f1/ f1_minimap2.bam
