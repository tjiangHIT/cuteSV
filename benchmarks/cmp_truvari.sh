#!/usr/bin/bash 

INVCF=$1
outdir=$2
EXP=/data2/tjiang/benchmarks/GiaB_trio/HG002_SVs_Tier1_v0.6.vcf.gz
HIGHBED=/data2/tjiang/benchmarks/GiaB_trio/HG002_SVs_Tier1_v0.6.bed
REF=/home/tjiang//Ref/hs37d5/hs37d5.fa

grep '^#' ${INVCF} > ${INVCF}.sorted.vcf && grep -v '^#' ${INVCF} | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> ${INVCF}.sorted.vcf
bgzip -c ${INVCF}.sorted.vcf > ${INVCF}.sorted.vcf.gz
tabix -p vcf ${INVCF}.sorted.vcf.gz

truvari.py -b ${EXP} -c ${INVCF}.sorted.vcf.gz -o ${outdir} -f ${REF} --giabreport --passonly --includebed ${HIGHBED} -r 1000 -p 0
