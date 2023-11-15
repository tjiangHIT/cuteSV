# Force Calling Benchmark

We provide a demo to help users of performing the regenotyping and evaluating performance. Here, we regard the "High confidence variant" of HG002 from the Genome in a Bottle (GiaB) ground truth set (SV v0.6) from the National Institute of Standards and Technology (NIST) as the target SVs. Then the HG002 alignments from PacBio HiFi sequencing are selected to complete regenotyping on the above target SVs.
The following procedures show the steps of applying regenotyping and reproducing the benchmark results. 

# Get tools

Information about how to install `conda` and add the `bioconda` channel is available on https://bioconda.github.io/.

```sh
conda create -n sniffles1_env python=3.10
conda activate sniffles1_env
conda install -c bioconda sniffles==1.0.12
```
```sh
conda create -n test_fc python=3.10
conda activate test_fc
conda install -c bioconda sniffles==2.0.7 cuteSV==2.1.0 svjedi==1.1.6 truvari==3.5.0 samtools tabix
```

# Get data
1) Create directory structure:
```sh
conda activate test_fc
mkdir -p ref alns tools giab
```

2) Download NIST ground truth:
```sh
FTPDIR=https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed > giab/HG002_SVs_Tier1_v0.6.bed
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz > giab/HG002_SVs_Tier1_v0.6.vcf.gz
```

3) Download hg19 reference with decoys and map non-ACGT characters to N:
```sh
curl -s https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > ref/human_hs37d5.fasta.gz
gunzip ref/human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ref/human_hs37d5.fasta
```

4) Download the alignment files:
```sh
curl -s https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam > alns/HG002.m84011_220902_175841_s1.GRCh38.bam
samtools fasta -@ 16 alns/HG002.m84011_220902_175841_s1.GRCh38.bam > alns/HG002.fasta
pbmm2 align --num-threads 16 --preset CCS --sample HG002 --log-level INFO --sort --unmapped -c 0 -y 70 ref/human_hs37d5.fasta alns/HG002.fasta alns/HG002_origin.bam
samtools calmd -b alns/HG002_origin.bam ref/human_hs37d5.fasta > alns/HG002_all.bam
samtools index alns/HG002_all.bam
```

5) Generate target SV files:
```sh
gzip -d giab/HG002_SVs_Tier1_v0.6.vcf.gz
grep '#' giab/HG002_SVs_Tier1_v0.6.vcf > giab/HG002_SVs_Tier1_v0.6.filter.vcf
grep -v '#' giab/HG002_SVs_Tier1_v0.6.vcf | awk -F '\t' '{split($10,X,":"); if(X[1]!="0/0"&&X[1]!="./.") print $0}' >> giab/HG002_SVs_Tier1_v0.6.filter.vcf
bgzip -c giab/HG002_SVs_Tier1_v0.6.filter.vcf > giab/HG002_SVs_Tier1_v0.6.filter.vcf.gz
tabix giab/HG002_SVs_Tier1_v0.6.filter.vcf.gz
```

# Run Sniffles1

6a) Run Sniffles1 (v1.0.12):
```sh
conda activate sniffles1_env
sniffles -m alns/HG002_all.bam -v tools/sniffles1.call.vcf --Ivcf giab/HG002_SVs_Tier1_v0.6.vcf
conda deactivate
```
6b) Prepare for truvari:
```sh
grep '#' tools/sniffles1.call.vcf > tools/sniffles1.sort.vcf
sed -i 'N;96 a ##FILTER=<ID=STRANDBIAS,Description="STRANDBIAS.">' tools/sniffles1.sort.vcf
grep -v '#' tools/sniffles1.call.vcf | sort -k 1,1 -k 2,2n | awk -F '\t' '{split($10,X,":"); if(X[1]!="0/0"&&X[1]!="./.") print $0}' >> tools/sniffles1.sort.vcf
bgzip -c tools/sniffles1.sort.vcf > tools/sniffles1.sort.vcf.gz
tabix tools/sniffles1.sort.vcf.gz
```

# Run Sniffles2

7a) Run Sniffles2 (v2.0.7):
```sh
sniffles --input alns/HG002_all.bam --vcf tools/sniffles2.call.vcf --genotype-vcf giab/HG002_SVs_Tier1_v0.6.vcf
```
7b) Prepare for truvari:
```sh
grep '#' tools/sniffles2.call.vcf > tools/sniffles2.sort.vcf
sed -i 'N;104 a ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# High-quality variant reads">' tools/sniffles2.sort.vcf
sed -i 'N;104 a ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="# Genotype quality">' tools/sniffles2.sort.vcf
sed -i 'N;104 a ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# High-quality reference reads">' tools/sniffles2.sort.vcf
grep -v '#' tools/sniffles2.call.vcf | sort -k 1,1 -k 2,2n | awk -F '\t' '{split($10,X,":"); if(X[1]!="0/0"&&X[1]!="./.") print $0}' >> tools/sniffles2.sort.vcf
bgzip -c tools/sniffles2.sort.vcf > tools/sniffles2.sort.vcf.gz
tabix tools/sniffles2.sort.vcf.gz
```

# Run cuteSV2

8a) Run cuteSV2 (v2.1.0):
```sh
cuteSV alns/HG002_all.bam ref/human_hs37d5.fasta tools/cutesv.call.vcf ./ --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -Ivcf giab/HG002_SVs_Tier1_v0.6.vcf -q 10 -L -1
```
8b) Prepare for truvari:
```sh
grep '#' tools/cutesv.call.vcf > tools/cutesv.sort.vcf
grep -v '#' tools/cutesv.call.vcf | sort -k 1,1 -k 2,2n | awk -F '\t' '{split($10,X,":"); if(X[1]!="0/0"&&X[1]!="./.") print $0}' >> tools/cutesv.sort.vcf
bgzip -c tools/cutesv.sort.vcf > tools/cutesv.sort.vcf.gz
tabix tools/cutesv.sort.vcf.gz
```

# Run SVJedi

9a) Run SVJedi (v1.1.6):
```sh
python3 svjedi.py -v giab/HG002_SVs_Tier1_v0.6.vcf -r ref/human_hs37d5.fasta -i alns/HG002.fasta -o tools/svjedi.call.vcf
```
9b) Prepare for truvari:
```sh
grep '#' tools/svjedi.call.vcf > tools/svjedi.sort.vcf
grep -v '#' tools/svjedi.call.vcf | sort -k 1,1 -k 2,2n | awk -F '\t' '{split($10,X,":"); if(X[1]!="0/0"&&X[1]!="./.") print $0}' >> tools/svjedi.sort.vcf
bgzip -c tools/svjedi.sort.vcf > tools/svjedi.sort.vcf.gz
tabix tools/svjedi.sort.vcf.gz
```

# Final comparison

10) Compare to NIST ground truth via truvari (v3.5.0):
```sh
for tools in {sniffles1, sniffles2, cutesv, svjedi}
do
        truvari bench -b giab/HG002_SVs_Tier1_v0.6.filter.vcf.gz -c tools/$i.sort.vcf.gz --includebed giab/HG002_SVs_Tier1_v0.6.bed -o cmp -p 0 -r 2 -P 1 --sizemax 1000000
done
```

# Down sample

11) Downsample the original alignment file:
```sh
samtools view -h -s 0.66 alns/HG002_all.bam | samtools view -bS > alns/HG002_20x.bam
samtools view -h -s 0.33 alns/HG002_all.bam | samtools view -bS > alns/HG002_10x.bam
samtools view -h -s 0.17 alns/HG002_all.bam | samtools view -bS > alns/HG002_5x.bam
```
