# Force Calling Benchmark

We give a demo to help users of performing the regenotyping and evaluating performance. Here, we apply the integration VCF file from the callsets of HG002, HG003, and HG004 (all from GiaB Ashkenazim Trio) as the cohort-level SV target, then the HG002 alignments are selected to complete regenotyping on the SV target mentioned before.
The following procedures show the steps of applying regenotyping and reproducing the benchmark results. 

# Get tools

Information about how to install `conda` and add the `bioconda` channel is available on https://bioconda.github.io/.

```sh
conda create -n sniffles1_env python=3
conda activate sniffles1_env
conda install sniffles==1.0.12
```
```sh
conda create -n test_fc python=3
conda activate test_fc
conda install sniffles==2.0.2 cuteSV==2.0.2 svjedi==1.1.6 truvari==3.2.0 samtools tabix
# It will cost approximately 30 seconds to install cuteSV2.
```

# Get data
1) Create directory structure:
```sh
conda activate test_fc
mkdir -p ref alns tools/{sniffles1,sniffles2,cutesv,svjedi} giab
```

2) Download NIST and CMRG ground truth:
```sh
FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed > giab/HG002_SVs_Tier1_v0.6.bed
curl -s ${FTPDIR}/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz > giab/HG002_SVs_Tier1_v0.6.vcf.gz
```
```sh
FTPDIR=ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant  
curl -s ${FTPDIR}/HG002_GRCh37_CMRG_SV_v1.00.bed > giab/HG002_GRCh37_CMRG_SV_v1.00.bed
curl -s ${FTPDIR}/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz > giab/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz
```

3) Download hg19 reference with decoys and map non-ACGT characters to N:
```sh
curl -s ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz > ref/human_hs37d5.fasta.gz
gunzip ref/human_hs37d5.fasta.gz
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ref/human_hs37d5.fasta
```

4) Download all `.bam` files:
```sh
curl -s ftp://trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/alignment/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam > alns/HG002_origin.bam
samtools calmd -b alns/HG002_origin.bam ref/human_hs37d5.fasta > alns/HG002_all.bam
samtools index alns/HG002_all.bam
```

5) Download the trio vcf file:
```sh
curl -s https://zenodo.org/record/7347467/files/trio.vcf > trio.vcf
```

# Run sniffles1

6a) Run sniffles1 (v1.0.12):
```sh
conda activate sniffles1_env
sniffles -m alns/HG002_all.bam -v tools/sniffles1/sniffles1.call.vcf --Ivcf trio.vcf
conda deactivate
```
6b) Prepare for truvari:
```sh
grep '#' tools/sniffles1/sniffles1.call.vcf > tools/sniffles1/sniffles1.sort.vcf
grep -v '#' tools/sniffles1/sniffles1.call.vcf | sort -k 1,1 -k 2,2n >> tools/sniffles1/sniffles1.sort.vcf
grep '#' tools/sniffles1/sniffles1.sort.vcf > tools/sniffles1/sniffles1.vcf
grep -v '#' tools/sniffles1/sniffles1.sort.vcf | grep -v '0/0' | grep -v "\./\." >> tools/sniffles1/sniffles1.vcf
bgzip -c tools/sniffles1/sniffles1.vcf > tools/sniffles1/sniffles1.vcf.gz
tabix tools/sniffles1/sniffles1.vcf.gz
```

# Run sniffles2

7a) Run sniffles2 (v2.0.2):
```sh
sniffles --input alns/HG002_all.bam --vcf tools/sniffles2/sniffles2.call.vcf --genotype-vcf trio.vcf
```
7b) Prepare for truvari:
```sh
grep '#' tools/sniffles2/sniffles2.call.vcf > tools/sniffles2/sniffles2.sort.vcf
grep -v '#' tools/sniffles2/sniffles2.call.vcf | sort -k 1,1 -k 2,2n >> tools/sniffles2/sniffles2.sort.vcf
awk -F '\t' '{if($1=="#CHROM") {for(i=1;i<10;i++) printf($i"\t"); print($10);} else print($0);}' tools/sniffles2/sniffles2.sort.vcf > temp.vcf
sed -i 'N;122 a ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# High-quality variant reads">' temp.vcf
sed -i 'N;122 a ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="# Genotype quality">' temp.vcf
grep '#' temp.vcf > tools/sniffles2/sniffles2.vcf
grep -v '#' temp.vcf | grep -v '0/0' | grep -v "\./\." >> tools/sniffle2/sniffles2.vcf
rm temp.vcf
bgzip -c tools/sniffles2/sniffles2.vcf > tools/sniffles2/sniffles2.vcf.gz
tabix tools/sniffles2/sniffles2.vcf.gz
```

# Run cuteSV2

8a) Run cuteSV2 (v2.0.2):
```sh
cuteSV alns/HG002_all.bam ref/human_hs37d5.fasta tools/cutesv/cutesv.call.vcf ./ --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.5 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -Ivcf trio.vcf -q 10
```
8b) Prepare for truvari:
```sh
grep '#' tools/cutesv/cutesv.call.vcf > tools/cutesv/cutesv.vcf
grep -v '#' tools/cutesv/cutesv.call.vcf | grep -v '0/0' | grep -v "\./\." >> tools/cutesv/cutesv.vcf
bgzip -c tools/cutesv/cutesv.vcf > tools/cutesv/cutesv.vcf.gz
tabix tools/cutesv/cutesv.vcf.gz
# The sample output of cuteSV2 is available at https://doi.org/10.5281/zenodo.7347467.
```

# Run SVJedi

9a) Run SVJedi (v1.1.6):
```sh
samtools fasta alns/HG002_all.bam > alns/HG002_all.fasta
python3 svjedi.py -v trio.vcf -r ref/human_hs37d5.fasta -i alns/HG002_all.fasta -o tools/svjedi/svjedi.call.vcf
```
9b) Prepare for truvari:
```sh
grep '#' tools/svjedi/svjedi.call.vcf > tools/svjedi/svjedi.sort.vcf
grep -v '#' tools/svjedi/svjedi.call.vcf | sort -k 1,1 -k 2,2n >> tools/svjedi/svjedi.sort.vcf
grep '#' tools/svjedi/svjedi.sort.vcf > tools/svjedi/svjedi.vcf
grep -v '#' tools/svjedi/svjedi.sort.vcf | grep -v '0/0' | grep -v "\./\." >> tools/svjedi/svjedi.vcf
bgzip -c tools/svjedi/svjedi.vcf > tools/svjedi/svjedi.vcf.gz
tabix tools/svjedi/svjedi.vcf.gz
```

# Final comparison

10a) Compare to NIST ground truth (v3.2.0):
```sh
truvari bench -b giab/HG002_SVs_Tier1_v0.6.vcf.gz -c tools/sniffles1/sniffles1.vcf.gz\
        --includebed giab/HG002_SVs_Tier1_v0.6.bed -o NIST-sniffles1 -p 0 -r 1000 --multimatch --passonly
truvari bench -b giab/HG002_SVs_Tier1_v0.6.vcf.gz -c tools/sniffles2/sniffles2.vcf.gz\
        --includebed giab/HG002_SVs_Tier1_v0.6.bed -o NIST-sniffles2 -p 0 -r 1000 --multimatch --passonly
truvari bench -b giab/HG002_SVs_Tier1_v0.6.vcf.gz -c tools/cutesv/cutesv.vcf.gz\
        --includebed giab/HG002_SVs_Tier1_v0.6.bed -o NIST-cutesv -p 0 -r 1000 --multimatch --passonly
truvari bench -b giab/HG002_SVs_Tier1_v0.6.vcf.gz -c tools/svjedi/svjedi.vcf.gz\
        --includebed giab/HG002_SVs_Tier1_v0.6.bed -o NIST-svjedi -p 0 -r 1000 --multimatch --passonly
```
10b) Compare to CMRG ground truth (v3.2.0):
```sh
truvari bench -b giab/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz -c tools/sniffles1/sniffles1.vcf.gz\
        --includebed giab/HG002_GRCh37_CMRG_SV_v1.00.bed -o CMRG-sniffles1 -p 0 -r 1000 --multimatch --passonly
truvari bench -b giab/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz -c tools/sniffles2/sniffles2.vcf.gz\
        --includebed giab/HG002_GRCh37_CMRG_SV_v1.00.bed -o CMRG-sniffles2 -p 0 -r 1000 --multimatch --passonly
truvari bench -b giab/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz -c tools/cutesv/cutesv.vcf.gz\
        --includebed giab/HG002_GRCh37_CMRG_SV_v1.00.bed -o CMRG-cutesv -p 0 -r 1000 --multimatch --passonly
truvari bench -b giab/HG002_GRCh37_CMRG_SV_v1.00.vcf.gz -c tools/svjedi/svjedi.vcf.gz\
        --includebed giab/HG002_GRCh37_CMRG_SV_v1.00.bed -o CMRG-svjedi -p 0 -r 1000 --multimatch --passonly
```

# Down sample

11) Downsample the original alignment file:
```sh
samtools view -h -s 0.66 alns/HG002_all.bam | samtools view -bS > alns/HG002_20x.bam
samtools view -h -s 0.33 alns/HG002_all.bam | samtools view -bS > alns/HG002_10x.bam
samtools view -h -s 0.17 alns/HG002_all.bam | samtools view -bS > alns/HG002_5x.bam
```
