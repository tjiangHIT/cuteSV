# cuteSV

[![PyPI version](https://badge.fury.io/py/cuteSV.svg)](https://badge.fury.io/py/cuteSV)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cutesv/badges/version.svg)](https://anaconda.org/bioconda/cutesv)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cutesv/badges/license.svg)](https://anaconda.org/bioconda/cutesv)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cutesv/badges/platforms.svg)](https://anaconda.org/bioconda/cutesv)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/cutesv/badges/latest_release_date.svg)](https://anaconda.org/bioconda/cutesv)

---
### Getting Start
	                                               __________    __       __
	                                              |   ____   |  |  |     |  |
	                          _                   |  |    |__|  |  |     |  |
	 _______    _     _   ___| |___     ______    |  |          |  |     |  |
	|  ___  |  | |   | | |___   ___|   / ____ \   |  |_______   |  |     |  |
	| |   |_|  | |   | |     | |      / /____\ \  |_______   |  |  |     |  |
	| |        | |   | |     | |      | _______|   __     |  |  \  \     /  /
	| |    _   | |   | |     | |  _   | |     _   |  |    |  |   \  \   /  /
	| |___| |  | |___| |     | |_| |  \ \____/ |  |  |____|  |    \  \_/  /
	|_______|  |_______|     |_____|   \______/   |__________|     \_____/


---	
### Installation
	$ pip install cuteSV
	or
	$ conda install -c bioconda cutesv
	or 
	$ git clone https://github.com/tjiangHIT/cuteSV.git && cd cuteSV/ && python setup.py install 

---	
### Introduction
Long-read sequencing enables the comprehensive discovery of structural variations (SVs). However, it is still non-trivial to achieve high sensitivity and performance simultaneously due to the complex SV characteristics implied by noisy long reads. Therefore, we propose cuteSV, a sensitive, fast and scalable long-read-based SV detection approach. cuteSV uses tailored methods to collect the signatures of various types of SVs and employs a clustering-and-refinement method to analyze the signatures to implement sensitive SV detection. Benchmarks on real Pacific Biosciences (PacBio) and Oxford Nanopore Technology (ONT) datasets demonstrate that cuteSV has better yields and scalability than state-of-the-art tools.

The benchmark results of cuteSV on the HG002 human sample are below:

<img src="hg2_assessments.png" width="800px">

BTW, we used [Truvari](https://github.com/spiralgenetics/truvari) to calculate the recall, precision, and f-measure.
For more detailed implementation of SV benchmarks, we show an example [here](https://github.com/tjiangHIT/sv-benchmark).

---
### Dependence
	
	1. python3
	2. pysam
	3. Biopython
	4. cigar
	5. numpy

---
### Usage
	cuteSV <sorted.bam> <output.vcf> <work_dir>
	
*Suggestions*

	> For PacBio CLR/ONT data:
		--max_cluster_bias_INS		100
		--diff_ratio_merging_INS	0.2
		--diff_ratio_filtering_INS	0.6
		--diff_ratio_filtering_DEL	0.7
	> For PacBio CCS(HIFI) data:
		--max_cluster_bias_INS		200
		--diff_ratio_merging_INS	0.65
		--diff_ratio_filtering_INS	0.65
		--diff_ratio_filtering_DEL	0.35
	
| Parameter | Description | Default |
| :------------ |:---------------|-------------:|
|--threads|Number of threads to use.| 16 |
|--batches| Batch of genome segmentation interval.|10,000,000|
|--sample| Sample name/id |NULL|
|--max_split_parts|Maximum number of split segments a read may be aligned before it is ignored.|7|
|--min_mapq|Minimum mapping quality value of alignment to be taken into account.|20|
|--min_read_len|Ignores reads that only report alignments with not longer than bp.|500|
|--merge_del_threshold|Maximum distance of deletion signals to be merged.|0|
|--merge_ins_threshold|Maximum distance of insertion signals to be merged.|100|
|--min_support|Minimum number of reads that support a SV to be reported.|10|
|--min_size|Minimum length of SV to be reported.|30|
|--max_size|Minimum length of SV to be reported.|100000|
|--genotype|Enable to generate genotypes.|False|
|--hom|Threshold on allele frequency for homozygous.|0.8|
|--het|Threshold on allele frequency for heterozygous.|0.2|
|--max_cluster_bias_INS|Maximum distance to cluster read together for insertion.|100|
|--diff_ratio_merging_INS|Do not merge breakpoints with basepair identity more than the ratio of *default* for insertion.|0.2|
|--diff_ratio_filtering_INS|Filter breakpoints with basepair identity less than the ratio of *default* for insertion.|0.6|
|--max_cluster_bias_DEL|Maximum distance to cluster read together for deletion.|200|
|--diff_ratio_merging_DEL|Do not merge breakpoints with basepair identity more than the ratio of *default* for deletion.|0.3|
|--diff_ratio_filtering_DEL|Filter breakpoints with basepair identity less than the ratio of *default* for deletion.|0.7|
|--max_cluster_bias_INV|Maximum distance to cluster read together for inversion.|500|
|--max_cluster_bias_DUP|Maximum distance to cluster read together for duplication.|500|
|--max_cluster_bias_TRA|Maximum distance to cluster read together for translocation.|50|
|--diff_ratio_filtering_TRA|Filter breakpoints with basepair identity less than the ratio of *default* for translocation.|0.6|

---
### Datasets generated from cuteSV
We provided the SV callsets of the HG002 human sample produced by cuteSV form three different long-read sequencing platforms (i.e. PacBio CLR, PacBio CCS, and ONT PromethION). 

You can download them at: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3556403.svg)](https://doi.org/10.5281/zenodo.3556403)

Please cite the manuscript of cuteSV before using these callsets.

---
### Changelog

	cuteSV (v1.0.5):
	1.Add new options for specificly setting the threshold of deletion/insertion signals merging in the same read. The default parameters are 0 bp for deletion and 100 bp for insertion.
	2.Remove parameter --merge_threshold.
	3.Fix bugs in inversion and translocation calling.
	4.Add new option for specificly setting the maximum size of SV to be discovered. The default value is 100,000 bp. 


	cuteSV (v1.0.4):
	1.Add a new option for specificly setting the threshold of SV signals merging in the same read. The default parameter is 500 bp. You can reduce it for high-quality sequencing datasets like PacBio HiFi (CCS).
	2.Make the genotyping function optional.
	3.Enable users to set the threshold of SV allele frequency of homozygous/heterozygous.
	4.Update the description of recommendation parameters in processing ONT data.

	cuteSV (v1.0.3):
	1.Refine the genotyping model.
	2.Adjust the threshold value of heterozygosis alleles.

	cuteSV (v1.0.2):
	1.Improve the genotyping performance and enable it to be default option.
	2.Make the description of parameters better.
	3.Modify the header description of vcf file.
	4.Add two new indicators, i.e., BREAKPOINT_STD and SVLEN_STD, to further characterise deletion and insertion.
	5.Remove a few redundant functions which will reduce code readability.


---
### Citation
Long Read based Human Genomic Structural Variation Detection with cuteSV. Tao Jiang, et al. bioRxiv 780700; doi: https://doi.org/10.1101/780700
	
---
### Contact
For advising, bug reporting and requiring help, please post on [Github Issue](https://github.com/tjiangHIT/cuteSV/issues) or contact tjiang@hit.edu.cn.
