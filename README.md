# cuteSV

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

	
	$ git clone https://github.com/tjiangHIT/cuteSV.git
	$ cd cuteSV/
	$ pip install .

---	
### Introduction
Long-read sequencing technologies enable to comprehensively discover structural variations (SVs). However, it is still non-trivial for state-of-the-art approaches to detect SVs with high sensitivity or high performance or both. Herein, we propose cuteSV, a sensitive, fast and lightweight SV detection approach. cuteSV uses tailored methods to comprehensively collect various types of SV signatures, and a clustering-and-refinement method to implement a stepwise SV detection, which enables to achieve high sensitivity without loss of accuracy. Benchmark results demonstrate that cuteSV has better yields on real datasets. Further, its speed and scalability are outstanding and promising to large-scale data analysis.

---
### Dependences
	
	1. python
	2. pysam
	3. Biopython
	4. cigar
	5. numpy

---
### Synopsis
	python cuteSV.py <alignment> <output.vcf> <work_dir>
	
---
### Contact
For advising, bug reporting and requiring help, please post on [Github Issue](https://github.com/tjiangHIT/cuteSV/issues) or contact tjiang@hit.edu.cn.
