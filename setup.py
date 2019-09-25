# coding=utf-8

from setuptools import setup, find_packages

LONG_DESCRIPTION = '''Long-read sequencing technologies enable to comprehensively \
discover structural variations (SVs). However, it is still non-trivial for \
state-of-the-art approaches to detect SVs with high sensitivity or high performance \
or both. Herein, we propose cuteSV, a sensitive, fast and lightweight SV detection \
approach. cuteSV uses tailored methods to comprehensively collect various types of \
SV signatures, and a clustering-and-refinement method to implement a stepwise SV \
detection, which enables to achieve high sensitivity without loss of accuracy. \
Benchmark results demonstrate that cuteSV has better yields on real datasets. \
Further, its speed and scalability are outstanding and promising to large-scale data analysis.'''

setup(
    name = "cuteSV",
    version = "1.0.1",
    description = "Long read based human genomic structural variation detection with cuteSV",
    author = "Jiang Tao",
    author_email = "tjiang@hit.edu.cn",
    url = "https://github.com/tjiangHIT/cuteSV",
    license = "MIT",
    packages = find_packages("src"),
    package_dir = {"": "src"},
    data_files = [("", ["LICENSE"])],
    scripts=['src/cuteSV/cuteSV'],
    long_description = LONG_DESCRIPTION,
    zip_safe = False,
    install_requires = ['pysam', 'Biopython', 'Cigar', 'numpy']
)
