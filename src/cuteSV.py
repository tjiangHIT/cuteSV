
import sys
import logging
import time

from cuteSV_utils import setupLogging, __version__, __author__, __email__
from cuteSV_workflow import wrk_flow

USAGE = '''
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

	cuteSV is a novel, user-friendly tool specifically designed for \
	the needs of the SV detection with long reads. Totally, it cotains
	3 steps for generation final SV call set.
	
	[STEP 1] Extract signatures of SV in each reads and store them into 
		files with similar physical locations and same Signature-type.

	[STEP 2] Clean and filter low-quality signatures using diversity and 
		density of signal distribution.

	[STEP 3] Divide chromosome into several parts, and then cluster each parts. 	 
	
	[STEP 4] Merge each clusters into super-cluster if they satisfy the 
		conditions that close to each other and other thresholds.

	[STEP 5] Replenish genotype info and write into hard-disk in VCF format.

	cuteSV V%s
	Auther %s
	Contact %s
'''%(__version__, __author__, __email__)

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="cuteSV", description=USAGE, 
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input", metavar="[BAM]", type=str, 
		help="Sorted .bam file form NGMLR or Minimap2.")
	parser.add_argument('output', type=str, help = "the path of [Output]")
	parser.add_argument('temp_dir', type=str, help = "temporary directory to use for distributed jobs")
	parser.add_argument('-s', '--min_support', 
		help = "Minimum number of reads that support a SV to be reported.[%(default)s]", 
		default = 5, type = int)
	parser.add_argument('-l', '--min_length', help = "Minimum length of SV to be reported.[%(default)s]", 
		default = 50, type = int)
	parser.add_argument('-p', '--max_split_parts', 
		help = "Maximum number of split segments a read may be aligned before it is ignored.[%(default)s]", 
		default = 7, type = int)
	parser.add_argument('-q', '--min_mapq', 
		help = "Minimum mapping quality value of alignment to be taken into account.[%(default)s]", 
		default = 20, type = int)
	parser.add_argument('-d', '--max_distance', help = "Maximum distance to group SV together..[%(default)s]", 
		default = 1000, type = int)
	parser.add_argument('-r', '--min_seq_size', 
		help = "Ignores reads that only report alignments with not longer then bp.[%(default)s]", 
		default = 2000, type = int)
	parser.add_argument('-t', '--threads', help = "Number of threads to use.[%(default)s]", 
		default = 16, type = int)
	parser.add_argument('-b', '--batches', help = "A batches of reads to load.[%(default)s]", 
		default = 10000000, type = int)
	parser.add_argument('-c', '--max_cluster_bias', 
		help = "Maximum distance to cluster read together.[%(default)s]", default = 50, type = int)
	args = parser.parse_args(argv)
	return args

def run(argv):
	'''run cuteSV workflow and generate SV call set'''
	args = parseArgs(argv)
	setupLogging(False)
	starttime = time.time()
	wrk_flow(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	run(sys.argv[1:])