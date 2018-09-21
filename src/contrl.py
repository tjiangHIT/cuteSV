import argparse
import sys
import logging
import pysam
import cigar
from multiprocessing import Pool
from __utils__ import *
import time, gc

MainCandidate = dict()
MainCandidate["DEL"] = dict()
MainCandidate["INS"] = dict()
MainCandidate["INV"] = dict()
MainCandidate["DUP"] = dict()
MainCandidate["TRA"] = dict()
MainCandidate["l_DUP"] = dict()
MainCandidate["DUP_r"] = dict()
MainCandidate["l_INV"] = dict()
MainCandidate["INV_r"] = dict()

def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
    logging.info("Running %s" % " ".join(sys.argv))

def single_pipe(sam_path, Chr_name, min_length, min_mapq, max_split_parts, min_seq_size, total_reads):
	samfile = pysam.AlignmentFile(sam_path)
	Chr_length = samfile.get_reference_length(Chr_name)
	# print Chr_name, Chr_length
	logging.info("Resolving the chromsome %s."%(Chr_name))
	logging.info("%d reads on %s."%(total_reads, Chr_name))

	for read in samfile.fetch(Chr_name):
		parse_read(read, Chr_name, min_length, min_mapq, max_split_parts, min_seq_size)
	samfile.close()
	return candidate
	# print MainCandidate

def multi_run_wrapper(args):
	return single_pipe(*args)

def main_ctrl(args):
	# starttime = time.time()
	samfile = pysam.AlignmentFile(args.input)

	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:x-[1])
	analysis_pools = Pool(processes=int(args.threads))

	result = list()
	for i in process_list:
		para = [(args.input, i[0], args.min_length, args.min_mapq, args.max_split_parts, args.min_seq_size, i[1])]
		result.append(analysis_pools.map_async(multi_run_wrapper, para))
	analysis_pools.close()
	analysis_pools.join()

	samfile.close()
	# process_list = list()
	# for i in samfile.get_index_statistics():
	# 	process_list.append([i[0], i[3]])
	# 	# #chr #read
	# process_list = sorted(process_list, key = lambda x:-x[1])
	# Chr_name = process_list[0][0]
	# print("[INFO]: %d reads on %s."%(process_list[0][1], process_list[0][0]))

	# for read in samfile.fetch(Chr_name):
	# 	parse_read(read, Chr_name, 50, 20, 7, 2000)
	# samfile.close()
	# print("[INFO]: Parse reads used %0.2f seconds."%(time.time() - starttime))

	for res in result:
		temp = res.get()[0]

		for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]:
			for chr in temp[sv_type]:
				if chr not in MainCandidate[sv_type]:
					MainCandidate[sv_type][chr] = list()
				MainCandidate[sv_type][chr].extend(temp[sv_type][chr])
		for sv_type in ["l_DUP", "DUP_r", "l_INV", "INV_r"]:
			for chr in temp[sv_type]:
				if chr not in MainCandidate[sv_type]:
					MainCandidate[sv_type][chr] = dict()
				for hash_1 in temp[sv_type][chr]:
					if hash_1 not in MainCandidate[sv_type][chr]:
						MainCandidate[sv_type][chr][hash_1] = dict()
					for hash_2 in temp[sv_type][chr][hash_1]:
						if hash_2 not in MainCandidate[sv_type][chr][hash_1]:
							MainCandidate[sv_type][chr][hash_1][hash_2] = list()
						MainCandidate[sv_type][chr][hash_1][hash_2].extend(temp[sv_type][chr][hash_1][hash_2])

	show_temp_result(args.min_support, args.min_length, 10, args.max_distance, args.output, MainCandidate)
	# print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Detection of all types of Structural Variants.
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="cuteSV", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input", metavar="[BAM]", type=str, help="Sorted .bam file form NGMLR.")
	parser.add_argument('output', type=str, help = "the prefix of novel sequence insertion pridections")

	parser.add_argument('-s', '--min_support', help = "Minimum number of reads that support a SV to be reported.[%(default)s]", default = 5, type = int)
	parser.add_argument('-l', '--min_length', help = "Minimum length of SV to be reported.[%(default)s]", default = 50, type = int)
	parser.add_argument('-p', '--max_split_parts', help = "Maximum number of split segments a read may be aligned before it is ignored.[%(default)s]", default = 7, type = int)
	# # parser.add_argument('-hom', '--homozygous', help = "The mininum score of a genotyping reported as a homozygous.[%(default)s]", default = 0.8, type = float)
	# # parser.add_argument('-het','--heterozygous', help = "The mininum score of a genotyping reported as a heterozygous.[%(default)s]", default = 0.3, type = float)
	parser.add_argument('-q', '--min_mapq', help = "Minimum mapping quality value of alignment to be taken into account.[%(default)s]", default = 20, type = int)
	parser.add_argument('-d', '--max_distance', help = "Maximum distance to group SV together..[%(default)s]", default = 1000, type = int)
	parser.add_argument('-r', '--min_seq_size', help = "Ignores reads that only report alignments with not longer then bp.[%(default)s]", default = 2000, type = int)
	parser.add_argument('-t', '--threads', help = "Number of threads to use.[%(default)s]", default = 16, type = int)
	args = parser.parse_args(argv)

	return args

def run(argv):
	# import time
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	# main_ctrl(sys.argv[1], sys.argv[2])
	run(sys.argv[1:])
