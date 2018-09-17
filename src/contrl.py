import argparse
import sys
import logging
import pysam
import cigar
from multiprocessing import Pool
from __utils__ import *
import time




def main_ctrl(sam_path, out_path):
	starttime = time.time()
	samfile = pysam.AlignmentFile(sam_path)
	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:-x[1])
	Chr_name = process_list[0][0]
	print("[INFO]: %d reads on %s."%(process_list[0][1], process_list[0][0]))
	# Chr_length = samfile.get_reference_length(Chr_name)
	# print Chr_name, Chr_length
	# logging.info("Resolving the chromsome %s."%(Chr_name))
	for read in samfile.fetch(Chr_name):
		parse_read(read, Chr_name, 50, 20, 7, 2000)
	samfile.close()
	print("[INFO]: Parse reads used %0.2f seconds."%(time.time() - starttime))

	show_temp_result(5, 50, 10, 1000, out_path)
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	main_ctrl(sys.argv[1], sys.argv[2])
