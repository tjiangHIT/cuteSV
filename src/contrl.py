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
	process_list = sorted(process_list, key = lambda x:x[1])

	for i in process_list:
		Chr_name = i[0]
		local_starttime = time.time()
		print("[INFO]: %d reads on %s."%(i[1], i[0]))
	# Chr_length = samfile.get_reference_length(Chr_name)
	# print Chr_name, Chr_length
	# logging.info("Resolving the chromsome %s."%(Chr_name))
		for read in samfile.fetch(Chr_name):
			parse_read(read, Chr_name, 50, 20, 4, 2000)
		print("[INFO]: Parse reads used %0.2f seconds."%(time.time() - local_starttime))
	samfile.close()

	show_temp_result(5, 50, 10, 1000, out_path)
	print("[INFO]: Finished in %0.2f seconds."%(time.time() - starttime))

if __name__ == '__main__':
	main_ctrl(sys.argv[1], sys.argv[2])
