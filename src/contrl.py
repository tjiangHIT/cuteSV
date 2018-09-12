import argparse
import sys
import logging
import pysam
import cigar
from multiprocessing import Pool
from __utils__ import *




def main_ctrl(sam_path):
	samfile = pysam.AlignmentFile(sam_path)
	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:-x[1])
	Chr_name = process_list[0][0]
	# Chr_length = samfile.get_reference_length(Chr_name)
	# print Chr_name, Chr_length
	# logging.info("Resolving the chromsome %s."%(Chr_name))
	for read in samfile.fetch(Chr_name):
		parse_read(read, Chr_name, 50, 0)
	samfile.close()

	show_temp_result()

if __name__ == '__main__':
	main_ctrl(sys.argv[1])
