#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  extract.py
 * @Package: argparse, pysam, sys, Bio, os, logging
 * @Description: Parse the ME signals from alignments
 * @author: tjiang
 * @date: Apr 24 2018
 * @version V1.0.0   
'''

VERSION = '1.0.0'

USAGE="""\
	Detection of Structural Variants using long reads.

	cuteSV V%s
"""%(VERSION)

import pysam
import cigar
from multiprocessing import Pool
from CommandRunner import *
from resolution_type import * 
import os
import argparse
import logging
import sys
import time
import gc

dic_starnd = {1:'+', \
		2: '-'}
signal = {1 << 2: 0, \
		1 >> 1: 1, \
		1 << 4: 2, \
		1 << 11: 3, \
		1 << 4 | 1 << 11: 4}
'''
	1 >> 1 means normal_foward read
	1 << 2 means unmapped read
	1 << 4 means reverse_complement read
	1 << 11 means supplementary alignment read
	1 << 4 | 1 << 11 means supplementary alignment with reverse_complement read
'''
def detect_flag(Flag):
	back_sig = signal[Flag] if Flag in signal else 0
	return back_sig

def store_info(pos_list, svtype, candidate):
	for ele in pos_list:
		if ele[0] not in candidate[svtype]:
			candidate[svtype][ele[0]] = dict()
		hash_1 = int(ele[1] /10000)
		mod = ele[1] % 10000
		hash_2 = int(mod / 50)
		if hash_1 not in candidate[svtype][ele[0]]:
			candidate[svtype][ele[0]][hash_1] = dict()
			candidate[svtype][ele[0]][hash_1][hash_2] = list()
			candidate[svtype][ele[0]][hash_1][hash_2].append(ele[1:])
		else:
			if hash_2 not in candidate[svtype][ele[0]][hash_1]:
				candidate[svtype][ele[0]][hash_1][hash_2] = list()
			candidate[svtype][ele[0]][hash_1][hash_2].append(ele[1:])

def analysis_split_read(split_read, SV_size, RLength, read_name, candidate):
	'''
	read_start	read_end	ref_start	ref_end	chr	strand
	#0			#1			#2			#3		#4	#5
	'''
	SP_list = sorted(split_read, key = lambda x:x[0])
	# print read_name
	# for i in SP_list:
	# 	print i
	DUP_flag = [0]*len(SP_list)
	for a in range(len(SP_list[:-1])):
		ele_1 = SP_list[a]
		ele_2 = SP_list[a+1]
		if ele_1[4] == ele_2[4]:
			# INV
			if ele_1[5] != ele_2[5]:
				# exist a bug 
				#print "this", ele_1[3], ele_2[3]
				if ele_2[3] - ele_1[3] >= SV_size:
					if ele_2[0] >= ele_1[1]:
						# print "that", ele_1[3], ele_2[3]
						# print "that", ele_1[2], ele_2[2]
						if ele_1[4] not in candidate["INV"]:
							candidate["INV"][ele_1[4]] = list()
						candidate["INV"][ele_1[4]].append([ele_1[3], ele_2[3], read_name])
						candidate["INV"][ele_1[4]].append([ele_1[2], ele_2[2], read_name])
				#	# candidate["INV"][a[4]].append([a[3] , call_inv[call_inv.index(a)+1][3], read_name])
				#	print "this\tINV", ele_1[4], ele_1[3], ele_2[3], read_name
				#	candidate["INV"][ele_1[4]].append([ele_1[3], ele_2[3], read_name])
				if ele_1[3] - ele_2[3] >= SV_size:
					if ele_2[0] >= ele_1[1]:
						if ele_1[4] not in candidate["INV"]:
							candidate["INV"][ele_1[4]] = list()
						candidate["INV"][ele_1[4]].append([ele_2[3], ele_1[3], read_name])
						candidate["INV"][ele_1[4]].append([ele_2[2], ele_1[2], read_name])

			else:
				# dup & ins & del 
				if ele_1[5] == '-':
					ele_1 = [RLength-SP_list[a+1][1], RLength-SP_list[a+1][0]]+SP_list[a+1][2:]
					ele_2 = [RLength-SP_list[a][1],RLength-SP_list[a][0]]+SP_list[a][2:]
					# print ele_1
					# print ele_2

				if ele_1[3] - ele_2[2] >= SV_size:
					# DUP_flag[a] = 1
					# DUP_flag[a+1] = 1
					'''
					if a == 0:
						store_info([[ele_1[4], ele_1[3], read_name]], "DUP_r", candidate)

					if a+1 == len(SP_list)-1:
						store_info([[ele_2[4], ele_2[2], read_name]], "l_DUP", candidate)
					else:
						# store_signal([[ele_2[4], ele_2[2], ele_2[3]]], "DUP")
						if ele_2[4] not in candidate["DUP"]:
							candidate["DUP"][ele_2[4]] = list()
						candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])
					'''
					if ele_2[4] not in candidate["DUP"]:
						candidate["DUP"][ele_2[4]] = list()
					candidate["DUP"][ele_2[4]].append([ele_2[2], ele_1[3], read_name])

				if ele_1[3] <= ele_2[2]:
					if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
						# store_signal([[ele_2[4], (ele_2[2]+ele_1[3])/2, ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1]]], "INS")
						if ele_2[4] not in candidate["INS"]:
							candidate["INS"][ele_2[4]] = list()
						candidate["INS"][ele_2[4]].append([(ele_2[2]+ele_1[3])/2, ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1], read_name])
					if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
						# store_signal([[ele_2[4], ele_1[3], ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3]]], "DEL")
						# print ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]
						if ele_2[4] not in candidate["DEL"]:
							candidate["DEL"][ele_2[4]] = list()
						candidate["DEL"][ele_2[4]].append([ele_1[3], ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3], read_name])
		else:
			# tra
			if ele_1[4] < ele_2[4]:
				# store_signal([[ele_1[4], ele_1[3], ele_2[4], ele_2[2]]], "TRA")
				if ele_1[4] not in candidate["TRA"]:
					candidate["TRA"][ele_1[4]] = list()
				candidate["TRA"][ele_1[4]].append([ele_1[3], ele_2[4], ele_2[2], read_name])
			else:
				# store_signal([[ele_2[4], ele_2[2], ele_1[4], ele_1[3]]], "TRA")
				if ele_2[4] not in candidate["TRA"]:
					candidate["TRA"][ele_2[4]] = list()
				candidate["TRA"][ele_2[4]].append([ele_2[2], ele_1[4], ele_1[3], read_name])

def acquire_clip_pos(deal_cigar):
	seq = list(cigar.Cigar(deal_cigar).items())
	if seq[0][1] == 'S':
		first_pos = seq[0][0]
	else:
		first_pos = 0
	if seq[-1][1] == 'S':
		last_pos = seq[-1][0]
	else:
		last_pos = 0

	bias = 0
	for i in seq:
		if i[1] == 'M' or i[1] == 'D' or i[1] == '=' or i[1] == 'X':
			bias += i[0]
	return [first_pos, last_pos, bias]

def organize_split_signal(chr, primary_info, Supplementary_info, total_L, low_bandary, min_mapq, max_split_parts, read_name, candidate):
	split_read = list()
	if len(primary_info) > 0:
		split_read.append(primary_info)
	for i in Supplementary_info:
		seq = i.split(',')
		local_chr = seq[0]
		local_start = int(seq[1])
		local_cigar = seq[3]
		local_strand = seq[2]
		local_mapq = int(seq[4])
		if local_mapq >= min_mapq:
			local_set = acquire_clip_pos(local_cigar)
			# split_read.append([local_set[0], total_L-local_set[1], local_start, local_start+local_set[2], local_chr, local_strand])
			# print read_name, local_set, local_strand
			if local_strand == '+':
			 	split_read.append([local_set[0], total_L-local_set[1], local_start, local_start+local_set[2], local_chr, local_strand])
			else:
				# if read_name == "bfe83948_158430_4316":
				# 	print "this", seq
				# 	print "this", [local_set[1], total_L-local_set[0], local_start, local_start+local_set[2], local_chr, local_strand]
				# 	# split_read.append([local_set[1], total_L-local_set[0], local_start, local_start+local_set[2], local_chr, local_strand])
				# 	continue
				try:
					split_read.append([local_set[1], total_L-local_set[0], local_start, local_start+local_set[2], local_chr, local_strand])
				except:
					pass
	if len(split_read) <= max_split_parts:
		analysis_split_read(split_read, low_bandary, total_L, read_name, candidate)

# def parse_read(read, Chr_name, low_bandary):
def parse_read(read, Chr_name, low_bandary, min_mapq, max_split_parts, min_seq_size, candidate):
	if read.query_length < min_seq_size:
		return 0

	process_signal = detect_flag(read.flag)
	if read.mapq >= min_mapq:
		pos_start = read.reference_start
		pos_end = read.reference_end
		shift_del = 0
		shift_ins = 0
		softclip_left = 0
		softclip_right = 0
		for element in read.cigar:
			if element[0] in [0, 7 ,8]:
				shift_del += element[1]
			if element[0] == 2 and element[1] < low_bandary:
				shift_del += element[1]
			if element[0] == 2 and element[1] >= low_bandary:
				# DEL_ME_pos.append([pos_start+shift, element[1]])
				# store_signal([Chr_name, pos_start+shift, element[1]], "DEL")
				if Chr_name not in candidate["DEL"]:
					candidate["DEL"][Chr_name] = list()
				candidate["DEL"][Chr_name].append([pos_start+shift_del, element[1], read.query_name])
				shift_del += element[1]

			if element[0] in [0, 2, 7, 8]:
				shift_ins += element[1]
			if element[0] == 1 and element[1] >= low_bandary:
				shift_ins += 1
				# MEI_contig = read.query_sequence[_shift_read_ - element[1]:_shift_read_]
				if Chr_name not in candidate["INS"]:
					candidate["INS"][Chr_name] = list()
				candidate["INS"][Chr_name].append([pos_start+shift_ins, element[1], read.query_name])

		if read.cigar[0][0] == 4:
			softclip_left = read.cigar[0][1]
		if read.cigar[-1][0] == 4:
			softclip_right = read.cigar[-1][1]

	# print read.query_name, process_signal
	if process_signal == 1 or process_signal == 2:
		Tags = read.get_tags()
		if read.mapq >= min_mapq:
			if process_signal == 1:
				primary_info = [softclip_left, read.query_length-softclip_right, pos_start, pos_end, Chr_name, dic_starnd[process_signal]]
			else:
				primary_info = [softclip_right, read.query_length-softclip_left, pos_start, pos_end, Chr_name, dic_starnd[process_signal]]
		else:
			primary_info = []

		for i in Tags:
			if i[0] == 'SA':
				Supplementary_info = i[1].split(';')[:-1]
				# print read.query_name
				organize_split_signal(Chr_name, primary_info, Supplementary_info, read.query_length, low_bandary, min_mapq, max_split_parts, read.query_name, candidate)

def single_pipe(sam_path, min_length, min_mapq, max_split_parts, min_seq_size, temp_dir, task):
	# args.input, args.min_length, args.min_mapq, args.max_split_parts, args.min_seq_size, temporary_dir, i
	candidate = dict()
	candidate["DEL"] = dict()
	candidate["INS"] = dict()
	candidate["INV"] = dict()
	candidate["DUP"] = dict()
	candidate["TRA"] = dict()
	candidate["l_DUP"] = dict()
	candidate["DUP_r"] = dict()
	Chr_name = task[0]
	samfile = pysam.AlignmentFile(sam_path)
	Chr_length = samfile.get_reference_length(Chr_name)
	#logging.info("Resolving the chromsome %s:%d-%d."%(Chr_name, task[1], task[2]))

	for read in samfile.fetch(Chr_name, task[1], task[2]):
		parse_read(read, Chr_name, min_length, min_mapq, max_split_parts, min_seq_size, candidate)

	logging.info("Finished %s:%d-%d."%(Chr_name, task[1], task[2]))
	samfile.close()

	# os.mkdir("%ssignatures"%temp_dir)
	output = "%ssignatures/_%s_%d_%d.bed"%(temp_dir, Chr_name, task[1], task[2])
	#print output 
	#print candidate["INV"], output
	file = open(output, 'w')
	for sv_type in ["DEL", "INS", "INV", "DUP", 'TRA']:
		try:
			for chr in candidate[sv_type]:
				# print candidate[sv_type][chr]
				for ele in candidate[sv_type][chr]:
					if len(ele) == 3:
						file.write("%s\t%s\t%d\t%d\t%s\n"%(sv_type, chr, ele[0], ele[1], ele[2]))
					elif len(ele) == 4:
						file.write("%s\t%s\t%d\t%s\t%d\t%s\n"%(sv_type, chr, ele[0], ele[1], ele[2], ele[3]))
					#print ele
		except:
			pass

	# for sv_type in ["l_DUP", "DUP_r", "l_INV", "INV_r"]:
	for sv_type in ["l_DUP", "DUP_r"]:
		try:
			for chr in candidate[sv_type]:
				# print candidate[sv_type][chr]
				for key_1 in candidate[sv_type][chr]:
					for key_2 in candidate[sv_type][chr][key_1]:
						for ele in candidate[sv_type][chr][key_1][key_2]:
							file.write("%s\t%s\t%d\t%d\t%d\t%s\n"%(sv_type, chr, key_1, key_2, ele[0], ele[1]))
							# for i in ele:
								# file.write("%s\t%s\t%d\t%d\t%d\t%s\n"%(sv_type, chr, key_1, key_2, i[0], i[1]))
		except:
			pass
	file.close()	
	# return candidate

def multi_run_wrapper(args):
	return single_pipe(*args)

def main_ctrl(args):
	samfile = pysam.AlignmentFile(args.input)
	contig_num = len(samfile.get_index_statistics())
	logging.info("The total number of chromsomes: %d"%(contig_num))

	Task_list = list()
	# process_list = list()
	chr_name_list = list()

	ref_ = samfile.get_index_statistics()
	for i in ref_:
		# process_list.append([i[0], i[3]])
		chr_name_list.append(i[0])
		# print i[0], samfile.get_reference_length(i[0])
		local_ref_len = samfile.get_reference_length(i[0])
		if local_ref_len < args.batches:
			Task_list.append([i[0], 0, local_ref_len])
		else:
			pos = 0
			task_round = int(local_ref_len/args.batches)
			for j in range(task_round):
				Task_list.append([i[0], pos, pos+args.batches])
				pos += args.batches
			if pos < local_ref_len:
				Task_list.append([i[0], pos, local_ref_len])

	analysis_pools = Pool(processes=int(args.threads))

	# # # result = list()
	if args.temp_dir[-1] == '/':
		temporary_dir = args.temp_dir
	else:
		temporary_dir = args.temp_dir+'/'
	# '''
	os.mkdir("%ssignatures"%temporary_dir)
	for i in Task_list:
		para = [(args.input, args.min_length, args.min_mapq, args.max_split_parts, args.min_seq_size, temporary_dir, i)]
		analysis_pools.map_async(multi_run_wrapper, para)
	analysis_pools.close()
	analysis_pools.join()
	# '''
	# '''
	logging.info("Rebuilding ssignatures of structural variants.")
	analysis_pools = Pool(processes=int(args.threads))
	cmd_del = ("cat %ssignatures/*.bed | grep DEL | sort -u | sort -k 2,2 -k 3,3n > %sDEL.sigs"%(temporary_dir, temporary_dir))
	cmd_ins = ("cat %ssignatures/*.bed | grep INS | sort -u | sort -k 2,2 -k 3,3n > %sINS.sigs"%(temporary_dir, temporary_dir))
	cmd_inv = ("cat %ssignatures/*.bed | grep INV | sort -u | sort -k 2,2 -k 3,3n > %sINV.sigs"%(temporary_dir, temporary_dir))
	cmd_tra = ("cat %ssignatures/*.bed | grep TRA | sort -u | sort -k 2,2 -k 4,4 -k 3,3n > %sTRA.sigs"%(temporary_dir, temporary_dir))
	cmd_dup = ("cat %ssignatures/*.bed | grep DUP | sort -u | sort -k 1,1r -k 2,2 -k 3,4n > %sDUP.sigs"%(temporary_dir, temporary_dir))
	for i in [cmd_ins, cmd_del, cmd_dup, cmd_tra, cmd_inv]:
		analysis_pools.map_async(exe, (i,))
	analysis_pools.close()
	analysis_pools.join()
	# '''

	chr_name_list.sort()
	process_tra = list()
	for i in range(len(chr_name_list)-1):
		# print i
		for j in chr_name_list[i+1:]:
			process_tra.append([chr_name_list[i], j])
			# print chr_name_list[i], j

	logging.info("Clustering structural variants.")
	analysis_pools = Pool(processes=int(args.threads))
	# print "%s%s.sigs"%(temporary_dir, "DEL")
	# local_res = resolution_INDEL("%s%s.sigs"%(temporary_dir, "DEL"), "1", "DEL", 10, 0.5, 20, 50)
	# print len(local_res)
	result = list()
	for chr in chr_name_list:
		for svtype in ["DEL", "INS", "INV"]:
			if svtype == "INV":
				# para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.5, 20, args.min_length)]
				# para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.5, args.max_cluster_bias, args.min_length, args.max_distance)]
				pass
			if svtype == 'DEL':
				para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.3, 200, 0.7, 5)]
				result.append(analysis_pools.map_async(run_del, para))
			if svtype == 'INS':
				para = [("%s%s.sigs"%(temporary_dir, svtype), chr, svtype, args.min_support, 0.2, 100, 0.6, 5)]
				result.append(analysis_pools.map_async(run_ins, para))
			# resolution_INDEL(sys.argv[1], "1", "INV", 10, 0.5, 20, 50)
			# result.append(analysis_pools.map_async(run_indel_inv, para))
			# print chr, svtype
		# DUP
		# resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size=50)
		# resolution_DUP(sys.argv[1], "1", 10, 50, 50)
		# para = [("%s%s.sigs"%(temporary_dir, "DUP"), chr, args.min_support, args.max_cluster_bias, args.min_length, args.max_distance)]
		# result.append(analysis_pools.map_async(run_dup, para))

	# for i in process_tra:
	# 	# TRA
	# 	# resolution_TRA(path, chr_1, chr_2, read_count, overlap_size, max_cluster_bias)
	# 	para = [("%s%s.sigs"%(temporary_dir, "TRA"), i[0], i[1], args.min_support, 0.6, args.max_cluster_bias)]
	# 	result.append(analysis_pools.map_async(run_tra, para))

	analysis_pools.close()
	analysis_pools.join()

	semi_result = list()
	for res in result:
		try:
			semi_result += res.get()[0]
		except:
			pass

	logging.info("Writing into disk.")
	file = open(args.output, 'w')
	for i in semi_result:
		file.write(i)
	file.close()

	# logging.info("Cleaning temporary files.")
	# cmd_remove_tempfile = ("rm -r %ssignatures %s*.sigs"%(temporary_dir, temporary_dir))
	# exe(cmd_remove_tempfile)
	
	samfile.close()

def cal_GT(a, b):
	if b == 0:
		return "1/1"
	if a*1.0/b < 0.3:
		return "0/0"
	elif a*1.0/b >= 0.3 and a*1.0/b < 0.8:
		return "0/1"
	elif a*1.0/b >= 0.8 and a*1.0/b < 1.0:
		return "1/1"
	else:
		return "1/1"

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

def run(argv):
	args = parseArgs(argv)
	setupLogging(False)
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

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

if __name__ == '__main__':
	run(sys.argv[1:])
