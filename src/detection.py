#!/usr/bin/env python

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  extract.py
 * @Package: argparse, pysam, sys, Bio, os, logging
 * @Description: Parse the ME signals from alignments
 * @author: tjiang
 * @date: Apr 24 2018
 * @version V1.0     
'''

import pysam
import cigar
from multiprocessing import Pool
import os
import argparse
import logging
import sys, time
import gc

dic_starnd = {1:'+', 2: '-'}
Normal_foward = 1 >> 1
Abnormal = 1 << 2
Reverse_complement = 1 << 4
Supplementary_map = 1 << 11
signal = {Abnormal: 0, Normal_foward: 1, Reverse_complement: 2, Supplementary_map:3, Reverse_complement | Supplementary_map:4}

candidate = dict()
semi_result = dict()
dup_result = dict()
candidate["DEL"] = dict()
candidate["INS"] = dict()
candidate["INV"] = dict()
candidate["DUP"] = dict()
candidate["TRA"] = dict()
candidate["l_DUP"] = dict()
candidate["DUP_r"] = dict()
candidate["l_INV"] = dict()
candidate["INV_r"] = dict()

def detect_flag(Flag):
	if Flag in signal:
		return signal[Flag]
	else:
		return 0

def store_info(pos_list, svtype):
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

def analysis_split_read(split_read, SV_size, RLength, read_name):
	SP_list = sorted(split_read, key = lambda x:x[0])
	# for i in SP_list:
	# 	print i
	DUP_flag = [0]*len(SP_list)
	for a in xrange(len(SP_list[:-1])):
		ele_1 = SP_list[a]
		ele_2 = SP_list[a+1]
		if ele_1[4] == ele_2[4]:
			# dup & ins & del 
			if ele_1[3] - ele_2[2] >= SV_size:
				# DUP_flag[a] = 1
				# DUP_flag[a+1] = 1
				if a == 0:
					store_info([[ele_1[4], ele_1[3], read_name]], "DUP_r")

				if a+1 == len(SP_list)-1:
					store_info([[ele_2[4], ele_2[2], read_name]], "l_DUP")
				else:
					# store_signal([[ele_2[4], ele_2[2], ele_2[3]]], "DUP")
					if ele_2[4] not in candidate["DUP"]:
						candidate["DUP"][ele_2[4]] = list()
					candidate["DUP"][ele_2[4]].append([ele_2[2], ele_2[3], read_name])

			if ele_1[3] <= ele_2[2]:
				if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
					# store_signal([[ele_2[4], (ele_2[2]+ele_1[3])/2, ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1]]], "INS")
					if ele_2[4] not in candidate["INS"]:
						candidate["INS"][ele_2[4]] = list()
					candidate["INS"][ele_2[4]].append([(ele_2[2]+ele_1[3])/2, ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1], read_name])
				if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
					# store_signal([[ele_2[4], ele_1[3], ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3]]], "DEL")
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

	call_inv = sorted(SP_list, key = lambda x:x[2])
	if len(call_inv) >= 3:
		for a in call_inv[:-2]:
			if a[5] != call_inv[call_inv.index(a)+1][5] and a[5] == call_inv[call_inv.index(a)+2][5] and a[4] == call_inv[call_inv.index(a)+1][4] and a[4] == call_inv[call_inv.index(a)+2][4]:
				if call_inv[call_inv.index(a)+1][3] - call_inv[call_inv.index(a)+1][2] >= SV_size:
					# store_signal([[a[4], a[3] , call_inv[call_inv.index(a)+1][3]]], "INV")
					if a[4] not in candidate["INV"]:
						candidate["INV"][a[4]] = list()
					candidate["INV"][a[4]].append([a[3] , call_inv[call_inv.index(a)+1][3], read_name])
	
	if len(call_inv) == 2:
		if call_inv[0][5] != call_inv[1][5] and call_inv[0][4] == call_inv[1][4]:
			ls_1 = call_inv[0][3] - call_inv[0][2]
			ls_2 = call_inv[1][3] - call_inv[1][2]
			if ls_1 > ls_2:
				if call_inv[1][2] > call_inv[0][3] and ls_2 >= SV_size:
					store_info([[call_inv[0][4], call_inv[1][3], read_name]], "INV_r")
					# pass
			else:
				if call_inv[1][2] > call_inv[0][3] and ls_1 >= SV_size:
					store_info([[call_inv[0][4], call_inv[0][2], read_name]], "l_INV")
					# pass

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
		if i[1] == 'M' or i[1] == 'D':
			bias += i[0]
	return [first_pos, last_pos, bias]

def organize_split_signal(chr, primary_info, Supplementary_info, total_L, low_bandary, min_mapq, max_split_parts, read_name):
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
			split_read.append([local_set[0], total_L-local_set[1], local_start, local_start+local_set[2], local_chr, local_strand])
	if len(split_read) <= max_split_parts:
		analysis_split_read(split_read, low_bandary, total_L, read_name)

# def parse_read(read, Chr_name, low_bandary):
def parse_read(read, Chr_name, low_bandary, min_mapq, max_split_parts, min_seq_size):
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
			if element[0] == 0:
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

			if element[0] == 0 or element[0] == 2:
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


	if process_signal == 1 or process_signal == 2:
		Tags = read.get_tags()
		if read.mapq >= min_mapq:
			primary_info = [softclip_left, read.query_length-softclip_right, pos_start, pos_end, Chr_name, dic_starnd[process_signal]]
		else:
			primary_info = []

		for i in Tags:
			if i[0] == 'SA':
				Supplementary_info = i[1].split(';')[:-1]
				# print read.query_name
				organize_split_signal(Chr_name, primary_info, Supplementary_info, read.query_length, low_bandary, min_mapq, max_split_parts, read.query_name)

def acquire_length(len_list, threshold):
	average_len = 0
	flag = 0
	temp = list()
	temp.append(len_list[0])
	for ele in len_list[1:]:
		if temp[-1] + threshold < ele:
			if flag < len(temp):
				flag = len(temp)
				average_len = int(sum(temp)/flag)
			temp = list()
			temp.append(ele)
		else:
			temp.append(ele)
	if flag < len(temp):
		flag = len(temp)
		average_len = int(sum(temp)/flag)
	return average_len

# @jit
def merge_pos_indel(pos_list, chr, evidence_read, SV_size, svtype):
	# print pos_list
	# result = list()
	if len(pos_list) >= evidence_read:
		start = list()
		size = list()
		diff = list()
		tag = dict()
		# max_L = -1
		# min_L = sys.maxint
		for ele in pos_list:
			start.append(ele[0])
			size.append(ele[1])
			if ele[2] not in tag:
				tag[ele[2]] = 0

		if len(tag) < evidence_read:
			return 0

		for k in xrange(len(pos_list[1:])):
			diff.append(pos_list[k+1][1]-pos_list[k][1])
			# if ele[1] > max_L:
			# 	max_L = ele[1]
			# if ele[1] < min_L:
			# 	min_L = ele[1]
		size.sort()
		breakpoint = sum(start)/len(start)
		# SV_len = sum(size[int(0.25*len(size)):int(0.75*len(size))])/len(size[int(0.25*len(size)):int(0.75*len(size))])
		SV_len = acquire_length(size, max(int(sum(diff)/len(diff)), 50))
		# concider the stability of data
		# need to repair the buf!!!!!!!!

		# if SV_len >= SV_size and size[int(0.75*len(size))] - size[int(0.25*len(size))] < int(0.5*SV_len):
		if SV_len >= SV_size:
			# result.append([chr, breakpoint, SV_len, len(pos_list)])
			if chr not in semi_result:
				semi_result[chr] = list()
			# semi_result[chr].append([breakpoint, SV_len, len(pos_list), svtype])
			semi_result[chr].append("%d\t%d\t%d\t%s\n"%(breakpoint, SV_len, len(tag), svtype))

	# return result
# @jit
def intergrate_indel(chr, evidence_read, SV_size, low_bandary, svtype, max_distance, MainCandidate):
	# _cluster_ = list()
	# print MainCandidate[svtype][chr]
	temp = list()
	temp.append(MainCandidate[svtype][chr][0])
	for pos in MainCandidate[svtype][chr][1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			# print temp[-1][0] - temp[0][0]
			if temp[-1][0] - temp[0][0] <= max_distance:
				merge_pos_indel(temp, chr, evidence_read, SV_size, svtype)
			# if len(result) != 0:
			# 	_cluster_.append(result[0])
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if temp[-1][0] - temp[0][0] <= max_distance:
		merge_pos_indel(temp, chr, evidence_read, SV_size, svtype)

def acquire_locus(down, up, keytype, chr, MainCandidate):
	re = list()
	re_id = list()
	if chr not in MainCandidate[keytype]:
		return re, re_id
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in MainCandidate[keytype][chr]:
			return re, re_id
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					re.append(ele[0])
					re_id.append(ele[1])
	else:
		key_1 = int(down/10000)
		if key_1 in MainCandidate[keytype][chr]:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in MainCandidate[keytype][chr][key_1]:
					continue
				for ele in MainCandidate[keytype][chr][key_1][key_2]:
					if ele[0] >= down and ele[0] <= up:
						re.append(ele[0])
						re_id.append(ele[1])
		key_1 += 1
		if key_1 not in MainCandidate[keytype][chr]:
			return re, re_id
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					re.append(ele[0])
					re_id.append(ele[1])
	return re, re_id
# @jit
def merge_pos_dup(pos_list, chr, evidence_read, SV_size, ll, lr, svtype, MainCandidate):
	start = list()
	end = list()
	tag = dict()
	diff = list()
	# up_B = int(len(pos_list)*0.75)
	# low_B = int(len(pos_list)*0.25)
	for ele in pos_list:
		start.append(ele[0])
	# for ele in pos_list[low_B:up_B+1]:
		end.append(ele[1])
		# tag.append(ele[2])
		if ele[2] not in tag:
			# tag.append(ele[2])
			tag[ele[2]] = 0

	re_l, rl_id = acquire_locus(min(start)-50, max(start)+50, ll, chr, MainCandidate)
	re_r, rr_id = acquire_locus(min(end)-50, max(end)+50, lr, chr, MainCandidate)

	for flag in rl_id:
		if flag not in tag:
			# tag.append(flag)
			tag[flag] = 0
	for flag in rr_id:
		if flag not in tag:
			# tag.append(flag)
			tag[flag] = 0

	# print pos_list
	# print re_l
	# print re_r

	breakpoint = sum(start+re_l)/len(start+re_l)
	# size = sum(end+re_r)/len(end+re_r) - breakpoint
	end = end + re_r
	end.sort()
	for k in xrange(len(end)-1):
		diff.append(end[k+1] - end[k])
	# print end
	# print diff
	if len(diff) > 0:
		size = acquire_length(end, max(int(sum(diff)/len(diff)), 50)) - breakpoint
	else:
		size = acquire_length(end, 50) - breakpoint
	# print breakpoint, size, re_r, re_l, start, end

	# result = list()
	# if len(pos_list+re_l+re_r) >= evidence_read and len(pos_list) >= 5:
	if len(tag) >= evidence_read:
		# need to add some semi-signals
		# if size >= SV_size and max(end+re_r) - min(end+re_r) < int(0.5*size):
		if size >= SV_size:
			# result.append([chr, breakpoint, size, len(pos_list+re_l+re_r)])

			# if chr not in semi_result:
			# 	semi_result[chr] = list()
			# # semi_result[chr].append([breakpoint, size, len(pos_list+re_l+re_r), svtype])
			# semi_result[chr].append("%d\t%d\t%d\t%s\n"%(breakpoint, size, len(tag), svtype))
			if chr not in dup_result:
				dup_result[chr] = list()
			dup_result[chr].append([breakpoint, size, len(tag)])

def polish_dup(chr):
	if chr in dup_result:
		dup_result[chr] = sorted(dup_result[chr], key = lambda x:x[:])
		# temp = list()
		# temp.append(dup_result[chr][0])
		temp = dup_result[chr][0]
		if chr not in semi_result:
			semi_result[chr] = list()
		for i in dup_result[chr][1:]:
			if temp[0] <= i[0] and i[0] <= temp[0]+temp[1]:
				if temp[2] < i[2]:
					temp = i
			else:
				semi_result[chr].append("%d\t%d\t%d\t%s\n"%(temp[0], temp[1], temp[2], "DUP"))
				temp = i
		semi_result[chr].append("%d\t%d\t%d\t%s\n"%(temp[0], temp[1], temp[2], "DUP"))


# @jit
def intergrate_dup(chr, evidence_read, SV_size, low_bandary, ll, lr, svtype, max_distance, MainCandidate):
	temp = list()
	temp.append(MainCandidate[svtype][chr][0])
	for pos in MainCandidate[svtype][chr][1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			if temp[-1][0] - temp[0][0] <= max_distance:
				merge_pos_dup(temp, chr, evidence_read, SV_size, ll, lr, svtype, MainCandidate)
				# print temp
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if temp[-1][0] - temp[0][0] <= max_distance:
		merge_pos_dup(temp, chr, evidence_read, SV_size, ll, lr, svtype, MainCandidate)

def merge_pos_tra(pos_list, chr, evidence_read, SV_size, svtype, low_bandary):
	if len(pos_list) >= evidence_read:
		start = list()
		tra_dic = dict()
		tag = dict()
		for ele in pos_list:
			if ele[3] not in tag:
				tag[ele[3]] = 0

			start.append(ele[0])
			# size.append(ele[1])
			if ele[1] not in tra_dic:
				tra_dic[ele[1]] = list()
			tra_dic[ele[1]].append(ele[2])

		if len(tag) < evidence_read:
			return 0

		max_chr_c = 0
		max_chr = ""
		for key in tra_dic:
			if len(tra_dic[key]) > max_chr_c:
				max_chr = key
				max_chr_c = len(tra_dic[key])

		if max_chr_c >= evidence_read:
			breakpoint_1 = int(sum(start)/len(start))
			breakpoint_2 = int(sum(tra_dic[max_chr])/max_chr_c)

			if chr not in semi_result:
				semi_result[chr] = list()
			# semi_result[chr].append([breakpoint_1, max_chr, breakpoint_2, svtype])
			semi_result[chr].append("%d\t%s\t%s\t%d\t%s\n"%(breakpoint_1, max_chr, breakpoint_2, len(tag), svtype))
# @jit
def intergrate_tra(chr, evidence_read, SV_size, low_bandary, svtype, max_distance, MainCandidate):
	temp = list()
	temp.append(MainCandidate[svtype][chr][0])
	for pos in MainCandidate[svtype][chr][1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			if temp[-1][0] - temp[0][0] <= max_distance:
				merge_pos_tra(temp, chr, evidence_read, SV_size, svtype, low_bandary)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if temp[-1][0] - temp[0][0] <= max_distance:
		merge_pos_tra(temp, chr, evidence_read, SV_size, svtype, low_bandary)

# def load_sam(sam_path):
# 	starttime = time.time()
# 	samfile = pysam.AlignmentFile(sam_path)
# 	contig_num = len(samfile.get_index_statistics())
# 	print("The total number of chromsomes: %d"%(contig_num))
# 	for _num_ in xrange(contig_num):
# 		Chr_name = samfile.get_reference_name(_num_)
# 		print("Resolving the chromsome %s."%(Chr_name))
# 		for read in samfile.fetch(Chr_name):
# 			parse_read(read, Chr_name, 50)
# 	samfile.close()

# 	for chr in candidate["INS"]:
# 		candidate["INS"][chr] = sorted(candidate["INS"][chr], key = lambda x:x[0])
# 		intergrate_indel(chr, 5, 50, 10, "INS", 50)
# 	for chr in candidate["DEL"]:
# 		candidate["DEL"][chr] = sorted(candidate["DEL"][chr], key = lambda x:x[0])
# 		intergrate_indel(chr, 5, 50, 10, "DEL", 50)

# 	for chr in candidate["DUP"]:
# 		candidate["DUP"][chr] = sorted(candidate["DUP"][chr], key = lambda x:x[0])
# 		intergrate_dup(chr, 5, 50, 10, "l_DUP", "DUP_r", "DUP", 50)
# 	for chr in candidate["INV"]:
# 		candidate["INV"][chr] = sorted(candidate["INV"][chr], key = lambda x:x[0])
# 		intergrate_dup(chr, 5, 50, 10, "l_INV", "INV_r", "INV", 50)

# 	for chr in candidate["TRA"]:
# 		candidate["TRA"][chr] = sorted(candidate["TRA"][chr], key = lambda x:x[0])
# 		intergrate_tra(chr, 5, 50, 10, "TRA", 50)

# 	for chr in result:
# 		result[chr] = sorted(result[chr], key = lambda x:x[0])
# 		for i in result[chr]:
# 			print i,
# 	print("Finished in %0.2f seconds."%(time.time() - starttime))

def single_pipe(sam_path, Chr_name, min_length, min_mapq, max_split_parts, min_seq_size, total_reads, temp_dir):
	samfile = pysam.AlignmentFile(sam_path)
	Chr_length = samfile.get_reference_length(Chr_name)
	logging.info("Resolving the chromsome %s."%(Chr_name))

	for read in samfile.fetch(Chr_name):
		parse_read(read, Chr_name, min_length, min_mapq, max_split_parts, min_seq_size)

	logging.info("%d reads on %s."%(total_reads, Chr_name))
	samfile.close()

	output = "%s_%s.txt"%(temp_dir, Chr_name)
	#print output 
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

	for sv_type in ["l_DUP", "DUP_r", "l_INV", "INV_r"]:
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
	
	samfile = pysam.AlignmentFile(args.input)
	contig_num = len(samfile.get_index_statistics())
	logging.info("The total number of chromsomes: %d"%(contig_num))

	process_list = list()
	for i in samfile.get_index_statistics():
		process_list.append([i[0], i[3]])
		# #chr #read
	process_list = sorted(process_list, key = lambda x:-x[1])
	analysis_pools = Pool(processes=int(args.threads))

	# result = list()
	temporary_dir = args.temp_dir
	for i in process_list:
		para = [(args.input, i[0], args.min_length, args.min_mapq, args.max_split_parts, args.min_seq_size, i[1], temporary_dir)]
		analysis_pools.map_async(multi_run_wrapper, para)
	analysis_pools.close()
	analysis_pools.join()
	samfile.close()

	# for res in result:
		# temp = res.get()[0]

	load_signals(temporary_dir, MainCandidate, process_list)

	# 	for sv_type in ["DEL", "INS", "INV", "DUP", "TRA"]:
	# 		for chr in temp[sv_type]:
	# 			if chr not in MainCandidate[sv_type]:
	# 				MainCandidate[sv_type][chr] = list()
	# 			MainCandidate[sv_type][chr].extend(temp[sv_type][chr])
	# 	for sv_type in ["l_DUP", "DUP_r", "l_INV", "INV_r"]:
	# 		for chr in temp[sv_type]:
	# 			if chr not in MainCandidate[sv_type]:
	# 				MainCandidate[sv_type][chr] = dict()
	# 			for hash_1 in temp[sv_type][chr]:
	# 				if hash_1 not in MainCandidate[sv_type][chr]:
	# 					MainCandidate[sv_type][chr][hash_1] = dict()
	# 				for hash_2 in temp[sv_type][chr][hash_1]:
	# 					if hash_2 not in MainCandidate[sv_type][chr][hash_1]:
	# 						MainCandidate[sv_type][chr][hash_1][hash_2] = list()
	# 					MainCandidate[sv_type][chr][hash_1][hash_2].extend(temp[sv_type][chr][hash_1][hash_2])
	# logging.info("Adjusting SV candidates...")
	show_temp_result(args.min_support, args.min_length, 30, args.max_distance, args.output, MainCandidate)

def load_signals(path, candidate, process_list):
	# file_path = list()
	for i in process_list:
		file_path = "%s_%s.txt"%(path, i[0])
		file = open(file_path, 'r')
		for line in file:
			seq = line.strip('\n').split('\t')
			svtype = seq[0]

			if svtype == "INS" or svtype == "DEL" or svtype == "DUP" or svtype == "INV":
				chr = seq[1]
				start_pos = int(seq[2])
				sv_len = int(seq[3])
				read_name = seq[4]
				if chr not in candidate[svtype]:
					candidate[svtype][chr] = list()
				candidate[svtype][chr].append([start_pos, sv_len, read_name])
			elif svtype == "TRA":
				chr = seq[1]
				start_pos = int(seq[2])
				chr2 = seq[3]
				pos2 = int(seq[4])
				read_name = seq[5]
				if chr not in candidate[svtype]:
					candidate[svtype][chr] = list()
				candidate[svtype][chr].append([start_pos, chr2, pos2, read_name])
			elif svtype == "l_DUP" or svtype == "DUP_r" or svtype == "l_INV" or svtype == "INV_r":
				chr = seq[1]
				key_1 = int(seq[2])
				key_2 = int(seq[3])
				pos = int(seq[4])
				read_name = seq[5]
				if chr not in candidate[svtype]:
					candidate[svtype][chr] = dict()
					candidate[svtype][chr][key_1] = dict()
					candidate[svtype][chr][key_1][key_2] = list()
				else:
					if key_1 not in candidate[svtype][chr]:
						candidate[svtype][chr][key_1] = dict()
						candidate[svtype][chr][key_1][key_2] = list()
					else:
						if key_2 not in candidate[svtype][chr][key_1]:
							candidate[svtype][chr][key_1][key_2] = list()

				candidate[svtype][chr][key_1][key_2].append([pos, read_name])

def show_temp_result(evidence_read, SV_size, low_bandary, max_distance, out_path, MainCandidate):
	# starttime = time.time()
	for chr in MainCandidate["DEL"]:
		MainCandidate["DEL"][chr] = sorted(MainCandidate["DEL"][chr], key = lambda x:x[0])
		# for ele in MainCandidate["DEL"][chr]:
		# 	print chr, ele
		intergrate_indel(chr, evidence_read, SV_size, 200, "DEL", max_distance, MainCandidate)
		MainCandidate["DEL"][chr] = []
		gc.collect()
	# print("[INFO]: Parse deletions used %0.2f seconds."%(time.time() - starttime))

	# starttime = time.time()
	for chr in MainCandidate["INS"]:
		MainCandidate["INS"][chr] = sorted(MainCandidate["INS"][chr], key = lambda x:x[0])
		# for i in MainCandidate["INS"][chr]:
		# 	print chr, i
		intergrate_indel(chr, evidence_read, SV_size, 200, "INS", max_distance, MainCandidate)
		MainCandidate["INS"][chr] = []
		gc.collect()
	# print("[INFO]: Parse insertions used %0.2f seconds."%(time.time() - starttime))

	# starttime = time.time()
	for chr in MainCandidate["DUP"]:
		MainCandidate["DUP"][chr] = sorted(MainCandidate["DUP"][chr], key = lambda x:x[:])
		# for i in MainCandidate["DUP"][chr]:
		# 	print chr, i
		# for i in MainCandidate["l_DUP"][chr]:
		# 	for j in MainCandidate["l_DUP"][chr][i]:
		# 		for k in MainCandidate["l_DUP"][chr][i][j]:
		# 			print k, 0
		# # 	print chr, i
		# for i in MainCandidate["DUP_r"][chr]:
		# 	for j in MainCandidate["DUP_r"][chr][i]:
		# 		for k in MainCandidate["DUP_r"][chr][i][j]:
		# 			print 0, k
		# # 	print chr, i
		intergrate_dup(chr, evidence_read, SV_size, 50, "l_DUP", "DUP_r", "DUP", max_distance, MainCandidate)
		MainCandidate["DUP"][chr] = []
		MainCandidate["l_DUP"][chr] = []
		MainCandidate["DUP_r"][chr] = []
		polish_dup(chr)
		gc.collect()
	# print("[INFO]: Parse duplications used %0.2f seconds."%(time.time() - starttime))
	
	# starttime = time.time()
	for chr in MainCandidate["INV"]:
		MainCandidate["INV"][chr] = sorted(MainCandidate["INV"][chr], key = lambda x:x[:])
		intergrate_dup(chr, evidence_read, SV_size, 10, "l_INV", "INV_r", "INV", max_distance, MainCandidate)
		MainCandidate["INV"][chr] = []
		MainCandidate["l_INV"][chr] = []
		MainCandidate["INV_r"][chr] = []
		gc.collect()
	# print("[INFO]: Parse inversions used %0.2f seconds."%(time.time() - starttime))

	# starttime = time.time()
	for chr in MainCandidate["TRA"]:
		MainCandidate["TRA"][chr] = sorted(MainCandidate["TRA"][chr], key = lambda x:x[:])
		# for i in MainCandidate["TRA"][chr]:
		# 	print chr, i
		intergrate_tra(chr, evidence_read, SV_size, 10, "TRA", max_distance, MainCandidate)
		MainCandidate["TRA"][chr] = []
		gc.collect()
	# print("[INFO]: Parse translocations used %0.2f seconds."%(time.time() - starttime))

	file = open(out_path, 'w')
	# for chr in sorted(semi_result.keys()):
	for chr in semi_result:
		semi_result[chr] = sorted(semi_result[chr], key = lambda x:x[:])
		# result merging
		# exist a bug
		for i in semi_result[chr]:
			# print i
			file.write("%s\t%s"%(chr, i))
	file.close()
	# print candidate['DEL']


def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

def run(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Detection of all types of Structural Variants.
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="cuteSV", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("input", metavar="[BAM]", type=str, help="Sorted .bam file form NGMLR.")
	parser.add_argument('output', type=str, help = "the prefix of novel sequence insertion pridections")
	parser.add_argument('temp_dir', type=str, help = "temporary directory to use for distributed jobs")
	parser.add_argument('-s', '--min_support', help = "Minimum number of reads that support a SV to be reported.[%(default)s]", default = 10, type = int)
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

if __name__ == '__main__':
	run(sys.argv[1:])
