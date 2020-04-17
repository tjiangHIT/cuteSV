import sys
import numpy as np
from collections import Counter

'''
*******************************************
				TO DO LIST
*******************************************
	1. Identify DP with samfile pointer;
	2. Add CIPOS, CILEN and/or CIEND;
	3. Determine (IM)PRECISE type.
	4. Filter DUP to improve INS FN rate.
*******************************************
'''

def resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size, 
	bam_path, action, hom, het, MaxSize):
	semi_dup_cluster = list()
	semi_dup_cluster.append([0,0,''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos_1 = int(seq[2])
		pos_2 = int(seq[3])
		read_id = seq[4]
		
		if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias:
			if len(semi_dup_cluster) >= read_count:
				generate_dup_cluster(semi_dup_cluster, 
										chr, 
										read_count, 
										max_cluster_bias, 
										sv_size, 
										candidate_single_SV,
										bam_path,
										action,
										hom,
										het,
										MaxSize)
			semi_dup_cluster = []
			semi_dup_cluster.append([pos_1, pos_2, read_id])
		else:
			semi_dup_cluster.append([pos_1, pos_2, read_id])

	if len(semi_dup_cluster) >= read_count:
		generate_dup_cluster(semi_dup_cluster, 
								chr, read_count, 
								max_cluster_bias, 
								sv_size, 
								candidate_single_SV,
								bam_path,
								action,
								hom,
								het,
								MaxSize)
	file.close()
	return candidate_single_SV

def generate_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, 
	sv_size, candidate_single_SV, bam_path, action, hom, het, MaxSize):
	# calculate support reads
	support_read = list(set([i[2] for i in semi_dup_cluster]))
	# print(support_read)
	if len(support_read) < read_count:
		return

	dup_cluster_len = sorted(semi_dup_cluster, key = lambda x:x[1])
	# for i in dup_cluster_len:
	# 	print(i)

	last_len = dup_cluster_len[0][1]
	temp_list = dict()
	temp_list[dup_cluster_len[0][2]] = 0
	temp_count = 1
	subcluster_bp1 = list()
	subcluster_bp1.append(dup_cluster_len[0][0])
	# temp_LEN_sum = last_len
	subcluster_bp2 = list()
	subcluster_bp2.append(last_len)

	predictions = list()
	for i in dup_cluster_len[1:]:
		if i[1] - last_len > max_cluster_bias:
			if temp_count >= read_count:
				if len(temp_list) >= read_count:
					# breakpoint = round(temp_sum / temp_count)
					breakpoint_1 = Counter(subcluster_bp1).most_common(1)[0]
					# del_len = round(temp_LEN_sum / temp_count)
					breakpoint_2 = Counter(subcluster_bp2).most_common(1)[0]
					predictions.append([breakpoint_1, breakpoint_2, len(temp_list)])
					# print(breakpoint_1, breakpoint_2, len(temp_list))

			temp_list = dict()
			temp_count = 1
			subcluster_bp1 = list()
			subcluster_bp1.append(i[0])
			subcluster_bp2 = list()
			subcluster_bp2.append(i[1])
			temp_list[i[2]] = 0
		else:
			if i[2] not in temp_list:
				temp_list[i[2]] = 0
			else:
				temp_list[i[2]] += 1
			temp_count += 1
			subcluster_bp1.append(i[0])
			subcluster_bp2.append(i[1])
		last_len = i[1]

	if temp_count >= read_count:
		if len(temp_list) >= read_count:
			breakpoint_1 = Counter(subcluster_bp1).most_common(1)[0]
			breakpoint_2 = Counter(subcluster_bp2).most_common(1)[0]
			predictions.append([breakpoint_1, breakpoint_2, len(temp_list)])

	predictions = sorted(predictions, key = lambda x:-x[2])

	for i in predictions:
		# print(i)
		if i[2] >= read_count and i[2] >= int(len(support_read) / 3) and i[1][0] - i[0][0] >= sv_size:
			if min(i[0][1], i[1][1]) * 2 >= i[2]:
				reliability = "PRECISION"
			else:
				reliability = "IMPRECISION"

			# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % (chr, 'DUP', i[0][0], i[1][0] - i[0][0], i[2], reliability))
			'''genotyping'''
			if i[1][0] - i[0][0] <= MaxSize:
				if action:
					DV, DR, GT = call_gt(bam_path, int(i[0][0]), int(i[1][0]), 
						chr, support_read, max_cluster_bias, hom, het)
				else:
					DR = '.'
					GT = './.'
				candidate_single_SV.append([chr,
											'DUP', 
											str(i[0][0]), 
											str(i[1][0] - i[0][0]), 
											str(i[2]),
											str(DR),
											str(GT)])

def run_dup(args):
	return resolution_DUP(*args)

def count_coverage(chr, s, e, f):
	read_count = set()
	for i in f.fetch(chr, s, e):
		# read_count.add(i.query_name)
		if i.flag not in [0,16]:
			continue
		if i.reference_start < s and i.reference_end > e:
			read_count.add(i.query_name)
	return read_count

def assign_gt(a, b, hom, het):
	if b == 0:
		return "1/1"
	if a*1.0/b < het:
		return "0/0"
	elif a*1.0/b >= het and a*1.0/b < hom:
		return "0/1"
	elif a*1.0/b >= hom and a*1.0/b < 1.0:
		return "1/1"
	else:
		return "1/1"

def call_gt(bam_path, pos_1, pos_2, chr, read_id_list, max_cluster_bias, hom, het):
	import pysam
	bamfile = pysam.AlignmentFile(bam_path)
	search_start = max(int(pos_1) - max_cluster_bias, 0)
	search_end = min(int(pos_2) + max_cluster_bias, bamfile.get_reference_length(chr))
	querydata = count_coverage(chr, search_start, search_end, bamfile)
	bamfile.close()
	DR = 0
	for query in querydata:
		if query not in read_id_list:
			DR += 1
	return len(read_id_list), DR, assign_gt(len(read_id_list), DR+len(read_id_list), hom, het)
