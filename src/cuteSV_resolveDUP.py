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

def resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size):
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
				generate_dup_cluster_ccs(semi_dup_cluster, 
										chr, 
										read_count, 
										max_cluster_bias, 
										sv_size, 
										candidate_single_SV)
			semi_dup_cluster = []
			semi_dup_cluster.append([pos_1, pos_2, read_id])
		else:
			semi_dup_cluster.append([pos_1, pos_2, read_id])

	if len(semi_dup_cluster) >= read_count:
		generate_dup_cluster_ccs(semi_dup_cluster, 
								chr, read_count, 
								max_cluster_bias, 
								sv_size, 
								candidate_single_SV)
	file.close()
	return candidate_single_SV

def generate_dup_cluster_ccs(semi_dup_cluster, chr, read_count, max_cluster_bias, sv_size, candidate_single_SV):
	# calculate support reads
	support_read = len(list(set([i[2] for i in semi_dup_cluster])))
	# print(support_read)
	if support_read < read_count:
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
		if i[2] >= read_count and i[2] >= int(support_read / 3) and i[1][0] - i[0][0] >= sv_size:
			if min(i[0][1], i[1][1]) * 2 >= i[2]:
				reliability = "PRECISION"
			else:
				reliability = "IMPRECISION"

			# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % (chr, 'DUP', i[0][0], i[1][0] - i[0][0], i[2], reliability))
			candidate_single_SV.append([chr, 'DUP', str(i[0][0]), str(i[1][0] - i[0][0]), str(i[2])])


def generate_semi_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, sv_size, dup_candidates, candidate_single_SV):
	'''
	generate DUP
	*******************************************
	overlap_size 	max_cluster_bias 	sv_size
	-------------------------------------------
	0.5				50 bp (<500 bp)		50 bp
	*******************************************
	'''
	start = list()
	end = list()
	start_id = list()
	end_id = list()
	read_tag = dict()
	for element in semi_dup_cluster:
		start.append(element[0])
		end.append(element[1])
		start_id.append([element[0], element[2]])
		end_id.append([element[1], element[2]])
		# if element[2] not in read_tag:
		# 	read_tag[element[2]] = 0

	# re_l, rl_id = acquire_locus(int(min(start)/1000)*1000, int(max(start)/1000+1)*1000, "l_DUP", chr, dup_candidates)
	# re_r, rr_id = acquire_locus(min(end)-1000, max(end)+1000, "DUP_r", chr, dup_candidates)
	re_l = acquire_locus(min(start)-max_cluster_bias, max(start)+max_cluster_bias, "l_DUP", chr, dup_candidates)
	re_r = acquire_locus(min(end)-max_cluster_bias, max(end)+max_cluster_bias, "DUP_r", chr, dup_candidates)
	info_1 = cal_cluster_breakpoint(start_id+re_l, max_cluster_bias)
	info_2 = cal_cluster_breakpoint(end_id+re_r, max_cluster_bias)

	cal_pos = 0
	for i in info_1:
		if i[1] not in read_tag:
			read_tag[i[1]] = 0
		cal_pos += i[0]
	breakpoint_1 = int(cal_pos/len(info_1))
	cal_pos = 0
	for i in info_2:
		if i[1] not in read_tag:
			read_tag[i[1]] = 0
		cal_pos += i[0]
	breakpoint_2 = int(cal_pos/len(info_2))	

	if breakpoint_2 - breakpoint_1 + 1 >= sv_size and len(read_tag) >= read_count:
		# print info_1
		# print info_2
		# print("%s\t%s\t%d\t%d\t%d"%(chr, "DUP", breakpoint_1, breakpoint_2-breakpoint_1+1, len(read_tag)))
		# candidate_single_SV.append("%s\t%s\t%d\t%d\t%d\n"%(chr, "DUP", breakpoint_1, breakpoint_2-breakpoint_1+1, len(read_tag)))
		candidate_single_SV.append([chr, "DUP", str(breakpoint_1), str(breakpoint_2-breakpoint_1+1), str(len(read_tag))])

def acquire_locus(down, up, keytype, chr, MainCandidate):
	re = list()
	# re_id = list()
	if chr not in MainCandidate[keytype]:
		return re
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in MainCandidate[keytype][chr]:
			return re
		for i in range(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					re.append([ele[0],ele[1]])
					# re_id.append(ele[1])
	else:
		key_1 = int(down/10000)
		if key_1 in MainCandidate[keytype][chr]:
			for i in range(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in MainCandidate[keytype][chr][key_1]:
					continue
				for ele in MainCandidate[keytype][chr][key_1][key_2]:
					if ele[0] >= down and ele[0] <= up:
						re.append([ele[0],ele[1]])
						# re_id.append(ele[1])
		key_1 += 1
		if key_1 not in MainCandidate[keytype][chr]:
			return re
		for i in range(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					re.append([ele[0],ele[1]])
					# re_id.append(ele[1])
	return re

def cal_cluster_breakpoint(pos_list, bias):
	pos_list = sorted(pos_list, key = lambda x:x[0])
	#print pos_list
	temp = list()
	temp.append(pos_list[0])
	best_bp =list()
	len_best_bp = 0
	for pos in pos_list[1:]:
		if temp[-1][0] + bias < pos[0]:
			if len(temp) > len_best_bp:
				len_best_bp = len(temp)
				# best_bp = sum(temp)/len(temp)
				best_bp = temp
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if len(temp) > len_best_bp:
		len_best_bp = len(temp)
		# best_bp = sum(temp)/len(temp)
		best_bp = temp
	return best_bp

def run_dup(args):
	return resolution_DUP(*args)
