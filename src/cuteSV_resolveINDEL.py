import sys
import numpy as np
from collections import Counter

def resolution_DEL(path, chr, svtype, read_count, threshold_gloab, max_cluster_bias, threshold_local, minimum_support_reads):
	'''
	cluster deletion
	*****************************
	path:	DEL.sigs
	chr:	chromosome id
	svtype:	<DEL>
	
	read_count 	max_cluster_bias 	 	
	----------------------------
	5/10		200 bp (<500 bp)
	*****************************
	'''
	semi_del_cluster = list()
	semi_del_cluster.append([0,0,''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos = int(seq[2])
		indel_len = int(seq[3])
		read_id = seq[4]
		
		if pos - semi_del_cluster[-1][0] > max_cluster_bias:
			if len(semi_del_cluster) >= read_count:
				generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
					threshold_gloab, threshold_local, minimum_support_reads, 
					candidate_single_SV)
			semi_del_cluster = []
			semi_del_cluster.append([pos, indel_len, read_id])
		else:
			semi_del_cluster.append([pos, indel_len, read_id])

	if len(semi_del_cluster) >= read_count:
		generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
			threshold_gloab, threshold_local, minimum_support_reads, 
			candidate_single_SV)
	file.close()
	return candidate_single_SV

def generate_del_cluster_ccs(semi_del_cluster, chr, svtype, read_count, threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):
	
	# calculate support reads
	support_read = len(list(set([i[2] for i in semi_del_cluster])))
	# print(support_read)
	if support_read < read_count:
		return

	del_cluster_len = sorted(semi_del_cluster, key = lambda x:x[1])

	last_len = del_cluster_len[0][1]
	temp_list = dict()
	temp_list[del_cluster_len[0][2]] = 0
	temp_count = 1
	subcluster_bp = list()
	subcluster_bp.append(del_cluster_len[0][0])
	# temp_LEN_sum = last_len
	subcluster_len = list()
	subcluster_len.append(last_len)

	predictions = list()
	for i in del_cluster_len[1:]:
		if last_len < 1000:
			bias = max(20, int(last_len*0.05))
		elif last_len < 5000:
			bias = max(50, int(last_len*0.02))
		elif last_len < 10000:
			bias = max(100,int(last_len*0.02))
		else:
			bias = 200
		if i[1] - last_len > bias:
			if temp_count >= read_count:
				if len(temp_list) >= read_count:
					# breakpoint = round(temp_sum / temp_count)
					breakpoint = Counter(subcluster_bp).most_common(1)[0]
					# del_len = round(temp_LEN_sum / temp_count)
					del_len = Counter(subcluster_len).most_common(1)[0]
					predictions.append([breakpoint, del_len, len(temp_list)])

			temp_list = dict()
			temp_count = 1
			subcluster_bp = list()
			subcluster_bp.append(i[0])
			subcluster_len = list()
			subcluster_len.append(i[1])
			temp_list[i[2]] = 0
		else:
			if i[2] not in temp_list:
				temp_list[i[2]] = 0
			else:
				temp_list[i[2]] += 1
			temp_count += 1
			subcluster_bp.append(i[0])
			subcluster_len.append(i[1])
		last_len = i[1]

	if temp_count >= read_count:
		if len(temp_list) >= read_count:
			breakpoint = Counter(subcluster_bp).most_common(1)[0]
			del_len = Counter(subcluster_len).most_common(1)[0]
			predictions.append([breakpoint, del_len, len(temp_list)])

	predictions = sorted(predictions, key = lambda x:-x[2])

	for i in predictions:
		if i[2] >= read_count and i[2] >= int(support_read / 3):
			if min(i[0][1], i[1][1]) * 2 >= i[2]:
				reliability = "PRECISION"
			else:
				reliability = "IMPRECISION"

			candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % (chr, svtype, i[0][0], -i[1][0], i[2], reliability))



def generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):
	'''
	generate deletion
	*************************************************************
	threshold_gloab 	threshold_local 	minimum_support_reads
	-------------------------------------------------------------
		0.3					0.7 					5
	*************************************************************
	'''
	# Remove duplicates
	read_tag = dict()
	for element in semi_del_cluster:
		if element[2] not in read_tag:
			read_tag[element[2]] = element
		else:
			if element[1] > read_tag[element[2]][1]:
				read_tag[element[2]] = element

	if len(read_tag) < read_count:
		return

	read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])
	global_len = [i[1] for i in read_tag2SortedList]
	DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = threshold_gloab * np.mean(global_len)

	last_len = -DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP
	max_conut = 0
	max_LEN_sum = 0
	max_bps_sum = 0
	temp_count = 0
	temp_LEN_sum = 0
	temp_bps_sum = 0
	for i in read_tag2SortedList:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
			if temp_count > max_conut:
				max_conut = temp_count
				max_LEN_sum = temp_LEN_sum
				max_bps_sum = temp_bps_sum
			temp_count = 1
			temp_LEN_sum = i[1]
			temp_bps_sum = i[0]
		else:
			temp_count += 1
			temp_LEN_sum += i[1]
			temp_bps_sum += i[0]
		last_len = i[1]
	if temp_count > max_conut:
		max_conut = temp_count
		max_LEN_sum = temp_LEN_sum
		max_bps_sum = temp_bps_sum

	breakpoint_starts = round(max_bps_sum / max_conut)
	signal_len = round(max_LEN_sum / max_conut)
	overlap_score = float(max_conut * 1.0 / len(read_tag))

	if max_conut >= minimum_support_reads and overlap_score >= threshold_local:
		# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%d\t%.3f\n'%
		# (chr, svtype, breakpoint_starts, -signal_len, len(read_tag), max_conut, overlap_score))
		candidate_single_SV.append([chr, svtype, str(int(breakpoint_starts)), str(-int(signal_len)), str(len(read_tag))])

def resolution_INS(path, chr, svtype, read_count, threshold_gloab, max_cluster_bias, threshold_local, minimum_support_reads):
	'''
	cluster insertion
	*****************************
	path:	INS.sigs
	chr:	chromosome id
	svtype:	<INS>
	
	read_count 	max_cluster_bias 	 	
	----------------------------
	5/10		100 bp (<500 bp)
	*****************************
	'''
	semi_ins_cluster = list()
	semi_ins_cluster.append([0,0,''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos = int(seq[2])
		indel_len = int(seq[3])
		read_id = seq[4]
		
		if pos - semi_ins_cluster[-1][0] > max_cluster_bias:
			if len(semi_ins_cluster) >= read_count:
				generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV)
			semi_ins_cluster = []
			semi_ins_cluster.append([pos, indel_len, read_id])
		else:
			semi_ins_cluster.append([pos, indel_len, read_id])

	if len(semi_ins_cluster) >= read_count:
		generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV)
	file.close()
	return candidate_single_SV

def generate_ins_cluster_ccs(semi_ins_cluster, chr, svtype, read_count, threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):

	# calculate support reads
	support_read = len(list(set([i[2] for i in semi_ins_cluster])))
	if support_read < read_count:
		return

	del_cluster_len = sorted(semi_ins_cluster, key = lambda x:x[1])

	last_len = del_cluster_len[0][1]
	temp_list = dict()
	temp_list[del_cluster_len[0][2]] = 0
	temp_count = 1
	subcluster_bp = list()
	subcluster_bp.append(del_cluster_len[0][0])
	# temp_LEN_sum = last_len
	subcluster_len = list()
	subcluster_len.append(last_len)

	predictions = list()
	for i in del_cluster_len[1:]:
		if last_len < 1000:
			bias = max(20, int(last_len*0.05))
		elif last_len < 5000:
			bias = max(50, int(last_len*0.02))
		elif last_len < 10000:
			bias = max(100,int(last_len*0.02))
		else:
			bias = 200
		if i[1] - last_len > bias:
			if temp_count >= read_count:
				if len(temp_list) >= read_count:
					# breakpoint = round(temp_sum / temp_count)
					breakpoint = Counter(subcluster_bp).most_common(1)[0]
					# del_len = round(temp_LEN_sum / temp_count)
					del_len = Counter(subcluster_len).most_common(1)[0]
					predictions.append([breakpoint, del_len, len(temp_list)])

			temp_list = dict()
			temp_count = 1
			subcluster_bp = list()
			subcluster_bp.append(i[0])
			subcluster_len = list()
			subcluster_len.append(i[1])
			temp_list[i[2]] = 0
		else:
			if i[2] not in temp_list:
				temp_list[i[2]] = 0
			else:
				temp_list[i[2]] += 1
			temp_count += 1
			subcluster_bp.append(i[0])
			subcluster_len.append(i[1])
		last_len = i[1]

	if temp_count >= read_count:
		if len(temp_list) >= read_count:
			breakpoint = Counter(subcluster_bp).most_common(1)[0]
			del_len = Counter(subcluster_len).most_common(1)[0]
			predictions.append([breakpoint, del_len, len(temp_list)])

	predictions = sorted(predictions, key = lambda x:-x[2])

	for i in predictions:
		# print(i)
		# print(support_read)
		# if i[2] >= read_count and i[2] >= int(support_read / 3):
		if i[2] >= read_count:
			if min(i[0][1], i[1][1]) * 2 >= i[2]:
				reliability = "PRECISION"
			else:
				reliability = "IMPRECISION"

			if i[0][1] == 1 and i[1][1] == 1:
				return

			candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % (chr, svtype, i[0][0], i[1][0], i[2], reliability))


def generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):
	'''
	generate insertion
	*************************************************************
	threshold_gloab 	threshold_local 	minimum_support_reads
	-------------------------------------------------------------
		0.2					0.6						5
	*************************************************************
	'''
	# Remove duplicates
	read_tag = dict()
	for element in semi_ins_cluster:
		if element[2] not in read_tag:
			read_tag[element[2]] = element
		else:
			if element[1] > read_tag[element[2]][1]:
				read_tag[element[2]] = element

	if len(read_tag) < read_count:
		return

	read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])
	# start&end breakpoint
	breakpoint_starts = [i[0] for i in read_tag2SortedList]
	global_len = [i[1] for i in read_tag2SortedList]

	DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)

	last_len = -DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP
	max_conut = 0
	max_LEN_sum = 0
	temp_count = 0
	temp_LEN_sum = 0
	for i in read_tag2SortedList:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
			if temp_count > max_conut:
				max_conut = temp_count
				max_LEN_sum = temp_LEN_sum
			temp_count = 1
			temp_LEN_sum = i[1]
		else:
			temp_count += 1
			temp_LEN_sum += i[1]
		last_len = i[1]
	if temp_count > max_conut:
		max_conut = temp_count
		max_LEN_sum = temp_LEN_sum

	breakpointStart = np.mean(breakpoint_starts)
	breakpointStart_STD = np.std(breakpoint_starts)
	signal_len = round(max_LEN_sum / max_conut)
	overlap_score = float(max_conut * 1.0 / len(read_tag))

	if max_conut >= minimum_support_reads and overlap_score >= threshold_local:
		# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%d\n' % 
		# 	(chr, svtype, breakpointStart, signal_len, len(read_tag), max_conut, overlap_score, breakpointStart_STD, np.mean(global_len)))
		candidate_single_SV.append([chr, svtype, str(int(breakpointStart)), str(int(signal_len)), str(len(read_tag))])

def run_del(args):
	return resolution_DEL(*args)

def run_ins(args):
	return resolution_INS(*args)
