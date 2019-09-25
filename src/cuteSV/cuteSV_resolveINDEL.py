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
*******************************************

'''

def resolution_DEL(path, chr, svtype, read_count, threshold_gloab, max_cluster_bias, 
	threshold_local, minimum_support_reads):

	'''
	cluster DEL
	********************************************************************************************
	path:	DEL.sigs
	chr:	chromosome id
	svtype:	<DEL>
	
	SEQTYPE		read_count 	max_cluster_bias 	sv_size		threshold_gloab 	threshold_local 
	--------------------------------------------------------------------------------------------
	CCS			3			200 bp (<500 bp)	30 bp 		0.4					0.5
	CLR			5/10		200 bp (<500 bp)	50 bp 		0.3					0.7
	--------------------------------------------------------------------------------------------
	
	Input file format
	--------------------------------------------------------------------------------------------
	column	#1	#2	#3	#4	#5
			DEL	CHR	BP	LEN	ID	
	#1	deletion type
	#2	chromosome number
	#3	breakpoint in each read
	#4	DEL_len in each read
	#5	read ID
	********************************************************************************************
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
				generate_del_cluster_2(semi_del_cluster, 
										chr, 
										svtype, 
										read_count, 
										threshold_gloab, 
										threshold_local, 
										minimum_support_reads, 
										candidate_single_SV)
			semi_del_cluster = []
			semi_del_cluster.append([pos, indel_len, read_id])
		else:
			semi_del_cluster.append([pos, indel_len, read_id])

	if len(semi_del_cluster) >= read_count:
		generate_del_cluster_2(semi_del_cluster, 
								chr, 
								svtype, 
								read_count, 
								threshold_gloab, 
								threshold_local, 
								minimum_support_reads, 
								candidate_single_SV)
	file.close()
	return candidate_single_SV

def generate_del_cluster_ccs(semi_del_cluster, chr, svtype, read_count, threshold_gloab,
 threshold_local, minimum_support_reads, candidate_single_SV):
	
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

			candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % 
				(chr, svtype, i[0][0], -i[1][0], i[2], reliability))



def generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):

	'''
	generate deletion
	*************************************************************
	threshold_gloab 	threshold_local 	minimum_support_reads
	-------------------------------------------------------------
		0.3					0.7 					5		CLR
		0.4					0.5 				  <=5		CCS
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
		# 	(chr, svtype, breakpoint_starts, -signal_len, len(read_tag), max_conut, overlap_score))
		candidate_single_SV.append([chr, 
									svtype, 
									str(int(breakpoint_starts)), 
									str(-int(signal_len)), 
									str(len(read_tag))])

def generate_del_cluster_2(semi_del_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):

	'''
	generate deletion
	*************************************************************
	threshold_gloab 	threshold_local 	minimum_support_reads
	-------------------------------------------------------------
		0.3					0.7 					5		CLR
		0.4					0.5 				  <=5		CCS
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

	last_len = read_tag2SortedList[0][1]

	alelle_collect = list()
	alelle_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[]])

	for i in read_tag2SortedList[1:]:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
			alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
			alelle_collect.append([[],[],[]])

		alelle_collect[-1][0].append(i[0])
		alelle_collect[-1][1].append(i[1])
		last_len = i[1]
	alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
	alelle_sort = sorted(alelle_collect, key = lambda x:x[2])

	if alelle_sort[-1][2][0] >= minimum_support_reads and float(alelle_sort[-1][2][0] * 1.0 / len(read_tag)) >= threshold_local:
		breakpointStart = np.mean(alelle_sort[-1][0])
		breakpointStart_STD = np.std(alelle_sort[-1][0])
		signalLen = np.mean(alelle_sort[-1][1])
		signalLen_STD = np.std(alelle_sort[-1][1])
		# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
		candidate_single_SV.append([chr, 
									svtype, 
									str(int(breakpointStart)), 
									str(int(-signalLen)), 
									str(alelle_sort[-1][2][0]), 
									"%.3f"%breakpointStart_STD, 
									"%.3f"%signalLen_STD])

		# extend to next alelle
		if (len(alelle_sort) > 1 and alelle_sort[-2][2][0] >= minimum_support_reads 
			and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag) 
			and alelle_sort[-2][2][0] >= 0.3*len(read_tag)):
			breakpointStart = np.mean(alelle_sort[-2][0])
			breakpointStart_STD = np.std(alelle_sort[-2][0])
			signalLen = np.mean(alelle_sort[-2][1])
			last_signalLen_STD = signalLen_STD
			signalLen_STD = np.std(alelle_sort[-2][1])
			# if signalLen_STD < last_signalLen_STD:
			# 	# pass
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(-signalLen)), 
										str(alelle_sort[-2][2][0]), 
										"%.3f"%breakpointStart_STD, 
										"%.3f"%signalLen_STD])

	elif alelle_sort[-2][2][0] >= minimum_support_reads and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag):
		if alelle_sort[-2][2][0] >= 0.4*len(read_tag):
			breakpointStart = np.mean(alelle_sort[-1][0])
			breakpointStart_STD = np.std(alelle_sort[-1][0])
			signalLen = np.mean(alelle_sort[-1][1])
			signalLen_STD = np.std(alelle_sort[-1][1])
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(-signalLen)), 
										str(alelle_sort[-1][2][0]), 
										"%.3f"%breakpointStart_STD, 
										"%.3f"%signalLen_STD])

			breakpointStart = np.mean(alelle_sort[-2][0])
			breakpointStart_STD = np.std(alelle_sort[-2][0])
			signalLen = np.mean(alelle_sort[-2][1])
			signalLen_STD = np.std(alelle_sort[-2][1])
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(-signalLen)), 
										str(alelle_sort[-2][2][0]), 
										"%.3f"%breakpointStart_STD, 
										"%.3f"%signalLen_STD])


	

def resolution_INS(path, chr, svtype, read_count, threshold_gloab, 
	max_cluster_bias, threshold_local, minimum_support_reads):
	
	'''
	cluster INS
	********************************************************************************************
	path:	INS.sigs
	chr:	chromosome id
	svtype:	<INS>
	
	SEQTYPE		read_count 	max_cluster_bias 	sv_size		threshold_gloab 	threshold_local 
	--------------------------------------------------------------------------------------------
	CCS			3			200 bp (<500 bp)	30 bp 		0.65				0.7
	CLR			5/10		100 bp (<500 bp)	50 bp 		0.2					0.6
	--------------------------------------------------------------------------------------------
	
	Input file format
	--------------------------------------------------------------------------------------------
	column	#1	#2	#3	#4	#5
			INS	CHR	BP	LEN	ID	
	#1	insertion type
	#2	chromosome number
	#3	breakpoint in each read
	#4	DEL_len in each read
	#5	read ID
	********************************************************************************************
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
				generate_ins_cluster_2(semi_ins_cluster, 
										chr, 
										svtype, 
										read_count, 
										threshold_gloab, 
										threshold_local, 
										minimum_support_reads, 
										candidate_single_SV)
			semi_ins_cluster = []
			semi_ins_cluster.append([pos, indel_len, read_id])
		else:
			semi_ins_cluster.append([pos, indel_len, read_id])

	if len(semi_ins_cluster) >= read_count:
		generate_ins_cluster_2(semi_ins_cluster, 
								chr, 
								svtype, 
								read_count, 
								threshold_gloab, 
								threshold_local, 
								minimum_support_reads, 
								candidate_single_SV)
	file.close()
	return candidate_single_SV

def generate_ins_cluster_ccs(semi_ins_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):

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

			candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % 
				(chr, svtype, i[0][0], i[1][0], i[2], reliability))


def generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):
		
	'''
	generate deletion
	*************************************************************
	threshold_gloab 	threshold_local 	minimum_support_reads
	-------------------------------------------------------------
		0.2					0.6 					5		CLR
		0.65				0.7 				  <=5		CCS
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
		candidate_single_SV.append([chr, 
									svtype, 
									str(int(breakpointStart)), 
									str(int(signal_len)), 
									str(len(read_tag))])

def generate_ins_cluster_2(semi_ins_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV):
		
	'''
	generate deletion
	*************************************************************
	threshold_gloab 	threshold_local 	minimum_support_reads
	-------------------------------------------------------------
		0.2					0.6 					5		CLR
		0.65				0.7 				  <=5		CCS
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
	global_len = [i[1] for i in read_tag2SortedList]
	DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP = threshold_gloab * np.mean(global_len)

	last_len = read_tag2SortedList[0][1]

	alelle_collect = list()
	alelle_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[]])

	for i in read_tag2SortedList[1:]:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
			alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
			alelle_collect.append([[],[],[]])

		alelle_collect[-1][0].append(i[0])
		alelle_collect[-1][1].append(i[1])
		last_len = i[1]
	alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
	alelle_sort = sorted(alelle_collect, key = lambda x:x[2])

	# print(len(read_tag), DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP)
	# for i in alelle_sort:
	# 	for j in i:
	# 		print(j)

	# print(alelle_sort[-1][2][0], minimum_support_reads)
	# print(float(alelle_sort[-1][2][0] * 1.0 / len(read_tag)), threshold_local)

	if alelle_sort[-1][2][0] >= minimum_support_reads and float(alelle_sort[-1][2][0] * 1.0 / len(read_tag)) >= threshold_local:
		breakpointStart = np.mean(alelle_sort[-1][0])
		breakpointStart_STD = np.std(alelle_sort[-1][0])
		signalLen = np.mean(alelle_sort[-1][1])
		signalLen_STD = np.std(alelle_sort[-1][1])
		# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
		candidate_single_SV.append([chr, 
									svtype, 
									str(int(breakpointStart)), 
									str(int(signalLen)), 
									str(alelle_sort[-1][2][0]), 
									"%.3f"%breakpointStart_STD, 
									"%.3f"%signalLen_STD])

		# extend to next alelle
		if (len(alelle_sort) > 1 and alelle_sort[-2][2][0] >= minimum_support_reads 
			and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag) 
			and alelle_sort[-2][2][0] >= 0.3*len(read_tag)):
			breakpointStart = np.mean(alelle_sort[-2][0])
			breakpointStart_STD = np.std(alelle_sort[-2][0])
			signalLen = np.mean(alelle_sort[-2][1])
			last_signalLen_STD = signalLen_STD
			signalLen_STD = np.std(alelle_sort[-2][1])
			if signalLen_STD < last_signalLen_STD:
				# pass
				# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
				candidate_single_SV.append([chr, 
											svtype, 
											str(int(breakpointStart)), 
											str(int(signalLen)), 
											str(alelle_sort[-2][2][0]), 
											"%.3f"%breakpointStart_STD, 
											"%.3f"%signalLen_STD])

	elif alelle_sort[-2][2][0] >= minimum_support_reads and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag):
		if alelle_sort[-2][2][0] >= 0.4*len(read_tag):
			breakpointStart = np.mean(alelle_sort[-1][0])
			breakpointStart_STD = np.std(alelle_sort[-1][0])
			signalLen = np.mean(alelle_sort[-1][1])
			signalLen_STD = np.std(alelle_sort[-1][1])
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(signalLen)), 
										str(alelle_sort[-1][2][0]), 
										"%.3f"%breakpointStart_STD, 
										"%.3f"%signalLen_STD])

			breakpointStart = np.mean(alelle_sort[-2][0])
			breakpointStart_STD = np.std(alelle_sort[-2][0])
			signalLen = np.mean(alelle_sort[-2][1])
			signalLen_STD = np.std(alelle_sort[-2][1])
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(signalLen)), 
										str(alelle_sort[-2][2][0]), 
										"%.3f"%breakpointStart_STD, 
										"%.3f"%signalLen_STD])

def run_del(args):
	return resolution_DEL(*args)

def run_ins(args):
	return resolution_INS(*args)
