import sys
import numpy as np
from collections import Counter
from cuteSV.cuteSV_genotype import cal_GL, cal_CIPOS, threshold_ref_count, count_coverage
import time

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
	threshold_local, minimum_support_reads, bam_path, action, gt_round):

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
				if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
					pass
				else:
					generate_del_cluster(semi_del_cluster, 
										chr, 
										svtype, 
										read_count, 
										threshold_gloab, 
										threshold_local, 
										minimum_support_reads, 
										candidate_single_SV,
										bam_path,
										max_cluster_bias,
										action,
										gt_round)
			semi_del_cluster = []
			semi_del_cluster.append([pos, indel_len, read_id])
		else:
			semi_del_cluster.append([pos, indel_len, read_id])

	if len(semi_del_cluster) >= read_count:
		if semi_del_cluster[-1][0] == semi_del_cluster[-1][1] == 0:
			pass
		else:
			generate_del_cluster(semi_del_cluster, 
								chr, 
								svtype, 
								read_count, 
								threshold_gloab, 
								threshold_local, 
								minimum_support_reads, 
								candidate_single_SV,
								bam_path,
								max_cluster_bias,
								action,
								gt_round)
	file.close()
	return candidate_single_SV

def generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV, 
	bam_path, max_cluster_bias, action, gt_round):

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
	'''
	*************************************************************
		#1 				#2			#3			#4
	-------------------------------------------------------------
		del-breakpoint	del-len		#support 	read-id
	*************************************************************
	'''
	alelle_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[],
		[read_tag2SortedList[0][2]]])

	for i in read_tag2SortedList[1:]:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
			alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
			alelle_collect.append([[],[],[],[]])

		alelle_collect[-1][0].append(i[0])
		alelle_collect[-1][1].append(i[1])
		alelle_collect[-1][3].append(i[2])
		last_len = i[1]
	alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
	alelle_sort = sorted(alelle_collect, key = lambda x:x[2])

	if alelle_sort[-1][2][0] >= minimum_support_reads and float(alelle_sort[-1][2][0] * 1.0 / len(read_tag)) >= threshold_local:
		breakpointStart = np.mean(alelle_sort[-1][0])
		# breakpointStart_STD = np.std(alelle_sort[-1][0])
		CIPOS = cal_CIPOS(np.std(alelle_sort[-1][0]), len(alelle_sort[-1][0]))
		search_threshold = np.min(alelle_sort[-1][0])
		signalLen = np.mean(alelle_sort[-1][1])
		CILEN = cal_CIPOS(np.std(alelle_sort[-1][1]), len(alelle_sort[-1][1]))
		signalLen_STD = np.std(alelle_sort[-1][1])
		# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)

		'''genotyping'''
		if action:
			# time_start = time.time()
			DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, search_threshold, chr, alelle_sort[-1][3], 
												max_cluster_bias, gt_round)
			# cost_time = time.time() - time_start
			# print("DEL", chr, int(breakpointStart), int(-signalLen), DR, DV, QUAL, "%.4f"%cost_time)
		else:
			DR = '.'
			GT = './.'
			GL = '.,.,.'
			GQ = "."
			QUAL = "."
		# print(DV, DR, GT)
		candidate_single_SV.append([chr, 
									svtype, 
									str(int(breakpointStart)), 
									str(int(-signalLen)), 
									str(alelle_sort[-1][2][0]), 
									str(CIPOS),
									str(CILEN),
									str(DR),
									str(GT),
									str(GL),
									str(GQ),
									str(QUAL),
									str(','.join(alelle_sort[-1][3]))])

		# extend to next alelle
		if (len(alelle_sort) > 1 and alelle_sort[-2][2][0] >= minimum_support_reads 
			and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag) 
			and alelle_sort[-2][2][0] >= 0.3*len(read_tag)):
			breakpointStart = np.mean(alelle_sort[-2][0])
			# breakpointStart_STD = np.std(alelle_sort[-2][0])
			CIPOS = cal_CIPOS(np.std(alelle_sort[-2][0]), len(alelle_sort[-2][0]))
			search_threshold = np.min(alelle_sort[-2][0])
			signalLen = np.mean(alelle_sort[-2][1])
			last_signalLen_STD = signalLen_STD
			signalLen_STD = np.std(alelle_sort[-2][1])
			CILEN = cal_CIPOS(np.std(alelle_sort[-2][1]), len(alelle_sort[-2][1]))
			# if signalLen_STD < last_signalLen_STD:
			# 	# pass
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
			'''genotyping'''
			if action:
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, search_threshold, chr, alelle_sort[-2][3], 
													max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("DEL", chr, int(breakpointStart), int(-signalLen), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(-signalLen)), 
										str(alelle_sort[-2][2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(alelle_sort[-2][3]))])

	elif alelle_sort[-2][2][0] >= minimum_support_reads and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag):
		if alelle_sort[-2][2][0] >= 0.4*len(read_tag):
			breakpointStart = np.mean(alelle_sort[-1][0])
			# breakpointStart_STD = np.std(alelle_sort[-1][0])
			CIPOS = cal_CIPOS(np.std(alelle_sort[-1][0]), len(alelle_sort[-1][0]))
			search_threshold = np.min(alelle_sort[-1][0])
			signalLen = np.mean(alelle_sort[-1][1])
			signalLen_STD = np.std(alelle_sort[-1][1])
			CILEN = cal_CIPOS(np.std(alelle_sort[-1][1]), len(alelle_sort[-1][1]))
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
			'''genotyping'''
			if action:
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, search_threshold, chr, alelle_sort[-1][3], 
													max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("DEL", chr, int(breakpointStart), int(-signalLen), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(-signalLen)), 
										str(alelle_sort[-1][2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(alelle_sort[-1][3]))])

			breakpointStart = np.mean(alelle_sort[-2][0])
			# breakpointStart_STD = np.std(alelle_sort[-2][0])
			CIPOS = cal_CIPOS(np.std(alelle_sort[-2][0]), len(alelle_sort[-2][0]))
			search_threshold = np.min(alelle_sort[-2][0])
			signalLen = np.mean(alelle_sort[-2][1])
			signalLen_STD = np.std(alelle_sort[-2][1])
			CILEN = cal_CIPOS(np.std(alelle_sort[-2][1]), len(alelle_sort[-2][1]))
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
			'''genotyping'''
			if action:
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, search_threshold, chr, alelle_sort[-2][3], 
												max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("DEL", chr, int(breakpointStart), int(-signalLen), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(-signalLen)), 
										str(alelle_sort[-2][2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(alelle_sort[-2][3]))])


	

def resolution_INS(path, chr, svtype, read_count, threshold_gloab, 
	max_cluster_bias, threshold_local, minimum_support_reads, bam_path, action, gt_round):
	
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
				if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
					pass
				else:
					generate_ins_cluster(semi_ins_cluster, 
										chr, 
										svtype, 
										read_count, 
										threshold_gloab, 
										threshold_local, 
										minimum_support_reads, 
										candidate_single_SV,
										bam_path,
										max_cluster_bias,
										action,
										gt_round)
			semi_ins_cluster = []
			semi_ins_cluster.append([pos, indel_len, read_id])
		else:
			semi_ins_cluster.append([pos, indel_len, read_id])

	if len(semi_ins_cluster) >= read_count:
		if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
			pass
		else:
			generate_ins_cluster(semi_ins_cluster, 
								chr, 
								svtype, 
								read_count, 
								threshold_gloab, 
								threshold_local, 
								minimum_support_reads, 
								candidate_single_SV,
								bam_path,
								max_cluster_bias,
								action,
								gt_round)
	file.close()
	return candidate_single_SV

def generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, 
	threshold_gloab, threshold_local, minimum_support_reads, candidate_single_SV, 
	bam_path, max_cluster_bias, action, gt_round):
		
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
	alelle_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[],
		[read_tag2SortedList[0][2]]])

	for i in read_tag2SortedList[1:]:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
			alelle_collect[-1][2].append(len(alelle_collect[-1][0]))
			alelle_collect.append([[],[],[],[]])

		alelle_collect[-1][0].append(i[0])
		alelle_collect[-1][1].append(i[1])
		alelle_collect[-1][3].append(i[2])
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
		# breakpointStart_STD = np.std(alelle_sort[-1][0])
		CIPOS = cal_CIPOS(np.std(alelle_sort[-1][0]), len(alelle_sort[-1][0]))
		signalLen = np.mean(alelle_sort[-1][1])
		signalLen_STD = np.std(alelle_sort[-1][1])
		CILEN = cal_CIPOS(np.std(alelle_sort[-1][1]), len(alelle_sort[-1][1]))
		# search_threshold = np.min(alelle_sort[-1][0])
		# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
		'''genotyping'''
		if action:
			# time_start = time.time()
			DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(breakpointStart), chr, alelle_sort[-1][3], 
											max_cluster_bias, gt_round)
			# cost_time = time.time() - time_start
			# print("INS", chr, int(breakpointStart), int(signalLen), DR, DV, QUAL, "%.4f"%cost_time)
		else:
			DR = '.'
			GT = './.'
			GL = '.,.,.'
			GQ = "."
			QUAL = "."
		candidate_single_SV.append([chr, 
									svtype, 
									str(int(breakpointStart)), 
									str(int(signalLen)), 
									str(alelle_sort[-1][2][0]), 
									str(CIPOS),
									str(CILEN),
									str(DR),
									str(GT),
									str(GL),
									str(GQ),
									str(QUAL),
									str(','.join(alelle_sort[-1][3]))])

		# extend to next alelle
		if (len(alelle_sort) > 1 and alelle_sort[-2][2][0] >= minimum_support_reads 
			and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag) 
			and alelle_sort[-2][2][0] >= 0.3*len(read_tag)):
			breakpointStart = np.mean(alelle_sort[-2][0])
			# breakpointStart_STD = np.std(alelle_sort[-2][0])
			CIPOS = cal_CIPOS(np.std(alelle_sort[-2][0]), len(alelle_sort[-2][0]))
			signalLen = np.mean(alelle_sort[-2][1])
			last_signalLen_STD = signalLen_STD
			signalLen_STD = np.std(alelle_sort[-2][1])
			CILEN = cal_CIPOS(np.std(alelle_sort[-2][1]), len(alelle_sort[-2][1]))
			# search_threshold = np.min(alelle_sort[-2][0])
			if signalLen_STD < last_signalLen_STD:
				# pass
				# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
				'''genotyping'''
				if action:
					# time_start = time.time()
					DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(breakpointStart), chr, alelle_sort[-2][3], 
													max_cluster_bias, gt_round)
					# cost_time = time.time() - time_start
					# print("INS", chr, int(breakpointStart), int(signalLen), DR, DV, QUAL, "%.4f"%cost_time)
				else:
					DR = '.'
					GT = './.'
					GL = '.,.,.'
					GQ = "."
					QUAL = "."
				candidate_single_SV.append([chr, 
											svtype, 
											str(int(breakpointStart)), 
											str(int(signalLen)), 
											str(alelle_sort[-2][2][0]), 
											str(CIPOS),
											str(CILEN),
											str(DR),
											str(GT),
											str(GL),
											str(GQ),
											str(QUAL),
											str(','.join(alelle_sort[-2][3]))])

	elif alelle_sort[-2][2][0] >= minimum_support_reads and alelle_sort[-2][2][0] + alelle_sort[-1][2][0] >= 0.95*len(read_tag):
		if alelle_sort[-2][2][0] >= 0.4*len(read_tag):
			breakpointStart = np.mean(alelle_sort[-1][0])
			# breakpointStart_STD = np.std(alelle_sort[-1][0])
			CIPOS = cal_CIPOS(np.std(alelle_sort[-1][0]), len(alelle_sort[-1][0]))
			signalLen = np.mean(alelle_sort[-1][1])
			signalLen_STD = np.std(alelle_sort[-1][1])
			CILEN = cal_CIPOS(np.std(alelle_sort[-1][1]), len(alelle_sort[-1][1]))
			# search_threshold = np.min(alelle_sort[-1][0])
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-1][2][0], breakpointStart_STD, signalLen_STD)
			'''genotyping'''
			if action:
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(breakpointStart), chr, alelle_sort[-1][3], 
												max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("INS", chr, int(breakpointStart), int(signalLen), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(signalLen)), 
										str(alelle_sort[-1][2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(alelle_sort[-1][3]))])

			breakpointStart = np.mean(alelle_sort[-2][0])
			# breakpointStart_STD = np.std(alelle_sort[-2][0])
			CIPOS = cal_CIPOS(np.std(alelle_sort[-2][0]), len(alelle_sort[-2][0]))
			signalLen = np.mean(alelle_sort[-2][1])
			signalLen_STD = np.std(alelle_sort[-2][1])
			CILEN = cal_CIPOS(np.std(alelle_sort[-2][1]), len(alelle_sort[-2][1]))
			# search_threshold = np.min(alelle_sort[-2][0])
			# print(chr, svtype, int(breakpointStart), int(signalLen), alelle_sort[-2][2][0], breakpointStart_STD, signalLen_STD)
			'''genotyping'''
			if action:
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(breakpointStart), chr, alelle_sort[-2][3], 
												max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("INS", chr, int(breakpointStart), int(signalLen), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpointStart)), 
										str(int(signalLen)), 
										str(alelle_sort[-2][2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(alelle_sort[-2][3]))])

def run_del(args):
	return resolution_DEL(*args)

def run_ins(args):
	return resolution_INS(*args)

def call_gt(bam_path, search_threshold, chr, read_id_list, max_cluster_bias, gt_round):
	import pysam
	querydata = set()
	bamfile = pysam.AlignmentFile(bam_path)
	search_start = max(int(search_threshold) - max_cluster_bias, 0)
	search_end = min(int(search_threshold) + max_cluster_bias, bamfile.get_reference_length(chr))

	up_bound = threshold_ref_count(len(read_id_list))

	status = count_coverage(chr, 
							search_start, 
							search_end, 
							bamfile, 
							querydata, 
							up_bound, 
							gt_round)
	bamfile.close()

	if status == -1:
		DR = '.'
		GT = "./."
		GL = ".,.,."
		GQ = "."
		QUAL = "."

	# elif status == 1:
	# 	pass
	else:
		DR = 0
		for query in querydata:
			if query not in read_id_list:
				DR += 1
		GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))
	return len(read_id_list), DR, GT, GL, GQ, QUAL
