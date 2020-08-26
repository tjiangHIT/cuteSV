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
				 minimum_support_reads, bam_path, action, gt_round):

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
										# threshold_local, 
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
								# threshold_local, 
								minimum_support_reads, 
								candidate_single_SV,
								bam_path,
								max_cluster_bias,
								action,
								gt_round)
	file.close()
	return candidate_single_SV

def generate_del_cluster(semi_del_cluster, chr, svtype, read_count, 
	threshold_gloab, minimum_support_reads, candidate_single_SV, 
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

	allele_collect = list()
	'''
	*************************************************************
		#1 				#2			#3			#4
	-------------------------------------------------------------
		del-breakpoint	del-len		#support 	read-id
	*************************************************************
	'''
	allele_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[],
		[read_tag2SortedList[0][2]]])

	for i in read_tag2SortedList[1:]:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP:
			allele_collect[-1][2].append(len(allele_collect[-1][0]))
			allele_collect.append([[],[],[],[]])

		allele_collect[-1][0].append(i[0])
		allele_collect[-1][1].append(i[1])
		allele_collect[-1][3].append(i[2])
		last_len = i[1]
	allele_collect[-1][2].append(len(allele_collect[-1][0]))
	allele_sort = sorted(allele_collect, key = lambda x:x[2])

	for allele in allele_sort:
		if allele[2][0] >= minimum_support_reads:
			breakpointStart = np.mean(allele[0])
			search_threshold = np.min(allele[0])
			CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
			signalLen = np.mean(allele[1])
			signalLen_STD = np.std(allele[1])
			CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))

			if action:
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, 
													int(search_threshold), 
													chr, 
													allele[3], 
													max_cluster_bias, 
													gt_round)
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
										str(allele[2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(allele[3]))])
	

def resolution_INS(path, chr, svtype, read_count, threshold_gloab, 
	max_cluster_bias, minimum_support_reads, bam_path, action, gt_round):
	
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
	semi_ins_cluster.append([0,0,'',''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos = int(seq[2])
		indel_len = int(seq[3])
		read_id = seq[4]
		ins_seq = seq[5]
		
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
										# threshold_local, 
										minimum_support_reads, 
										candidate_single_SV,
										bam_path,
										max_cluster_bias,
										action,
										gt_round)
			semi_ins_cluster = []
			semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])
		else:
			semi_ins_cluster.append([pos, indel_len, read_id, ins_seq])

	if len(semi_ins_cluster) >= read_count:
		if semi_ins_cluster[-1][0] == semi_ins_cluster[-1][1] == 0:
			pass
		else:
			generate_ins_cluster(semi_ins_cluster, 
								chr, 
								svtype, 
								read_count, 
								threshold_gloab, 
								# threshold_local, 
								minimum_support_reads, 
								candidate_single_SV,
								bam_path,
								max_cluster_bias,
								action,
								gt_round)
	file.close()
	return candidate_single_SV

def generate_ins_cluster(semi_ins_cluster, chr, svtype, read_count, 
	threshold_gloab, minimum_support_reads, candidate_single_SV, 
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

	allele_collect = list()
	allele_collect.append([[read_tag2SortedList[0][0]],
							[read_tag2SortedList[0][1]],
							[], 
							[read_tag2SortedList[0][2]],
							[read_tag2SortedList[0][3]]])

	for i in read_tag2SortedList[1:]:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER_INS_TEMP:
			allele_collect[-1][2].append(len(allele_collect[-1][0]))
			allele_collect.append([[],[],[],[],[]])

		allele_collect[-1][0].append(i[0])
		allele_collect[-1][1].append(i[1])
		allele_collect[-1][3].append(i[2])
		allele_collect[-1][4].append(i[3])
		last_len = i[1]
	allele_collect[-1][2].append(len(allele_collect[-1][0]))
	allele_sort = sorted(allele_collect, key = lambda x:x[2])

	for allele in allele_sort:
		if allele[2][0] >= minimum_support_reads:
			breakpointStart = np.mean(allele[0])
			CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
			signalLen = np.mean(allele[1])
			signalLen_STD = np.std(allele[1])
			CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))
			ideal_ins_seq = '<INS>'
			for i in allele[4]:
				if len(i) >= int(signalLen):
					ideal_ins_seq = i[0:int(signalLen)]
					break
			if ideal_ins_seq == '<INS>':
				continue
				
			if action:
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, 
													int(breakpointStart), 
													chr, 
													allele[3], 
													# max_cluster_bias, 
													1000,
													gt_round)
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
										str(allele[2][0]), 
										str(CIPOS),
										str(CILEN),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(allele[3])),
										ideal_ins_seq])


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
