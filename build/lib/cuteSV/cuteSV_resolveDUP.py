import sys
import numpy as np
from collections import Counter
from cuteSV.cuteSV_genotype import cal_GL, threshold_ref_count, count_coverage

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
	bam_path, action, MaxSize, gt_round):
	semi_dup_cluster = list()
	semi_dup_cluster.append([0, 0, ''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos_1 = int(seq[2])
		pos_2 = int(seq[3])
		read_id = seq[4]
		
		if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias or pos_2 - semi_dup_cluster[-1][1] > max_cluster_bias:
			if len(semi_dup_cluster) >= read_count:
				if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
					pass
				else:
					# if chr == "1":
					# 	print("DUP")
					# 	print(int(sum([i[0] for i in semi_dup_cluster])/len(semi_dup_cluster)), 
					# 		int(sum([i[1] for i in semi_dup_cluster])/len(semi_dup_cluster)))
					generate_dup_cluster(semi_dup_cluster, 
										chr, 
										read_count, 
										max_cluster_bias, 
										sv_size, 
										candidate_single_SV,
										bam_path,
										action,
										MaxSize,
										gt_round)
			semi_dup_cluster = []
			semi_dup_cluster.append([pos_1, pos_2, read_id])
		else:
			semi_dup_cluster.append([pos_1, pos_2, read_id])

	if len(semi_dup_cluster) >= read_count:
		if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
			pass
		else:
			# if chr == "1":
			# 	print("DUP")	
			# 	print(int(sum([i[0] for i in semi_dup_cluster])/len(semi_dup_cluster)), 
			# 		int(sum([i[1] for i in semi_dup_cluster])/len(semi_dup_cluster)))
			generate_dup_cluster(semi_dup_cluster, 
								chr, 
								read_count, 
								max_cluster_bias, 
								sv_size, 
								candidate_single_SV,
								bam_path,
								action,
								MaxSize,
								gt_round)
	file.close()
	return candidate_single_SV

def generate_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, 
	sv_size, candidate_single_SV, bam_path, action, MaxSize, gt_round):
	# calculate support reads
	support_read = list(set([i[2] for i in semi_dup_cluster]))
	# print(support_read)
	if len(support_read) < read_count:
		return

	low_b = int(len(semi_dup_cluster)*0.4)
	up_b = int(len(semi_dup_cluster)*0.6)
	# print(low_b, up_b)

	if low_b == up_b:
		breakpoint_1 = semi_dup_cluster[low_b][0]
		breakpoint_2 = semi_dup_cluster[low_b][1]
	else:
		breakpoint_1 = [i[0] for i in semi_dup_cluster[low_b:up_b]]
		breakpoint_2 = [i[1] for i in semi_dup_cluster[low_b:up_b]]
		breakpoint_1 = int(sum(breakpoint_1)/len(semi_dup_cluster[low_b:up_b]))
		breakpoint_2 = int(sum(breakpoint_2)/len(semi_dup_cluster[low_b:up_b]))

	# print(chr, breakpoint_1, breakpoint_2)

	if sv_size <= breakpoint_2 - breakpoint_1 <= MaxSize:
		if action:
			import time
			# time_start = time.time()
			DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, 
												breakpoint_1, 
												breakpoint_2, 
												chr, 
												support_read, 
												min(max_cluster_bias, breakpoint_2 - breakpoint_1),
												gt_round)
			# print(DV, DR, GT, GL, GQ, QUAL)
			# cost_time = time.time() - time_start
			# print("DUP", chr, int(breakpoint_1), int(breakpoint_2), DR, DV, QUAL, "%.4f"%cost_time)
		else:
			DR = '.'
			GT = './.'
			GL = '.,.,.'
			GQ = "."
			QUAL = "."
		candidate_single_SV.append([chr,
									'DUP', 
									str(breakpoint_1), 
									str(breakpoint_2 - breakpoint_1), 
									str(len(support_read)),
									str(DR),
									str(GT),
									str(GL),
									str(GQ),
									str(QUAL),
									str(','.join(support_read))])


# def generate_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, 
# 	sv_size, candidate_single_SV, bam_path, action, MaxSize):
# 	# calculate support reads
# 	support_read = list(set([i[2] for i in semi_dup_cluster]))
# 	# print(support_read)
# 	if len(support_read) < read_count:
# 		return

# 	dup_cluster_len = sorted(semi_dup_cluster, key = lambda x:x[1])
# 	# for i in dup_cluster_len:
# 	# 	print(i)

# 	last_p2 = dup_cluster_len[0][1]
# 	temp_list = dict()
# 	temp_list[dup_cluster_len[0][2]] = 0
# 	temp_count = 1
# 	subcluster_bp1 = list()
# 	subcluster_bp1.append(dup_cluster_len[0][0])
# 	# temp_LEN_sum = last_p2
# 	subcluster_bp2 = list()
# 	subcluster_bp2.append(last_p2)

# 	predictions = list()
# 	for i in dup_cluster_len[1:]:
# 		if i[1] - last_p2 > max_cluster_bias:
# 			if temp_count >= read_count:
# 				if len(temp_list) >= read_count:
# 					# breakpoint_1 = Counter(subcluster_bp1).most_common(1)[0]
# 					# breakpoint_2 = Counter(subcluster_bp2).most_common(1)[0]
# 					breakpoint_1 = (int(sum(subcluster_bp1[int(len(subcluster_bp1)*0.4):int(len(subcluster_bp1)*0.6)])
# 									 / len(subcluster_bp1[int(len(subcluster_bp1)*0.4):int(len(subcluster_bp1)*0.6)])))
# 					breakpoint_2 = (int(sum(subcluster_bp2[int(len(subcluster_bp2)*0.4):int(len(subcluster_bp2)*0.6)])
# 									 / len(subcluster_bp2[int(len(subcluster_bp2)*0.4):int(len(subcluster_bp2)*0.6)])))
# 					predictions.append([breakpoint_1, breakpoint_2, len(temp_list)])
# 					print(chr, breakpoint_1, breakpoint_2, len(temp_list))

# 			temp_list = dict()
# 			temp_count = 1
# 			subcluster_bp1 = list()
# 			subcluster_bp1.append(i[0])
# 			subcluster_bp2 = list()
# 			subcluster_bp2.append(i[1])
# 			temp_list[i[2]] = 0
# 		else:
# 			if i[2] not in temp_list:
# 				temp_list[i[2]] = 0
# 			else:
# 				temp_list[i[2]] += 1
# 			temp_count += 1
# 			subcluster_bp1.append(i[0])
# 			subcluster_bp2.append(i[1])
# 		last_p2 = i[1]

# 	if temp_count >= read_count:
# 		if len(temp_list) >= read_count:
# 			# breakpoint_1 = Counter(subcluster_bp1).most_common(1)[0]
# 			# breakpoint_2 = Counter(subcluster_bp2).most_common(1)[0]
# 			breakpoint_1 = (int(sum(subcluster_bp1[int(len(subcluster_bp1)*0.4):int(len(subcluster_bp1)*0.6)])
# 							 / len(subcluster_bp1[int(len(subcluster_bp1)*0.4):int(len(subcluster_bp1)*0.6)])))
# 			breakpoint_2 = (int(sum(subcluster_bp2[int(len(subcluster_bp2)*0.4):int(len(subcluster_bp2)*0.6)])
# 							 / len(subcluster_bp2[int(len(subcluster_bp2)*0.4):int(len(subcluster_bp2)*0.6)])))
# 			predictions.append([breakpoint_1, breakpoint_2, len(temp_list)])
# 			print(chr, breakpoint_1, breakpoint_2, len(temp_list))

# 	predictions = sorted(predictions, key = lambda x:-x[2])

# 	for i in predictions:
# 		# print(i)
# 		# if i[2] >= read_count and i[2] >= int(len(support_read) / 3) and i[1] - i[0] >= sv_size:
# 		if i[2] >= read_count and i[1] - i[0] >= sv_size:
# 			# if min(i[0][1], i[1][1]) * 2 >= i[2]:
# 			# 	reliability = "PRECISION"
# 			# else:
# 			# 	reliability = "IMPRECISION"

# 			# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%s\n' % (chr, 'DUP', i[0][0], i[1][0] - i[0][0], i[2], reliability))
# 			'''genotyping'''
# 			if i[1] - i[0] <= MaxSize:
# 				if action:
# 					DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(i[0]), int(i[1]), 
# 												chr, support_read, min(max_cluster_bias, i[1] - i[0]))
# 				else:
# 					DR = '.'
# 					GT = './.'
# 					GL = '.,.,.'
# 					GQ = "."
# 					QUAL = "."
# 				candidate_single_SV.append([chr,
# 											'DUP', 
# 											str(i[0]), 
# 											str(i[1] - i[0]), 
# 											str(i[2]),
# 											str(DR),
# 											str(GT),
# 											str(GL),
# 											str(GQ),
# 											str(QUAL)])

def run_dup(args):
	return resolution_DUP(*args)

def call_gt(bam_path, pos_1, pos_2, chr, read_id_list, max_cluster_bias, gt_round):
	import pysam
	bamfile = pysam.AlignmentFile(bam_path)
	querydata = set()
	search_start = max(int(pos_1 - max_cluster_bias/2), 0)
	search_end = min(int(pos_1 + max_cluster_bias/2), bamfile.get_reference_length(chr))

	up_bound = threshold_ref_count(len(read_id_list))
	status = count_coverage(chr, 
							search_start, 
							search_end, 
							bamfile, 
							querydata, 
							up_bound, 
							gt_round)

	if status == -1:
		DR = '.'
		GT = "./."
		GL = ".,.,."
		GQ = "."
		QUAL = "."

	elif status == 1:
		DR = 0
		for query in querydata:
			if query not in read_id_list:
				DR += 1
		GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))

	else:
		search_start = max(int(pos_2 - max_cluster_bias/2), 0)
		search_end = min(int(pos_2 + max_cluster_bias/2), bamfile.get_reference_length(chr))
		status_2 = count_coverage(chr, 
									search_start, 
									search_end, 
									bamfile, 
									querydata, 
									up_bound, 
									gt_round)
		# status_2 judgement
		DR = 0
		for query in querydata:
			if query not in read_id_list:
				DR += 1
		GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))

	bamfile.close()
	return len(read_id_list), DR, GT, GL, GQ, QUAL