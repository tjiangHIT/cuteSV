import sys
import numpy as np

def resolution_INV(path, chr, svtype, read_count, max_cluster_bias, sv_size):
	'''
	cluster INV
	************************************************************************
	path:	INV.sigs
	chr:	chromosome id
	svtype:	<INV>
	
	SEQTYPE		read_count 	max_cluster_bias 	sv_size 
	------------------------------------------------------------------------
	CCS			5			10 bp (<500 bp)		20 bp 	
	CLR			5			20 bp (<500 bp)		50 bp 	
	------------------------------------------------------------------------
	
	Input file format
	------------------------------------------------------------------------
	column	#1	#2	#3	#4	#5
			INV	CHR	BP1	BP2	ID	
	#1	inversion type
	#2	chromosome number
	#3	breakpoint_1 in each read
	#4	breakpoint_2 in each read
	#5	read ID
	************************************************************************
	'''

	# Initialization of some temporary variables
	semi_inv_cluster = list()
	semi_inv_cluster.append([0,0,''])
	candidate_single_SV = list()

	# Load inputs & cluster breakpoint from each signature read 
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		breakpoint_1_in_read = int(seq[2])
		breakpoint_2_in_read = int(seq[3])
		read_id = seq[4]

		if breakpoint_1_in_read - semi_inv_cluster[-1][0] > max_cluster_bias:
			if len(semi_inv_cluster) >= read_count:
				generate_semi_inv_cluster_2(semi_inv_cluster, 
											chr, 
											svtype, 
											read_count, 
											sv_size, 
											candidate_single_SV, 
											max_cluster_bias)
			semi_inv_cluster = []
			semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id])
		else:
			semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id])

	if len(semi_inv_cluster) >= read_count:
		generate_semi_inv_cluster_2(semi_inv_cluster, 
									chr, svtype, 
									read_count, 
									sv_size, 
									candidate_single_SV, 
									max_cluster_bias)
	file.close()
	return candidate_single_SV

def generate_semi_inv_cluster(semi_inv_cluster, chr, svtype, read_count, sv_size, candidate_single_SV, max_cluster_bias):

	read_id = [i[2] for i in semi_inv_cluster]
	support_read = len(list(set(read_id)))
	if support_read < read_count:
		return

	breakpoint_1_candidate = [i[0] for i in semi_inv_cluster]
	breakpoint_2_candidate_with_read_id = list()
	for i in semi_inv_cluster:
		breakpoint_2_candidate_with_read_id.append([i[1], i[2]])
	breakpoint_2_candidate = sorted(list(breakpoint_2_candidate_with_read_id), key = lambda x:x[0])

	# print("#")
	# for i in breakpoint_2_candidate:
	# 	print(i)

	breakpoint_1 = np.mean(breakpoint_1_candidate)
	last_bp = breakpoint_2_candidate[0][0]
	temp_count = 1
	max_count = 0
	max_count_id = 0
	temp_sum = last_bp
	max_sum = 0
	temp_id = dict()
	temp_id[breakpoint_2_candidate[0][1]] = 0

	for i in breakpoint_2_candidate[1:]:
		if i[0] - last_bp > max_cluster_bias:
			if temp_count > max_count:
				max_count = temp_count
				max_sum = temp_sum
				max_count_id = len(temp_id)
			temp_id = dict()
			temp_count = 1
			temp_sum = i[0]
			temp_id[i[1]] = 0
		else:
			if i[1] not in temp_id:
				temp_id[i[1]] = 0
			else:
				temp_id[i[1]] += 1
			temp_count += 1
			temp_sum += i[0]
		last_bp = i[0]
	if temp_count > max_count:
		max_count = temp_count
		max_sum = temp_sum
		max_count_id = len(temp_id)

	# print(max_count, len(temp_id), support_read)
	breakpoint_2 = int(max_sum / max_count)
	inv_len = breakpoint_2 - breakpoint_1

	if inv_len >= sv_size and max_count_id >= read_count:
		candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))


def generate_semi_inv_cluster_2(semi_inv_cluster, chr, svtype, read_count, sv_size, candidate_single_SV, max_cluster_bias):

	read_id = [i[2] for i in semi_inv_cluster]
	support_read = len(list(set(read_id)))
	if support_read < read_count:
		return

	inv_cluster_b2 = sorted(semi_inv_cluster, key = lambda x:x[1])


	# breakpoint_1 = np.mean(breakpoint_1_candidate)
	last_bp = inv_cluster_b2[0][1]
	temp_count = 1
	# max_count = 0
	temp_sum_b1 = inv_cluster_b2[0][0]
	temp_sum_b2 = last_bp

	# max_sum = 0
	temp_id = dict()
	temp_id[inv_cluster_b2[0][2]] = 0

	for i in inv_cluster_b2[1:]:
		if i[1] - last_bp > max_cluster_bias:
			if temp_count >= read_count:
				max_count_id = len(temp_id)

				breakpoint_1 = round(temp_sum_b1 / temp_count)
				breakpoint_2 = round(temp_sum_b2 / temp_count)
				inv_len = breakpoint_2 - breakpoint_1
				if inv_len >= sv_size and max_count_id >= read_count:
					# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))
					candidate_single_SV.append([chr, 
												svtype, 
												str(int(breakpoint_1)), 
												str(int(inv_len)), 
												str(max_count_id)])

			temp_id = dict()
			temp_count = 1
			temp_sum_b1 = i[0]
			temp_sum_b2 = i[1]
			temp_id[i[2]] = 0
		else:
			if i[2] not in temp_id:
				temp_id[i[2]] = 0
			else:
				temp_id[i[2]] += 1
			temp_count += 1
			temp_sum_b1 += i[0]
			temp_sum_b2 += i[1]
		last_bp = i[1]
	if temp_count >= read_count:
		max_count_id = len(temp_id)
		breakpoint_1 = round(temp_sum_b1 / temp_count)
		breakpoint_2 = round(temp_sum_b2 / temp_count)
		inv_len = breakpoint_2 - breakpoint_1
		if inv_len >= sv_size and max_count_id >= read_count:
			# candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\n'%(chr, svtype, breakpoint_1, breakpoint_2, max_count_id))
			candidate_single_SV.append([chr, 
										svtype, 
										str(int(breakpoint_1)), 
										str(int(inv_len)), 
										str(max_count_id)])

def run_inv(args):
	return resolution_INV(*args)