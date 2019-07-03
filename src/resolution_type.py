import sys
import numpy as np
# from scipy.stats import kurtosis

def run_indel_inv(args):
	return resolution_INDEL(*args)

def resolution_INDEL(path, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, max_distance):
	'''
	cluster INDEL
	************************************************************************
	path:	INDEL.sigs
	chr:	chromosome id
	svtype:	<INDEL>
	
	read_count 	overlap_size 	max_cluster_bias 	sv_size 	max_distance
	------------------------------------------------------------------------
	5/10		0.5				50 bp (<500 bp)		50 bp 		1000 bp
	************************************************************************
	'''
	semi_indel_cluster = list()
	semi_indel_cluster.append([0,0,''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos = int(seq[2])
		indel_len = int(seq[3])
		read_id = seq[4]
		
		if pos - semi_indel_cluster[-1][0] > max_cluster_bias:
			if len(semi_indel_cluster) >= read_count:
				generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, candidate_single_SV)
			semi_indel_cluster = []
			semi_indel_cluster.append([pos, indel_len, read_id])
		else:
			semi_indel_cluster.append([pos, indel_len, read_id])

	if len(semi_indel_cluster) >= read_count:
		generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, candidate_single_SV)
	file.close()
	# return polish_indel(candidate_single_SV, max_distance, read_count)
	return candidate_single_SV

def polish_indel(candidate_single_SV, max_distance, read_count):
	polish_indel_candidate = list()
	if len(candidate_single_SV) <= 1:
		# return candidate_single_SV
		for temp in candidate_single_SV:
			encode_temp = "%s\t%s\t%d\t%d\t%d\n"%(temp[0], temp[1], temp[2], temp[3], temp[4])
			if temp[4] >= read_count:
				polish_indel_candidate.append(encode_temp)
		return polish_indel_candidate
	
	temp = candidate_single_SV[0]
	for i in candidate_single_SV[1:]:
		if temp[2] <= i[2] and i[2] <= temp[2]+temp[3] and i[2]-temp[2] <= max_distance:
			# if temp[4] < i[4]:
			# 	temp = i
			new_read_count = temp[4] + i[4]
			new_break_point = int((temp[2]*temp[4]+i[2]*i[4])/new_read_count)
			new_indel_length = int((temp[3]*temp[4]+i[3]*i[4])/new_read_count)
			temp = [i[0], i[1], new_break_point, new_indel_length, new_read_count]
		else:
			encode_temp = "%s\t%s\t%d\t%d\t%d\n"%(temp[0], temp[1], temp[2], temp[3], temp[4])
			# polish_indel_candidate.append(encode_temp)
			if temp[4] >= read_count:
				polish_indel_candidate.append(encode_temp)
			temp = i
	encode_temp = "%s\t%s\t%d\t%d\t%d\n"%(temp[0], temp[1], temp[2], temp[3], temp[4])
	# polish_indel_candidate.append(encode_temp)
	if temp[4] >= read_count:
		polish_indel_candidate.append(encode_temp)
	return polish_indel_candidate

def generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, candidate_single_SV):
	'''
	generate INDEL
	*******************************************
	overlap_size 	max_cluster_bias 	sv_size
	-------------------------------------------
	0.5				50 bp (<500 bp)		50 bp
	*******************************************
	'''
	# unique read id
	read_tag = dict()
	for element in semi_indel_cluster:
		if element[2] not in read_tag:
			read_tag[element[2]] = element
		else:
			if element[1] > read_tag[element[2]][1]:
				read_tag[element[2]] = element

	if len(read_tag) < read_count:
		return

	read_tag2SortedList = sorted(list(read_tag.values()), key = lambda x:x[1])

	last_len = -DISCRETE_THRESHOLD_LEN_CLUSTER
	max_conut = 0
	max_LEN_sum = 0
	max_bps_sum = 0
	temp_count = 0
	temp_LEN_sum = 0
	temp_bps_sum = 0
	for i in read_tag2SortedList:
		if i[1] - last_len > DISCRETE_THRESHOLD_LEN_CLUSTER:
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

	breakpoint_starts = int(max_bps_sum / max_conut)
	signal_len = int(max_LEN_sum / max_conut)

	candidate_single_SV.append('%s\t%s\t%d\t%d\t%d\t%d\t%.3f\n'%
		(chr, svtype, breakpoint_starts, signal_len, len(read_tag), max_conut, float(max_conut * 1.0 / len(read_tag))))

def polish_dup(candidate_single_SV, max_distance):
	polish_dup_candidate = list()
	if len(candidate_single_SV) <= 1:
		# return candidate_single_SV
		for temp in candidate_single_SV:
			encode_temp = "%s\t%s\t%d\t%d\t%d\n"%(temp[0], temp[1], temp[2], temp[3], temp[4])
			polish_dup_candidate.append(encode_temp)
		return polish_dup_candidate

	temp = candidate_single_SV[0]
	for i in candidate_single_SV[1:]:
		if temp[2] <= i[2] and i[2] <= temp[2]+temp[3] and i[2]-temp[2] <= max_distance:
			if temp[4] < i[4]:
				temp = i
		else:
			encode_temp = "%s\t%s\t%d\t%d\t%d\n"%(temp[0], temp[1], temp[2], temp[3], temp[4])
			polish_dup_candidate.append(encode_temp)
			temp = i
	encode_temp = "%s\t%s\t%d\t%d\t%d\n"%(temp[0], temp[1], temp[2], temp[3], temp[4])
	polish_dup_candidate.append(encode_temp)
	return polish_dup_candidate





#if __name__ == '__main__':
	# candidate_single_SV = []
	# resolution_INDEL(sys.argv[1], "1", "DEL", 10, 0.5, 20, 50)
	# file_out

	# resolution_TRA(sys.argv[1], "1", "hs37d5", 10, 0.6, 50)
	# resolution_DUP(sys.argv[1], "1", 10, 50, 50)
