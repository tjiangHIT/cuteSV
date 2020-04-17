import sys
import numpy as np

def resolution_INV(path, chr, svtype, read_count, max_cluster_bias, sv_size, 
	bam_path, action, hom, het, MaxSize):
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
				generate_semi_inv_cluster(semi_inv_cluster, 
											chr, 
											svtype, 
											read_count, 
											sv_size, 
											candidate_single_SV, 
											max_cluster_bias,
											bam_path,
											action,
											hom,
											het,
											MaxSize)
			semi_inv_cluster = []
			semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id])
		else:
			semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id])

	if len(semi_inv_cluster) >= read_count:
		generate_semi_inv_cluster(semi_inv_cluster, 
									chr, svtype, 
									read_count, 
									sv_size, 
									candidate_single_SV, 
									max_cluster_bias,
									bam_path,
									action,
									hom,
									het,
									MaxSize)
	file.close()
	return candidate_single_SV

def generate_semi_inv_cluster(semi_inv_cluster, chr, svtype, read_count, sv_size, 
	candidate_single_SV, max_cluster_bias, bam_path, action, hom, het, MaxSize):

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
					if inv_len <= MaxSize:
						if action:
							DV, DR, GT = call_gt(bam_path, int(breakpoint_1), 
								int(breakpoint_2), chr, list(temp_id.keys()), 
								max_cluster_bias, hom, het)
						else:
							DR = '.'
							GT = './.'
						candidate_single_SV.append([chr, 
													svtype, 
													str(int(breakpoint_1)), 
													str(int(inv_len)), 
													str(max_count_id),
													str(DR),
													str(GT)])

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
			if inv_len <= MaxSize:
				if action:
					DV, DR, GT = call_gt(bam_path, int(breakpoint_1), 
						int(breakpoint_2), chr, list(temp_id.keys()), 
						max_cluster_bias, hom, het)
				else:
					DR = '.'
					GT = './.'
				candidate_single_SV.append([chr, 
											svtype, 
											str(int(breakpoint_1)), 
											str(int(inv_len)), 
											str(max_count_id),
											str(DR),
											str(GT)])

def run_inv(args):
	return resolution_INV(*args)

def count_coverage(chr, s, e, f):
	read_count = set()
	for i in f.fetch(chr, s, e):
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
	