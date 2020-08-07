import sys
import numpy as np
from cuteSV.cuteSV_genotype import cal_GL, threshold_ref_count, count_coverage

'''
*******************************************
				TO DO LIST
*******************************************
	1. Adjust format in BND format;
	2. Make TRA calls more sensitivity;
	3. Collect DP info.
	4. Determine (IM)PRECISE type.

	NOTICE:
		For INS larger than read,
		named it with UNRESOLVED.
*******************************************
'''

'''
			*********Description*********
			*	TYPE A:		N[chr:pos[	*
			*	TYPE B:		N]chr:pos]	*
			*	TYPE C:		[chr:pos[N	*
			*	TYPE D:		]chr:pos]N	*
			*****************************
			'''

def resolution_TRA(path, chr_1, chr_2, read_count, overlap_size, max_cluster_bias, bam_path, action, gt_round):
	semi_tra_cluster = list()
	semi_tra_cluster.append([0,0,'','N'])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr_1:
			continue
		if seq[4] != chr_2:
			continue

		pos_1 = int(seq[3])
		pos_2 = int(seq[5])
		read_id = seq[6]
		BND_type = seq[2]
	
		if pos_1 - semi_tra_cluster[-1][0] > max_cluster_bias or BND_type != semi_tra_cluster[-1][3]:
			if len(semi_tra_cluster) >= read_count:
				if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
					pass
				else:
					generate_semi_tra_cluster(semi_tra_cluster, 
											chr_1, 
											chr_2, 
											read_count, 
											overlap_size, 
											max_cluster_bias, 
											candidate_single_SV,
											bam_path,
											action,
											gt_round)
			semi_tra_cluster = []
			semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type])
		else:
			semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type])

	if len(semi_tra_cluster) >= read_count:
		if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
			pass
		else:
			generate_semi_tra_cluster(semi_tra_cluster, 
									chr_1, 
									chr_2, 
									read_count, 
									overlap_size, 
									max_cluster_bias, 
									candidate_single_SV,
									bam_path,
									action,
									gt_round)
	file.close()
	return candidate_single_SV

def generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, 
	max_cluster_bias, candidate_single_SV, bam_path, action, gt_round):
	BND_type = semi_tra_cluster[0][3]
	semi_tra_cluster = sorted(semi_tra_cluster, key = lambda x:x[1])
	read_tag = dict()
	temp = list()
	# p1, p2, count
	last_len = 0
	temp.append([0,0,set()])
	for element in semi_tra_cluster:
		if element[1] - last_len > max_cluster_bias:
			temp.append([element[0],element[1],{element[2]}])
			last_len = element[1]
		else:
			temp[-1][0] += element[0]
			temp[-1][1] += element[1]
			temp[-1][2].add(element[2])
			last_len = element[1]

		if element[2] not in read_tag:
			read_tag[element[2]] = 0
	if len(read_tag) < read_count:
		return

	temp = sorted(temp, key = lambda x:-len(x[2]))

	if len(temp[1][2]) >= 0.5*read_count:
		if len(temp[0][2])+len(temp[1][2]) >= len(semi_tra_cluster)*overlap_size:
			# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[1][0]/temp[1][2]), chr_2, int(temp[1][1]/temp[1][2]), len(read_tag)))
			BND_pos_1 = "%s:%s"%(chr_2, int(temp[0][1]/len(temp[0][2])))
			BND_pos_2 = "%s:%s"%(chr_2, int(temp[1][1]/len(temp[1][2])))
			if BND_type == 'A':
				TRA_1 = "N[%s["%(BND_pos_1)
				TRA_2 = "N[%s["%(BND_pos_2)
			elif BND_type == 'B':
				TRA_1 = "N]%s]"%(BND_pos_1)
				TRA_2 = "N]%s]"%(BND_pos_2)
			elif BND_type == 'C':
				TRA_1 = "[%s[N"%(BND_pos_1)
				TRA_2 = "[%s[N"%(BND_pos_2)
			elif BND_type == 'D':
				TRA_1 = "]%s]N"%(BND_pos_1)
				TRA_2 = "]%s]N"%(BND_pos_2)
			else:
				return

			if action:
				import time
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(temp[0][0]/len(temp[0][2])), 
											int(temp[0][1]/len(temp[0][2])), chr_1, chr_2, temp[0][2], 
											max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("BND", chr_1, chr_2, int(temp[0][0]/len(temp[0][2])), int(temp[0][1]/len(temp[0][2])), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr_1, 
										TRA_1, 
										str(int(temp[0][0]/len(temp[0][2]))), 
										chr_2, 
										str(int(temp[0][1]/len(temp[0][2]))), 
										str(len(temp[0][2])),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(temp[0][2]))])

			if action:
				import time
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(temp[1][0]/len(temp[1][2])), 
											int(temp[1][1]/len(temp[1][2])), chr_1, chr_2, temp[1][2], 
											max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("BND", chr_1, chr_2, int(temp[1][0]/len(temp[1][2])), int(temp[1][1]/len(temp[1][2])), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr_1, 
										TRA_2, 
										str(int(temp[1][0]/len(temp[1][2]))), 
										chr_2, 
										str(int(temp[1][1]/len(temp[1][2]))), 
										str(len(temp[1][2])),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(temp[1][2]))])
	else:
		if len(temp[0][2]) >= len(semi_tra_cluster)*overlap_size:
			# print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			BND_pos = "%s:%s"%(chr_2, int(temp[0][1]/len(temp[0][2])))
			if BND_type == 'A':
				TRA = "N[%s["%(BND_pos)
			elif BND_type == 'B':
				TRA = "N]%s]"%(BND_pos)
			elif BND_type == 'C':
				TRA = "[%s[N"%(BND_pos)
			elif BND_type == 'D':
				TRA = "]%s]N"%(BND_pos)
			else:
				return

			if action:
				import time
				# time_start = time.time()
				DV, DR, GT, GL, GQ, QUAL = call_gt(bam_path, int(temp[0][0]/len(temp[0][2])), 
											int(temp[0][1]/len(temp[0][2])), chr_1, chr_2, temp[0][2], 
											max_cluster_bias, gt_round)
				# cost_time = time.time() - time_start
				# print("BND", chr_1, chr_2, int(temp[0][0]/len(temp[0][2])), int(temp[0][1]/len(temp[0][2])), DR, DV, QUAL, "%.4f"%cost_time)
			else:
				DR = '.'
				GT = './.'
				GL = '.,.,.'
				GQ = "."
				QUAL = "."
			candidate_single_SV.append([chr_1, 
										TRA, 
										str(int(temp[0][0]/len(temp[0][2]))), 
										chr_2, 
										str(int(temp[0][1]/len(temp[0][2]))), 
										str(len(temp[0][2])),
										str(DR),
										str(GT),
										str(GL),
										str(GQ),
										str(QUAL),
										str(','.join(temp[0][2]))])


def run_tra(args):
	return resolution_TRA(*args)

def call_gt(bam_path, pos_1, pos_2, chr_1, chr_2, read_id_list, max_cluster_bias, gt_round):
	import pysam
	bamfile = pysam.AlignmentFile(bam_path)
	querydata = set()
	search_start = max(int(pos_1) - max_cluster_bias, 0)
	search_end = min(int(pos_1) + max_cluster_bias, bamfile.get_reference_length(chr_1))

	up_bound = threshold_ref_count(len(read_id_list))

	status = count_coverage(chr_1, 
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
		search_start = max(int(pos_2) - max_cluster_bias, 0)
		search_end = min(int(pos_2) + max_cluster_bias, bamfile.get_reference_length(chr_2))
		status_2 = count_coverage(chr_2, 
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
