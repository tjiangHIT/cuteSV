import sys
import numpy as np

def resolution_TRA(path, chr_1, chr_2, read_count, overlap_size, max_cluster_bias):
	semi_tra_cluster = list()
	semi_tra_cluster.append([0,0,''])
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr_1:
			continue
		if seq[3] != chr_2:
			continue

		pos_1 = int(seq[2])
		pos_2 = int(seq[4])
		read_id = seq[5]
		
		if pos_1 - semi_tra_cluster[-1][0] > max_cluster_bias:
			if len(semi_tra_cluster) >= read_count:
				generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias, candidate_single_SV)
			semi_tra_cluster = []
			semi_tra_cluster.append([pos_1, pos_2, read_id])
		else:
			semi_tra_cluster.append([pos_1, pos_2, read_id])

	if len(semi_tra_cluster) >= read_count:
		generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias, candidate_single_SV)
	file.close()
	return candidate_single_SV

def generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias, candidate_single_SV):
	semi_tra_cluster = sorted(semi_tra_cluster, key = lambda x:x[1])
	read_tag = dict()
	temp = list()
	# p1, p2, count
	last_len = 0
	temp.append([0,0,0])
	for element in semi_tra_cluster:
		if element[1] - last_len > max_cluster_bias:
			temp.append([element[0],element[1],1])
			last_len = element[1]
		else:
			temp[-1][0] += element[0]
			temp[-1][1] += element[1]
			temp[-1][2] += 1
			last_len = element[1]

		if element[2] not in read_tag:
			read_tag[element[2]] = 0
	if len(read_tag) < read_count:
		return

	temp = sorted(temp, key = lambda x:-x[2])

	if temp[1][2] >= 0.5*read_count:
		if temp[0][2]+temp[1][2] >= len(semi_tra_cluster)*overlap_size:
			# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[1][0]/temp[1][2]), chr_2, int(temp[1][1]/temp[1][2]), len(read_tag)))
			candidate_single_SV.append([chr_1, "TRA", str(int(temp[0][0]/temp[0][2])), chr_2, str(int(temp[0][1]/temp[0][2])), str(len(read_tag))])
			candidate_single_SV.append([chr_1, "TRA", str(int(temp[1][0]/temp[1][2])), chr_2, str(int(temp[1][1]/temp[1][2])), str(len(read_tag))])
	else:
		if temp[0][2] >= len(semi_tra_cluster)*overlap_size:
			# print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			candidate_single_SV.append([chr_1, "TRA", str(int(temp[0][0]/temp[0][2])), chr_2, str(int(temp[0][1]/temp[0][2])), str(len(read_tag))])


def run_tra(args):
	return resolution_TRA(*args)
