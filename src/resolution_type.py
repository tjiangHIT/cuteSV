import sys

candidate_single_SV = list()

def resolution_INDEL(path, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size=50):
	semi_indel_cluster = list()
	semi_indel_cluster.append([0,0,''])

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		pos = int(seq[2])
		indel_len = int(seq[3])
		read_id = seq[4]
		
		if pos - semi_indel_cluster[-1][0] > max_cluster_bias:
			if len(semi_indel_cluster) > read_count/2:
				generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size)
			semi_indel_cluster = []
			semi_indel_cluster.append([pos, indel_len, read_id])
		else:
			semi_indel_cluster.append([pos, indel_len, read_id])

	if len(semi_indel_cluster) > read_count/2:
		generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size)

def generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size):
	read_tag = dict()
	phase_indel_len = list()
	breakpoint_ele = 0
	for element in semi_indel_cluster:
		phase_indel_len.append(element[1])
		breakpoint_ele += element[0]
		if element[2] not in read_tag:
			read_tag[element[2]] = 0

	if len(read_tag) < read_count/2:
		return
	phase_indel_len.sort()
	last_len = -max_cluster_bias
	max_conut = 0
	max_conut_sum = 0
	temp_count = 0
	temp_count_sum = 0
	for i in phase_indel_len:
		if i - last_len > max_cluster_bias:
			if temp_count > max_conut:
				max_conut = temp_count
				max_conut_sum = temp_count_sum
			temp_count = 1
			temp_count_sum = i
		else:
			temp_count += 1
			temp_count_sum += i
		last_len = i
	if temp_count > max_conut:
		max_conut = temp_count
		max_conut_sum = temp_count_sum

	# result
	# chr, "DEL", breakpoint, len, read_count
	breakpoint = int(breakpoint_ele / len(semi_indel_cluster))
	indel_len = int(max_conut_sum/max_conut)
	if max_conut < int(len(read_tag)*overlap_size):
		return
	if svtype == "INV" and indel_len-breakpoint >= sv_size:
		print("%s\t%s\t%d\t%d\t%d"%(chr, svtype, breakpoint, indel_len, len(read_tag)))
		candidate_single_SV.append("%s\t%s\t%d\t%d\t%d\n"%(chr, svtype, breakpoint, indel_len, len(read_tag)))
	# print semi_indel_cluster
	# print max_conut_sum, max_conut

def resolution_TRA(path, chr_1, chr_2, read_count, overlap_size, max_cluster_bias):
	semi_tra_cluster = list()
	semi_tra_cluster.append([0,0,''])

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
			if len(semi_tra_cluster) > read_count:
				generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias)
			semi_tra_cluster = []
			semi_tra_cluster.append([pos_1, pos_2, read_id])
		else:
			semi_tra_cluster.append([pos_1, pos_2, read_id])

	if len(semi_tra_cluster) > read_count:
		generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias)

def generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias):
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
		print element[0], element[1]
	if len(read_tag) < read_count:
		return

	temp = sorted(temp, key = lambda x:-x[2])
	print temp

	if temp[1][2] > len(read_tag)*0.4:
		print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
		print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[1][0]/temp[1][2]), chr_2, int(temp[1][1]/temp[1][2]), len(read_tag)))
	else:
		if temp[0][2] > len(read_tag)*0.6:
			print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))

	# # result
	# # chr, "TRA", breakpoint_1, breakpoint_2, read_count
	# breakpoint_1 = int(breakpoint_ele / len(semi_tra_cluster))
	# breakpoint_2 = int(max_conut_sum/max_conut)
	# if max_conut < int(len(read_tag)*overlap_size):
	# 	return
	# print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, breakpoint_1, chr_2, breakpoint_2, len(read_tag)))
	# candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, breakpoint_1, chr_2, breakpoint_2, len(read_tag)))
	# # print semi_tra_cluster
	# # print max_conut_sum, max_conut

if __name__ == '__main__':
	candidate_single_SV = []
	# resolution_INDEL(sys.argv[1], "1", "INV", 10, 0.5, 20, 50)
	# file_out

	resolution_TRA(sys.argv[1], "1", "hs37d5", 10, 0.5, 50)

