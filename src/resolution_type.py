import sys

# candidate_single_SV = list()
# dup_candidates = dict()
# dup_candidates["l_DUP"] = dict()
# dup_candidates["DUP_r"] = dict()

def run_indel_inv(args):
	return resolution_INDEL(*args)

def run_tra(args):
	return resolution_TRA(*args)

def run_dup(args):
	return resolution_DUP(*args)

def resolution_INDEL(path, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size):
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
			if len(semi_indel_cluster) > read_count:
				generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, candidate_single_SV)
			semi_indel_cluster = []
			semi_indel_cluster.append([pos, indel_len, read_id])
		else:
			semi_indel_cluster.append([pos, indel_len, read_id])

	if len(semi_indel_cluster) > read_count:
		generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, candidate_single_SV)
	file.close()

	# file = open("%s%s_%s.bed"%(file_out, chr, svtype), 'w')
	# for i in candidate_single_SV:
	# 	file.write(i)
	# file.close()

	return candidate_single_SV

def generate_semi_indel_cluster(semi_indel_cluster, chr, svtype, read_count, overlap_size, max_cluster_bias, sv_size, candidate_single_SV):
	read_tag = dict()
	phase_indel_len = list()
	breakpoint_ele = 0
	for element in semi_indel_cluster:
		phase_indel_len.append(element[1])
		breakpoint_ele += element[0]
		if element[2] not in read_tag:
			read_tag[element[2]] = 0

	if len(read_tag) < read_count:
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
		#print("%s\t%s\t%d\t%d\t%d"%(chr, svtype, breakpoint, indel_len, len(read_tag)))
		candidate_single_SV.append("%s\t%s\t%d\t%d\t%d\n"%(chr, svtype, breakpoint, indel_len, len(read_tag)))
	else:
		#print("%s\t%s\t%d\t%d\t%d"%(chr, svtype, breakpoint, indel_len, len(read_tag)))
		candidate_single_SV.append("%s\t%s\t%d\t%d\t%d\n"%(chr, svtype, breakpoint, indel_len, len(read_tag)))
	# print semi_indel_cluster
	# print max_conut_sum, max_conut

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
			if len(semi_tra_cluster) > read_count:
				generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, max_cluster_bias, candidate_single_SV)
			semi_tra_cluster = []
			semi_tra_cluster.append([pos_1, pos_2, read_id])
		else:
			semi_tra_cluster.append([pos_1, pos_2, read_id])

	if len(semi_tra_cluster) > read_count:
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
		#print element[0], element[1]
	if len(read_tag) < read_count:
		return

	temp = sorted(temp, key = lambda x:-x[2])
	# print temp

	# if temp[1][2]+temp[0][2] > len(read_tag)*0.5:
	if temp[1][2] >= 0.5*read_count:
		if temp[0][2]+temp[1][2] >= len(semi_tra_cluster)*overlap_size:
			#print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			#print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[1][0]/temp[1][2]), chr_2, int(temp[1][1]/temp[1][2]), len(read_tag)))
			candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[1][0]/temp[1][2]), chr_2, int(temp[1][1]/temp[1][2]), len(read_tag)))
	else:
		if temp[0][2] >= len(semi_tra_cluster)*overlap_size:
			#print("%s\tTRA\t%d\t%s\t%d\t%d"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))
			candidate_single_SV.append("%s\tTRA\t%d\t%s\t%d\t%d\n"%(chr_1, int(temp[0][0]/temp[0][2]), chr_2, int(temp[0][1]/temp[0][2]), len(read_tag)))


def resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size=50):
	semi_dup_cluster = list()
	semi_dup_cluster.append([0,0,''])
	dup_candidates = dict()
	dup_candidates["l_DUP"] = dict()
	dup_candidates["DUP_r"] = dict()
	candidate_single_SV = list()

	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[1] != chr:
			continue

		if seq[0] != "DUP":
			# chr = seq[1]
			key_1 = int(seq[2])
			key_2 = int(seq[3])
			pos = int(seq[4])
			read_id = seq[5]
			if chr not in dup_candidates[seq[0]]:
				dup_candidates[seq[0]][chr] = dict()
				dup_candidates[seq[0]][chr][key_1] = dict()
				dup_candidates[seq[0]][chr][key_1][key_2] = list()
			else:
				if key_1 not in dup_candidates[seq[0]][chr]:
					dup_candidates[seq[0]][chr][key_1] = dict()
					dup_candidates[seq[0]][chr][key_1][key_2] = list()
				else:
					if key_2 not in dup_candidates[seq[0]][chr][key_1]:
						dup_candidates[seq[0]][chr][key_1][key_2] = list()
			dup_candidates[seq[0]][chr][key_1][key_2].append([pos, read_id])
		else:
			pos_1 = int(seq[2])
			pos_2 = int(seq[3])
			read_id = seq[4]
		
			if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias:
				if len(semi_dup_cluster) > 1:
					generate_semi_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, sv_size, dup_candidates, candidate_single_SV)
				semi_dup_cluster = []
				semi_dup_cluster.append([pos_1, pos_2, read_id])
			else:
				semi_dup_cluster.append([pos_1, pos_2, read_id])

	if len(semi_dup_cluster) > 1:
		generate_semi_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, sv_size, dup_candidates, candidate_single_SV)
	file.close()
	return candidate_single_SV

def acquire_locus(down, up, keytype, chr, MainCandidate):
	re = list()
	# re_id = list()
	if chr not in MainCandidate[keytype]:
		return re
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in MainCandidate[keytype][chr]:
			return re
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					re.append([ele[0],ele[1]])
					# re_id.append(ele[1])
	else:
		key_1 = int(down/10000)
		if key_1 in MainCandidate[keytype][chr]:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in MainCandidate[keytype][chr][key_1]:
					continue
				for ele in MainCandidate[keytype][chr][key_1][key_2]:
					if ele[0] >= down and ele[0] <= up:
						re.append([ele[0],ele[1]])
						# re_id.append(ele[1])
		key_1 += 1
		if key_1 not in MainCandidate[keytype][chr]:
			return re
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in MainCandidate[keytype][chr][key_1]:
				continue
			for ele in MainCandidate[keytype][chr][key_1][key_2]:
				if ele[0] >= down and ele[0] <= up:
					re.append([ele[0],ele[1]])
					# re_id.append(ele[1])
	return re

def cal_cluster_breakpoint(pos_list, bias):
	pos_list = sorted(pos_list, key = lambda x:x[0])
	#print pos_list
	temp = list()
	temp.append(pos_list[0])
	best_bp =list()
	len_best_bp = 0
	for pos in pos_list[1:]:
		if temp[-1][0] + bias < pos[0]:
			if len(temp) > len_best_bp:
				len_best_bp = len(temp)
				# best_bp = sum(temp)/len(temp)
				best_bp = temp
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if len(temp) > len_best_bp:
		len_best_bp = len(temp)
		# best_bp = sum(temp)/len(temp)
		best_bp = temp
	return best_bp

def generate_semi_dup_cluster(semi_dup_cluster, chr, read_count, max_cluster_bias, sv_size, dup_candidates, candidate_single_SV):
	start = list()
	end = list()
	start_id = list()
	end_id = list()
	read_tag = dict()
	for element in semi_dup_cluster:
		start.append(element[0])
		end.append(element[1])
		start_id.append([element[0], element[2]])
		end_id.append([element[1], element[2]])
		# if element[2] not in read_tag:
		# 	read_tag[element[2]] = 0

	# re_l, rl_id = acquire_locus(int(min(start)/1000)*1000, int(max(start)/1000+1)*1000, "l_DUP", chr, dup_candidates)
	# re_r, rr_id = acquire_locus(min(end)-1000, max(end)+1000, "DUP_r", chr, dup_candidates)
	re_l = acquire_locus(min(start)-max_cluster_bias, max(start)+max_cluster_bias, "l_DUP", chr, dup_candidates)
	re_r = acquire_locus(min(end)-max_cluster_bias, max(end)+max_cluster_bias, "DUP_r", chr, dup_candidates)
	info_1 = cal_cluster_breakpoint(start_id+re_l, max_cluster_bias)
	info_2 = cal_cluster_breakpoint(end_id+re_r, max_cluster_bias)

	cal_pos = 0
	for i in info_1:
		if i[1] not in read_tag:
			read_tag[i[1]] = 0
		cal_pos += i[0]
	breakpoint_1 = int(cal_pos/len(info_1))
	cal_pos = 0
	for i in info_2:
		if i[1] not in read_tag:
			read_tag[i[1]] = 0
		cal_pos += i[0]
	breakpoint_2 = int(cal_pos/len(info_2))	

	if breakpoint_2 - breakpoint_1 + 1 >= sv_size and len(read_tag) >= read_count:
		# print info_1
		# print info_2
		#print("%s\t%s\t%d\t%d\t%d"%(chr, "DUP", breakpoint_1, breakpoint_2-breakpoint_1+1, len(read_tag)))
		candidate_single_SV.append("%s\t%s\t%d\t%d\t%d\n"%(chr, "DUP", breakpoint_1, breakpoint_2-breakpoint_1+1, len(read_tag)))

#if __name__ == '__main__':
	# candidate_single_SV = []
	# resolution_INDEL(sys.argv[1], "1", "DEL", 10, 0.5, 20, 50)
	# file_out

	# resolution_TRA(sys.argv[1], "1", "hs37d5", 10, 0.6, 50)
	# resolution_DUP(sys.argv[1], "1", 10, 50, 50)
