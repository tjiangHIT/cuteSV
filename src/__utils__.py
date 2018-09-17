import cigar
import gc, sys
import time
# INS_flag = {1:'I'}
# DEL_flag = {2:'D'}

candidate = dict()
candidate["DEL"] = dict()
candidate["INS"] = dict()
candidate["INV"] = dict()
candidate["DUP"] = dict()
candidate["TRA"] = dict()
candidate["l_DUP"] = dict()
candidate["DUP_r"] = dict()
candidate["l_INV"] = dict()
candidate["INV_r"] = dict()

transfer = {"Indel_Ins":"INS", "INS":"INS", "TRA":"TRA", "DUP":"DUP", "l_DUP":"l_DUP", "DUP_r":"DUP_r", "INV":"INV", "Indel_Del":"DEL", "l_INV":"l_INV", "INV_r":"INV_r"}

semi_result = dict()

def detect_flag(Flag):
	# Signal
	Normal_foward = 1 >> 1
	Abnormal = 1 << 2
	Reverse_complement = 1 << 4
	Supplementary_map = 1 << 11

	signal = {Abnormal: 0, Normal_foward: 1, Reverse_complement: 2, Supplementary_map:3, Reverse_complement | Supplementary_map:4}
	if Flag in signal:
		return signal[Flag]
	else:
		return 0

def search_indel_str(deal_cigar, pos_start, SV_size, Chr_name, RLength):
	seq = list(cigar.Cigar(deal_cigar).items())
	# Ins_list = list()
	# Del_list = list()
	shift_ins = 0
	shift_del = 0
	for element in seq:
		if element[1] == 'M' or element[1] == 'D':
			shift_ins += element[0]
		if element[1] == 'I' and element[0] > SV_size:
			shift_ins += 1
			# Ins_list.append([Chr_name, pos_start + shift_ins, element[0]])
			store_signal([[Chr_name, pos_start + shift_ins, element[0]]], "INS")

		if element[1] == 'M':
			shift_del += element[0]
		if element[1] == 'D' and element[0] < SV_size:
			shift_del += element[0]
		if element[1] == 'D' and element[0] >= SV_size:
			# Del_list.append([Chr_name, pos_start + shift_del, element[0]])
			store_signal([[Chr_name, pos_start + shift_del, element[0]]], "DEL")
			shift_del += element[0]

	if seq[0][1] == 'S':
		softclip_left = seq[0][0]
	else:
		softclip_left = 0
	if seq[-1][1] == 'S':
		softclip_right = seq[-1][0]
	else:
		softclip_right = 0
	clip_list = [softclip_left, RLength - softclip_right, pos_start, pos_start + shift_ins, Chr_name]

	return clip_list

def search_indel_list(deal_cigar, pos_start, SV_size, Chr_name, RLength):
	# Ins_list = list()
	# Del_list = list()
	shift_ins = 0
	shift_del = 0
	# _shift_read_ = 0
	for element in deal_cigar:
		if element[0] == 0 or element[0] == 2:
			shift_ins += element[1]
		# if element[0] != 2:
		# 	_shift_read_ += element[1]
		if element[0] == 1 and element[1] > SV_size:
			shift_ins += 1
			# Ins_list.append([Chr_name, pos_start + shift_ins, element[1]])
			store_signal([[Chr_name, pos_start + shift_ins, element[1]]], "INS")
			
		if element[0] == 0:
			shift_del += element[1]
		if element[0] == 2 and element[1] < SV_size:
			shift_del += element[1]
		if element[0] == 2 and element[1] >= SV_size:
			# Del_list.append([Chr_name, pos_start + shift_del, element[1]])
			store_signal([[Chr_name, pos_start + shift_del, element[1]]], "DEL")
			shift_del += element[1]

	if deal_cigar[0][0] == 4:
		softclip_left = deal_cigar[0][1]
	else:
		softclip_left = 0
	if deal_cigar[-1][0] == 4:
		softclip_right = deal_cigar[-1][1]
	else:
		softclip_right = 0

	clip_list = [softclip_left, RLength - softclip_right, pos_start + 1, pos_start + shift_ins + 1, Chr_name]
	return clip_list

def store_info(pos_list, svtype):
	for ele in pos_list:
		if ele[0] not in candidate[svtype]:
			candidate[svtype][ele[0]] = dict()

		hash_1 = int(ele[1] /10000)
		mod = ele[1] % 10000
		hash_2 = int(mod / 50)

		if hash_1 not in candidate[svtype][ele[0]]:
			candidate[svtype][ele[0]][hash_1] = dict()
			candidate[svtype][ele[0]][hash_1][hash_2] = list()
			candidate[svtype][ele[0]][hash_1][hash_2].append(ele[1:])
		else:
			if hash_2 not in candidate[svtype][ele[0]][hash_1]:
				candidate[svtype][ele[0]][hash_1][hash_2] = list()

			candidate[svtype][ele[0]][hash_1][hash_2].append(ele[1:])

def store_signal(pos_list, svtype):
	for ele in pos_list:
		if ele[0] not in candidate[svtype]:
			candidate[svtype][ele[0]] = list()

		candidate[svtype][ele[0]].append(ele[1:])

def analysis_split_read(split_read, SV_size, RLength):
	# +indel++indel++indel++indel++indel++indel+
	SP_list = list()
	for read in split_read:
		if isinstance(read[3], str):
			sp_list = search_indel_str(read[3], read[1], SV_size, read[0], RLength)
		else:
			sp_list = search_indel_list(read[3], read[1], SV_size, read[0], RLength)
		# store_signal(Del_list, "DEL")
		# store_signal(Ins_list, "INS")

		sp_list += read[2]
		SP_list.append(sp_list)

	# split alignment
	if len(SP_list) == 1:
		return 0
	else:
		SP_list = sorted(SP_list, key = lambda x:x[0])
	# for i in SP_list:
	#  	print i

	DUP_flag = [0]*len(SP_list)
	# INVDUP_flag = [0]*len(SP_list)
	# INS_flag = [0]*len(SP_list)

	# for a in SP_list:
	# 	for b in SP_list:
	# 		if a[4] == b[4]:
	# 			# dup & INV & TRA & INS & DEL
	# 			if b[3] - a[2] >= SV_size and SP_list.index(a) > SP_list.index(b):
	# 				DUP_flag[SP_list.index(a)] = 1
	# 				DUP_flag[SP_list.index(b)] = 1
	# 			if a[0] + b[3] - a[2] - b[1] >= SV_size and b[3] <= a[2] and SP_list.index(a) == SP_list.index(b) + 1:
	# 				store_signal([[a[4], (a[2]+b[3])/2, a[0]+b[3]-a[2]-b[1]]], "INS")
	# 			if a[2] - a[0] + b[1] - b[3] >= SV_size and b[3] <= a[2] and SP_list.index(a) == SP_list.index(b) + 1:
	# 				store_signal([[a[4], b[3], a[2]-a[0]+b[1]-b[3]]], "DEL")
	# 		else:
	# 			# tra
	# 			if SP_list.index(a) > SP_list.index(b):
	# 				if b[4] < a[4]:
	# 					store_signal([[b[4], b[3], a[4], a[2]]], "TRA")
	# 				else:
	# 					store_signal([[a[4], a[2], b[4], b[3]]], "TRA")

	for a in xrange(len(SP_list[:-1])):
		ele_1 = SP_list[a]
		ele_2 = SP_list[a+1]
		if ele_1[4] == ele_2[4]:
			# dup & ins & del 
			if ele_1[3] - ele_2[2] >= SV_size:
				# DUP_flag[a] = 1
				# DUP_flag[a+1] = 1
				if a == 0:
					store_info([[ele_1[4], ele_1[3]]], "DUP_r")

				if a+1 == len(SP_list)-1:
					store_info([[ele_2[4], ele_2[2]]], "l_DUP")
				else:
					store_signal([[ele_2[4], ele_2[2], ele_2[3]]], "DUP")

			if ele_1[3] <= ele_2[2]:
				if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
					store_signal([[ele_2[4], (ele_2[2]+ele_1[3])/2, ele_2[0]+ele_1[3]-ele_2[2]-ele_1[1]]], "INS")
				if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
					store_signal([[ele_2[4], ele_1[3], ele_2[2]-ele_2[0]+ele_1[1]-ele_1[3]]], "DEL")
		else:
			# tra
			if ele_1[4] < ele_2[4]:
				store_signal([[ele_1[4], ele_1[3], ele_2[4], ele_2[2]]], "TRA")
			else:
				store_signal([[ele_2[4], ele_2[2], ele_1[4], ele_1[3]]], "TRA")



	# for k in xrange(len(SP_list)):
	# 	if DUP_flag[k] == 1:
	# 		if k == 0:
	# 			store_info([[SP_list[k][4], SP_list[k][3]]], "DUP_r")
	# 		elif k == len(SP_list)-1:
	# 			store_info([[SP_list[k][4], SP_list[k][2]]], "l_DUP")
	# 		else:
	# 			store_signal([[SP_list[k][4], SP_list[k][2], SP_list[k][3]]], "DUP")

	# +TRA++TRA++TRA++TRA++TRA++TRA++TRA++TRA+

	# +INV++INV++INV++INV++INV++INV++INV++INV+
	call_inv = sorted(SP_list, key = lambda x:x[2])
	if len(call_inv) >= 3:
		for a in call_inv[:-2]:
			if a[5] != call_inv[call_inv.index(a)+1][5] and a[5] == call_inv[call_inv.index(a)+2][5] and a[4] == call_inv[call_inv.index(a)+1][4] and a[4] == call_inv[call_inv.index(a)+2][4]:
				if call_inv[call_inv.index(a)+1][3] - call_inv[call_inv.index(a)+1][2] >= SV_size:
					store_signal([[a[4], a[3] , call_inv[call_inv.index(a)+1][3]]], "INV")

	
	if len(call_inv) == 2:
		if call_inv[0][5] != call_inv[1][5] and call_inv[0][4] == call_inv[1][4]:
			ls_1 = call_inv[0][3] - call_inv[0][2]
			ls_2 = call_inv[1][3] - call_inv[1][2]
			if ls_1 > ls_2:
				if call_inv[1][2] > call_inv[0][3] and ls_2 >= SV_size:
					store_info([[call_inv[0][4], call_inv[1][3]]], "INV_r")
					# pass
			else:
				if call_inv[1][2] > call_inv[0][3] and ls_1 >= SV_size:
					store_info([[call_inv[0][4], call_inv[0][2]]], "l_INV")
					# pass

def parse_read(read, Chr_name, SV_size, MQ_threshold, max_num_splits, min_seq_size):
	process_signal = detect_flag(read.flag)
	if process_signal == 0:
		# unmapped reads
		# return INS_ME_pos
		pass
	# split alignment phasing
	if process_signal == 1 or process_signal == 2:
		split_read = list()
		if read.mapq > MQ_threshold:
			if read.is_reverse:
				strand = '-'
			else:
				strand = '+'
			split_read.append([Chr_name, read.reference_start, strand, read.cigar])
		Tags = read.get_tags()
		for tag in Tags:
			if tag[0] == 'SA':
				split_alignment = tag[1].split(';')[:-1]
				for mapping in split_alignment:
					parse_mapping = mapping.split(',')
					if int(parse_mapping[4]) < MQ_threshold:
						continue
					split_read.append([parse_mapping[0], int(parse_mapping[1]), parse_mapping[2], parse_mapping[3]])
				break

		# print read.query_name
		if len(split_read) <= max_num_splits and read.query_length >= min_seq_size:
			analysis_split_read(split_read, SV_size, read.query_length)
		gc.collect()

def merge_pos_indel(pos_list, chr, evidence_read, SV_size, svtype):
	# result = list()
	if len(pos_list) >= evidence_read:
		start = list()
		size = list()
		# max_L = -1
		# min_L = sys.maxint
		for ele in pos_list:
			start.append(ele[0])
			size.append(ele[1])
			# if ele[1] > max_L:
			# 	max_L = ele[1]
			# if ele[1] < min_L:
			# 	min_L = ele[1]
		size.sort()
		breakpoint = sum(start)/len(start)
		SV_len = sum(size[int(0.25*len(size)):int(0.75*len(size))])/len(size[int(0.25*len(size)):int(0.75*len(size))])
		# concider the stability of data
		# need to repair the buf!!!!!!!!

		if SV_len >= SV_size and size[int(0.75*len(size))] - size[int(0.25*len(size))] < int(0.5*SV_len):
			# result.append([chr, breakpoint, SV_len, len(pos_list)])
			if chr not in semi_result:
				semi_result[chr] = list()
			# semi_result[chr].append([breakpoint, SV_len, len(pos_list), svtype])
			semi_result[chr].append("%d\t%d\t%d\t%s\n"%(breakpoint, SV_len, len(pos_list), svtype))

	# return result

def intergrate_indel(chr, evidence_read, SV_size, low_bandary, svtype, max_distance):
	# _cluster_ = list()
	temp = list()
	temp.append(candidate[svtype][chr][0])
	for pos in candidate[svtype][chr][1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			if temp[-1][0] - temp[0][0] <= max_distance:
				merge_pos_indel(temp, chr, evidence_read, SV_size, svtype)
			# if len(result) != 0:
			# 	_cluster_.append(result[0])
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if temp[-1][0] - temp[0][0] <= max_distance:
		merge_pos_indel(temp, chr, evidence_read, SV_size, svtype)
	# if len(result) != 0:
	# 	_cluster_.append(result[0])
	# return _cluster_

def acquire_locus(down, up, keytype, chr):
	re = list()
	if chr not in candidate[keytype]:
		return re
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in candidate[keytype][chr]:
			return re
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in candidate[keytype][chr][key_1]:
				continue
			for ele in candidate[keytype][chr][key_1][key_2]:
				if ele >= down and ele <= up:
					re.append(ele)
	else:
		key_1 = int(down/10000)
		if key_1 in candidate[keytype][chr]:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in candidate[keytype][chr][key_1]:
					continue
				for ele in candidate[keytype][chr][key_1][key_2]:
					if ele >= down and ele <= up:
						re.append(ele)
		key_1 += 1
		if key_1 not in candidate[keytype][chr]:
			return re
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in candidate[keytype][chr][key_1]:
				continue
			for ele in candidate[keytype][chr][key_1][key_2]:
				if ele >= down and ele <= up:
					re.append(ele)
	return re

def merge_pos_dup(pos_list, chr, evidence_read, SV_size, ll, lr, svtype):
	start = list()
	end = list()
	up_B = int(len(pos_list)*0.75)
	low_B = int(len(pos_list)*0.25)
	for ele in pos_list:
		start.append(ele[0])
	for ele in pos_list[low_B:up_B+1]:
		end.append(ele[1])

	re_l = acquire_locus(min(start), max(start), ll, chr)
	re_r = acquire_locus(min(end), max(end), lr, chr)

	breakpoint = sum(start+re_l)/len(start+re_l)
	size = sum(end+re_r)/len(end+re_r) - breakpoint

	# print breakpoint, size, re_r, re_l, start, end

	# result = list()
	if len(pos_list+re_l+re_r) >= evidence_read:
		# need to add some semi-signals
		if size >= SV_size and max(end+re_r) - min(end+re_r) < int(0.5*size):
			# result.append([chr, breakpoint, size, len(pos_list+re_l+re_r)])
			if chr not in semi_result:
				semi_result[chr] = list()
			# semi_result[chr].append([breakpoint, size, len(pos_list+re_l+re_r), svtype])
			semi_result[chr].append("%d\t%d\t%d\t%s\n"%(breakpoint, size, len(pos_list+re_l+re_r), svtype))

def intergrate_dup(chr, evidence_read, SV_size, low_bandary, ll, lr, svtype, max_distance):
	temp = list()
	temp.append(candidate[svtype][chr][0])
	for pos in candidate[svtype][chr][1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			if temp[-1][0] - temp[0][0] <= max_distance:
				merge_pos_dup(temp, chr, evidence_read, SV_size, ll, lr, svtype)
				# print temp
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if temp[-1][0] - temp[0][0] <= max_distance:
		merge_pos_dup(temp, chr, evidence_read, SV_size, ll, lr, svtype)
		# print temp

def mkindex(data, mask):
	data_struc = dict()
	for i in data:
		# if i[mask]
		hash_1 = int(i[mask] /10000)
		mod = i[mask] % 10000
		hash_2 = int(mod / 50)
		# element = [i[0], seq, flag]
		if hash_1 not in data_struc:
			data_struc[hash_1] = dict()
			data_struc[hash_1][hash_2] = list()
			data_struc[hash_1][hash_2].append(i[mask])
		else:
			if hash_2 not in data_struc[hash_1]:
				data_struc[hash_1][hash_2] = list()
				data_struc[hash_1][hash_2].append(i[mask])
			else:
				data_struc[hash_1][hash_2].append(i[mask])
	return data_struc

def merge_pos_tra(pos_list, chr, evidence_read, SV_size, svtype, low_bandary):
	if len(pos_list) >= evidence_read:
		start = list()
		tra_dic = dict()
		for ele in pos_list:
			start.append(ele[0])
			# size.append(ele[1])
			if ele[1] not in tra_dic:
				tra_dic[ele[1]] = list()
			tra_dic[ele[1]].append(ele[2])

		max_chr_c = 0
		max_chr = ""
		for key in tra_dic:
			if len(tra_dic[key]) > max_chr_c:
				max_chr = key
				max_chr_c = len(tra_dic[key])

		if max_chr_c >= evidence_read:
			breakpoint_1 = int(sum(start)/len(start))
			breakpoint_2 = int(sum(tra_dic[max_chr])/max_chr_c)

			if chr not in semi_result:
				semi_result[chr] = list()
			# semi_result[chr].append([breakpoint_1, max_chr, breakpoint_2, svtype])
			semi_result[chr].append("%d\t%s\t%s\t%s\n"%(breakpoint_1, max_chr, breakpoint_2, svtype))

def intergrate_tra(chr, evidence_read, SV_size, low_bandary, svtype, max_distance):
	temp = list()
	temp.append(candidate[svtype][chr][0])
	for pos in candidate[svtype][chr][1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			if temp[-1][0] - temp[0][0] <= max_distance:
				merge_pos_tra(temp, chr, evidence_read, SV_size, svtype, low_bandary)
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	if temp[-1][0] - temp[0][0] <= max_distance:
		merge_pos_tra(temp, chr, evidence_read, SV_size, svtype, low_bandary)

def show_temp_result(evidence_read, SV_size, low_bandary, max_distance, out_path):
	starttime = time.time()
	for chr in candidate["DEL"]:
		candidate["DEL"][chr] = sorted(candidate["DEL"][chr], key = lambda x:x[0])
		# for ele in candidate["DEL"][chr]:
		# 	print chr, ele
		intergrate_indel(chr, evidence_read, SV_size, low_bandary, "DEL", max_distance)
		candidate["DEL"][chr] = []
		gc.collect()
	print("[INFO]: Parse deletions used %0.2f seconds."%(time.time() - starttime))

	starttime = time.time()
	for chr in candidate["INS"]:
		candidate["INS"][chr] = sorted(candidate["INS"][chr], key = lambda x:x[0])
		# for i in candidate["INS"][chr]:
		# 	print chr, i
		intergrate_indel(chr, evidence_read, SV_size, low_bandary, "INS", max_distance)
		candidate["INS"][chr] = []
		gc.collect()
	print("[INFO]: Parse insertions used %0.2f seconds."%(time.time() - starttime))

	starttime = time.time()
	for chr in candidate["DUP"]:
		candidate["DUP"][chr] = sorted(candidate["DUP"][chr], key = lambda x:x[:])
		# for i in candidate["DUP"][chr]:
		# 	print chr, i
		intergrate_dup(chr, evidence_read, SV_size, low_bandary, "l_DUP", "DUP_r", "DUP", max_distance)
		candidate["DUP"][chr] = []
		candidate["l_DUP"][chr] = []
		candidate["DUP_r"][chr] = []
		gc.collect()
	print("[INFO]: Parse duplications used %0.2f seconds."%(time.time() - starttime))
	
	starttime = time.time()
	for chr in candidate["INV"]:
		candidate["INV"][chr] = sorted(candidate["INV"][chr], key = lambda x:x[:])
		intergrate_dup(chr, evidence_read, SV_size, low_bandary, "l_INV", "INV_r", "INV", max_distance)
		candidate["INV"][chr] = []
		candidate["l_INV"][chr] = []
		candidate["INV_r"][chr] = []
		gc.collect()
	print("[INFO]: Parse inversions used %0.2f seconds."%(time.time() - starttime))

	starttime = time.time()
	for chr in candidate["TRA"]:
		candidate["TRA"][chr] = sorted(candidate["TRA"][chr], key = lambda x:x[:])
		# for i in candidate["TRA"][chr]:
		# 	print chr, i
		intergrate_tra(chr, evidence_read, SV_size, low_bandary, "TRA", max_distance)
		candidate["TRA"][chr] = []
		gc.collect()
	print("[INFO]: Parse translocations used %0.2f seconds."%(time.time() - starttime))

	file = open(out_path, 'w')
	for chr in semi_result:
		semi_result[chr] = sorted(semi_result[chr], key = lambda x:x[0])
		for i in semi_result[chr]:
			# print i
			file.write("%s\t%s"%(chr, i))
	file.close()
	# print candidate['DEL']
