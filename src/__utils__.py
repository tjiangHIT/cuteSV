import cigar
import gc, sys
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
	Ins_list = list()
	Del_list = list()
	shift_ins = 0
	shift_del = 0
	for element in seq:
		if element[1] == 'M' or element[1] == 'D':
			shift_ins += element[0]
		if element[1] == 'I' and element[0] > SV_size:
			shift_ins += 1
			Ins_list.append([Chr_name, pos_start + shift_ins, element[0]])

		if element[1] == 'M':
			shift_del += element[0]
		if element[1] == 'D' and element[0] < SV_size:
			shift_del += element[0]
		if element[1] == 'D' and element[0] >= SV_size:
			Del_list.append([Chr_name, pos_start + shift_del, element[0]])
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

	return Ins_list, Del_list, clip_list

def search_indel_list(deal_cigar, pos_start, SV_size, Chr_name, RLength):
	Ins_list = list()
	Del_list = list()
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
			Ins_list.append([Chr_name, pos_start + shift_ins, element[1]])
			
		if element[0] == 0:
			shift_del += element[1]
		if element[0] == 2 and element[1] < SV_size:
			shift_del += element[1]
		if element[0] == 2 and element[1] >= SV_size:
			Del_list.append([Chr_name, pos_start + shift_del, element[1]])
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
	return Ins_list, Del_list, clip_list

def analysis_split_read(split_read, SV_size, RLength):
	candidate_SV = dict()
	# +indel++indel++indel++indel++indel++indel+
	candidate_SV["Indel_Ins"] = list()
	candidate_SV["Indel_Del"] = list()
	SP_list = list()

	for read in split_read:
		if isinstance(read[3], str):
			Ins_list, Del_list, sp_list = search_indel_str(read[3], read[1], SV_size, read[0], RLength)
		else:
			Ins_list, Del_list, sp_list = search_indel_list(read[3], read[1], SV_size, read[0], RLength)
		candidate_SV["Indel_Del"] += Del_list
		candidate_SV["Indel_Ins"] += Ins_list
		sp_list += read[2]
		SP_list.append(sp_list)
	'''
	for i in candidate_SV["Indel_Ins"]:
		print "INS", i
	for i in candidate_SV["Indel_Del"]:
		print "DEL", i
	'''
	# split alignment
	SP_list = sorted(SP_list, key = lambda x:x[0])
	# for i in SP_list:
	#  	print i

	candidate_SV["DUP"] = list()
	candidate_SV["l_DUP"] = list()
	candidate_SV["DUP_r"] = list()
	# candidate_SV["INVDUP"] = list()
	candidate_SV["INS"] = list()
	candidate_SV["INV"] = list()
	candidate_SV["l_INV"] = list()
	candidate_SV["INV_r"] = list()
	candidate_SV["TRA"] = list()
	DUP_flag = [0]*len(SP_list)
	INVDUP_flag = [0]*len(SP_list)
	# INS_flag = [0]*len(SP_list)

	for a in SP_list:
		for b in SP_list:
			if a[4] == b[4]:
				# dup & INV & TRA & INS
				if b[3] - a[2] >= SV_size and SP_list.index(a) > SP_list.index(b):
					# dup
					# if a[5] == b[5]:
					# 	DUP_flag[SP_list.index(a)] = 1
					# 	DUP_flag[SP_list.index(b)] = 1
					# else:
					# 	INVDUP_flag[SP_list.index(a)] = 1
					# 	INVDUP_flag[SP_list.index(b)] = 1
					DUP_flag[SP_list.index(a)] = 1
					DUP_flag[SP_list.index(b)] = 1

				if a[0] + b[3] - a[2] - b[1] >= SV_size and b[3] <= a[2] and SP_list.index(a) == SP_list.index(b) + 1:
					candidate_SV["INS"].append([a[4], (a[2]+b[3])/2, a[0]+b[3]-a[2]-b[1]])

				# if b[1] <= a[0] and b[3] <= a[2] and b[5] != a[5] and SP_list.index(a) == SP_list.index(b) + 1:
				# 	if SP_list.index(a) + 2 <= len(SP_list):
				# 		# print SP_list[SP_list.index(a)+1]
				# 		if a[1] <= SP_list[SP_list.index(a)+1][0] and a[3] <= SP_list[SP_list.index(a)+1][2] and a[5] != SP_list[SP_list.index(a)+1][5]:
				# 			candidate_SV["INV"].append([a[4], (a[2]+b[3])/2, (SP_list[SP_list.index(a)+1][2]+b[3]-a[2]-a[3])/2])

			else:
				# tra
				if SP_list.index(a) > SP_list.index(b):
					if b[4] < a[4]:
						candidate_SV["TRA"].append([b[4], b[3], a[4], a[2]])
					else:
						candidate_SV["TRA"].append([a[4], a[2], b[4], b[3]])
	# for k in xrange(len(SP_list[:-1])):
	# 	for l in xrange(len(SP_list[1:])):
	# 		if SP_list[k][4] == SP_list[l][4] and SP_list[k][2] < SP_list[l][3]:
	# 			# chr & ovlapping cmp
	# 			if SP_list[k][5] == SP_list[l][5]:
	# 				DUP_flag[k] = 1
	# 				DUP_flag[l] = 1
	# 			else:
	# 				INVDUP_flag[k] = 1
	# 				INVDUP_flag[l] = 1

			# if SP_list[k+1][2] < SP_list[k][2]:
			# 	# overlap cmp
			# 	if SP_list[k+1][5] == SP_list[k][5]:
			# 		# strand cmp
			# 		if k == len(SP_list)-2:
			# 			local_candidate_DUP.append([SP_list[k][4], SP_list[k+1][2], -1])
			# 		else:
			# 			local_candidate_DUP.append([SP_list[k][4], SP_list[k+1][2], SP_list[k+1][3]])
					
			# 		if k == 0:
			# 			local_candidate_DUP.append([SP_list[k][4], -1, SP_list[k][3]])
			# 	else:
			# 		if k == len(SP_list)-2:
			# 			local_candidate_INVDUP.append([SP_list[k][4], SP_list[k+1][2], -1])
			# 		else:
			# 			local_candidate_INVDUP.append([SP_list[k][4], SP_list[k+1][2], SP_list[k+1][3]])
			# 		if k == 0:
			# 			local_candidate_INVDUP.append([SP_list[k][4], -1, SP_list[k][3]])

	# +DUP++DUP++DUP++DUP++DUP++DUP++DUP++DUP+
	# temp_dup_s = dict()
	# temp_invdup_s = dict()
	# temp_dup_e = dict()
	# temp_invdup_e = dict()
	for k in xrange(len(SP_list)):
		if DUP_flag[k] == 1:
			# if SP_list[k][4] not in temp_dup_e:
			# 	temp_dup_e[SP_list[k][4]] = list()
			# if SP_list[k][4] not in temp_dup_s:
			# 	temp_dup_s[SP_list[k][4]] = list()
			if k == 0:
				candidate_SV["DUP_r"].append([SP_list[k][4], -1, SP_list[k][3]])
				# pass
				# temp_dup_e.append(SP_list[k][3])
				# temp_dup_e[SP_list[k][4]].append(SP_list[k][3])
			elif k == len(SP_list)-1:
				candidate_SV["l_DUP"].append([SP_list[k][4], SP_list[k][2], -1])
				# pass
				# temp_dup_s.append(SP_list[k][2])
				# temp_dup_s[SP_list[k][4]].append(SP_list[k][2])
			else:
				candidate_SV["DUP"].append([SP_list[k][4], SP_list[k][2], SP_list[k][3]])
				# temp_dup_s.append(SP_list[k][2])
				# temp_dup_s[SP_list[k][4]].append(SP_list[k][2])
				# temp_dup_e.append(SP_list[k][3])
				# temp_dup_e[SP_list[k][4]].append(SP_list[k][3])
		if INVDUP_flag[k] == 1:
			# if SP_list[k][4] not in temp_invdup_e:
			# 	temp_invdup_e[SP_list[k][4]] = list()
			# if SP_list[k][4] not in temp_invdup_s:
			# 	temp_invdup_s[SP_list[k][4]] = list()
			if k == 0:
				candidate_SV["INVDUP"].append([SP_list[k][4], -1, SP_list[k][3]])
				# temp_invdup_e.append(SP_list[k][3])
				# temp_invdup_e[SP_list[k][4]].append(SP_list[k][3])
			elif k == len(SP_list)-1:
				candidate_SV["INVDUP"].append([SP_list[k][4], SP_list[k][2], -1])
				# temp_invdup_s.append(SP_list[k][2])
				# temp_invdup_s[SP_list[k][4]].append(SP_list[k][2])
			else:
				candidate_SV["INVDUP"].append([SP_list[k][4], SP_list[k][2], SP_list[k][3]])
				# temp_invdup_s.append(SP_list[k][2])
				# temp_invdup_s[SP_list[k][4]].append(SP_list[k][2])
				# temp_invdup_e.append(SP_list[k][3])
				# temp_invdup_e[SP_list[k][4]].append(SP_list[k][3])

	# for key in temp_dup_s:
	# 	candidate_SV["DUP"].append([key, int(sum(temp_dup_s[key])/len(temp_dup_s[key])), int(sum(temp_dup_e[key])/len(temp_dup_e[key])), max(len(temp_dup_e[key]), len(temp_dup_s[key]))])
	# for key in temp_invdup_s:
	# 	candidate_SV["INVDUP"].append([key, int(sum(temp_invdup_s[key])/len(temp_invdup_s[key])), int(sum(temp_invdup_e[key])/len(temp_invdup_e[key])), max(len(temp_invdup_e[key]), len(temp_invdup_s[key]))])
	# +TRA++TRA++TRA++TRA++TRA++TRA++TRA++TRA+

	# +INV++INV++INV++INV++INV++INV++INV++INV+
	temp_inv_s = dict()
	temp_inv_e = dict()
	call_inv = sorted(SP_list, key = lambda x:x[2])
	if len(call_inv) >= 3:
		for a in call_inv[:-2]:
			if a[5] != call_inv[call_inv.index(a)+1][5] and a[5] == call_inv[call_inv.index(a)+2][5] and a[4] == call_inv[call_inv.index(a)+1][4] and a[4] == call_inv[call_inv.index(a)+2][4]:
				if call_inv[call_inv.index(a)+1][3] - call_inv[call_inv.index(a)+1][2] >= SV_size:
					candidate_SV["INV"].append([a[4], a[3] , call_inv[call_inv.index(a)+1][3]])
					# if a[4] not in temp_inv_s:
					# 	temp_inv_s[a[4]] = list()
					# if a[4] not in temp_inv_e:
					# 	temp_inv_e[a[4]] = list()
					# temp_inv_s[a[4]].append(a[3])
					# temp_inv_e[a[4]].append(call_inv[call_inv.index(a)+1][3])
	
	if len(call_inv) == 2:
		# if call_inv[0][4] not in temp_inv_s:
		# 	temp_inv_s[call_inv[0][4]] = list()
		# if call_inv[0][4] not in temp_inv_e:
		# 	temp_inv_e[call_inv[0][4]] = list()
		if call_inv[0][5] != call_inv[1][5] and call_inv[0][4] == call_inv[1][4]:
			ls_1 = call_inv[0][3] - call_inv[0][2]
			ls_2 = call_inv[1][3] - call_inv[1][2]
			if ls_1 > ls_2:
				if call_inv[1][2] > call_inv[0][3] and ls_2 >= SV_size:
					# temp_inv_e[call_inv[0][4]].append(call_inv[1][3])
					candidate_SV["INV_r"].append([call_inv[0][4], -1 , call_inv[1][3]])
					# pass
			else:
				if call_inv[1][2] > call_inv[0][3] and ls_1 >= SV_size:
					# temp_inv_s[call_inv[0][4]].append(call_inv[0][2])
					candidate_SV["l_INV"].append([call_inv[0][4], call_inv[0][2], -1])
					# pass
	

	# for key in temp_inv_s:
	# 	try:
	# 		candidate_SV["INV"].append([key, int(sum(temp_inv_s[key])/len(temp_inv_s[key])), int(sum(temp_inv_e[key])/len(temp_inv_e[key]))])
	# 	except:
	# 		# pass
	# 		if len(temp_inv_s[key]) == 0 :
	# 			candidate_SV["INV"].append([key, -1, int(sum(temp_inv_e[key])/len(temp_inv_e[key]))])
	# 		if len(temp_inv_e[key]) == 0:
	# 			candidate_SV["INV"].append([key, int(sum(temp_inv_s[key])/len(temp_inv_s[key])), -1])

	# print candidate_SV["DUP"]
	return candidate_SV

def parse_read(read, Chr_name, SV_size, MQ_threshold):
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
		# for key in split_read:
			# print isinstance(key[3], str), key[3]
		data = analysis_split_read(split_read, SV_size, read.query_length)
		gc.collect()

		for key in data:
			# if key not in candidate:
			# 	candidate[key] = list()
			for i in data[key]:
				# print [key, i, transfer[key]]
				if i[0] not in candidate[transfer[key]]:
					candidate[transfer[key]][i[0]] = list()
				candidate[transfer[key]][i[0]].append(i[1:])
				# print [key, i, transfer[key]]

def merge_pos_indel(pos_list, chr, evidence_read, SV_size):
	start = list()
	end = list()
	max_L = -1
	min_L = sys.maxint
	for ele in pos_list:
		start.append(ele[0])
		end.append(ele[0] + ele[1])
		if ele[1] > max_L:
			max_L = ele[1]
		if ele[1] < min_L:
			min_L = ele[1]

	breakpoint = sum(start)/len(start)
	size = sum(end)/len(end) - breakpoint

	result = list()
	if len(pos_list) < evidence_read:
		return result
	else:
		if size >= SV_size and max_L - min_L < int(0.5*size):
			result.append([chr, breakpoint, size, len(pos_list)])

	return result

def intergrate_indel(pos_list, chr, evidence_read, SV_size, low_bandary):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			result = merge_pos_indel(temp, chr, evidence_read, SV_size)
			if len(result) != 0:
				_cluster_.append(result[0])
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos_indel(temp, chr, evidence_read, SV_size)
	if len(result) != 0:
		_cluster_.append(result[0])
	return _cluster_

def acquire_locus(down, up, data_struc):
	re = list()
	if int(up/10000) == int(down/10000):
		key_1 = int(down/10000)
		if key_1 not in data_struc:
			return re
		for i in xrange(int((up%10000)/50)-int((down%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = int((down%10000)/50)+i
			if key_2 not in data_struc[key_1]:
				continue
			for ele in data_struc[key_1][key_2]:
				if ele >= down and ele <= up:
					re.append(ele)
	else:
		key_1 = int(down/10000)
		if key_1 in data_struc:
			for i in xrange(200-int((down%10000)/50)):
				# exist a bug ***********************************
				key_2 = int((down%10000)/50)+i
				if key_2 not in data_struc[key_1]:
					continue
				for ele in data_struc[key_1][key_2]:
					if ele >= down and ele <= up:
						re.append(ele)
		key_1 += 1
		if key_1 not in data_struc:
			return re
		for i in xrange(int((up%10000)/50)+1):
			# exist a bug ***********************************
			key_2 = i
			if key_2 not in data_struc[key_1]:
				continue
			for ele in data_struc[key_1][key_2]:
				if ele >= down and ele <= up:
					re.append(ele)
	return re

def merge_pos_dup(pos_list, chr, evidence_read, SV_size, ll, lr):
	start = list()
	end = list()
	up_B = int(len(pos_list)*0.75)
	low_B = int(len(pos_list)*0.25)
	for ele in pos_list:
		start.append(ele[0])
	for ele in pos_list[low_B:up_B+1]:
		end.append(ele[1])

	re_l = acquire_locus(min(start), max(start), ll)
	re_r = acquire_locus(min(end), max(end), lr)

	breakpoint = sum(start+re_l)/len(start+re_l)
	size = sum(end+re_r)/len(end+re_r) - breakpoint

	result = list()
	if len(pos_list+re_l+re_r) < evidence_read:
		# need to add some semi-signals
		return result
	else:
		if size >= SV_size and max(end+re_r) - min(end+re_r) < int(0.5*size):
			result.append([chr, breakpoint, size, len(pos_list+re_l+re_r)])

	return result

def intergrate_dup(pos_list, chr, evidence_read, SV_size, low_bandary, ll, lr):
	_cluster_ = list()
	temp = list()
	temp.append(pos_list[0])
	for pos in pos_list[1:]:
		if temp[-1][0] + low_bandary < pos[0]:
			result = merge_pos_dup(temp, chr, evidence_read, SV_size, ll, lr)
			if len(result) != 0:
				_cluster_.append(result[0])
			temp = list()
			temp.append(pos)
		else:
			temp.append(pos)
	result = merge_pos_dup(temp, chr, evidence_read, SV_size, ll, lr)
	if len(result) != 0:
		_cluster_.append(result[0])
	return _cluster_

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

def show_temp_result():
	for key in candidate["DEL"]:
		candidate["DEL"][key] = sorted(candidate["DEL"][key], key = lambda x:x[:])
		# for i in candidate["DEL"][key]:
		#	print i
		result = intergrate_indel(candidate["DEL"][key], key, 5, 50, 10)
		# for i in result:
		#	print "DEL", i
		if key not in semi_result:
			semi_result[key] = list()
		for i in result:
			i.append("DEL")
			semi_result[key].append(i)

	for key in candidate["INS"]:
		candidate["INS"][key] = sorted(candidate["INS"][key], key = lambda x:x[:])
		# for i in candidate["INS"][key]:
		# 	print key, i
		result = intergrate_indel(candidate["INS"][key], key, 5, 50, 10)
		# for i in result:
		#	print "INS", i
		if key not in semi_result:
			semi_result[key] = list()
		for i in result:
			i.append("INS")
			semi_result[key].append(i)

	for key in candidate["DUP"]:
		candidate["DUP"][key] = sorted(candidate["DUP"][key], key = lambda x:x[:])
		# candidate["l_DUP"][key] = sorted(candidate["l_DUP"][key], key = lambda x:x[:])
		# candidate["DUP_r"][key] = sorted(candidate["DUP_r"][key], key = lambda x:x[:])
		try:
			left_dup = mkindex(candidate["l_DUP"][key], 0)
		except:
			left_dup = dict()
		try:
			right_dup = mkindex(candidate["DUP_r"][key], 1)
		except:
			right_dup = dict()

		# for i in candidate["l_DUP"][key]:
		# 	print key, i
		# for i in candidate["DUP"][key]:
		# 	print key, i
		# for i in candidate["DUP_r"][key]:
		# 	print key, i
		result = intergrate_dup(candidate['DUP'][key], key, 5, 50, 10, left_dup, right_dup)
		if key not in semi_result:
			semi_result[key] = list()
		for i in result:
			i.append("DUP")
			semi_result[key].append(i)
	
	for key in candidate["INV"]:
		candidate["INV"][key] = sorted(candidate["INV"][key], key = lambda x:x[:])
		try:
			left_inv = mkindex(candidate["l_INV"][key], 0)
		except:
			left_inv = dict()
		try:
			right_inv = mkindex(candidate["INV_r"][key], 1)
		except:
			right_inv = dict()
		result = intergrate_dup(candidate['INV'][key], key, 5, 50, 10, left_dup, right_dup)
		if key not in semi_result:
			semi_result[key] = list()
		for i in result:
			i.append("INV")
			semi_result[key].append(i)

	for key in semi_result:
		for i in semi_result[key]:
			print i
	# print candidate['DEL']
