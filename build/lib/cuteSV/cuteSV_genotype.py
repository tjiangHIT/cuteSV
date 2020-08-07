
from cuteSV.cuteSV_Description import Generation_VCF_header
from math import log10
import numpy as np

err = 0.1
prior = float(1/3)
Genotype = ["0/0", "0/1", "1/1"]

def log10sumexp(log10_probs):
	# Normalization of Genotype likelihoods
	m = max(log10_probs)
	return m + log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
	# Adjust the Genotype likelihoods
	log10_probs = np.array(log10_probs)
	lse = log10sumexp(log10_probs)
	return np.minimum(log10_probs - lse, 0.0)

def rescale_read_counts(c0, c1, max_allowed_reads=100):
	"""Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
	Total = c0 + c1
	if Total > max_allowed_reads:
		c0 = int(max_allowed_reads * float(c0/Total))
		c1 = max_allowed_reads - c0
	return c0, c1

def cal_GL(c0, c1):
	# Approximate adjustment of events with larger read depth
	c0, c1 = rescale_read_counts(c0, c1)
	# original genotype likelihood
	# ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*comb(c0+c1,c0)*(1-prior)/2)
	# ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*comb(c0+c1,c0)*(1-prior)/2)
	# ori_GL01 = np.float64(pow(0.5, c0+c1)*comb(c0+c1,c0)*prior)
	
	ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior)/2)
	ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*(1-prior)/2)
	ori_GL01 = np.float64(pow(0.5, c0+c1)*prior)

	# normalized genotype likelihood
	prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
	GL_P = [pow(10, i) for i in prob]
	PL = [int(np.around(-10*log10(i))) for i in GL_P]
	GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
	QUAL = abs(np.around(-10*log10(GL_P[0]), 1))

	return Genotype[prob.index(max(prob))], "%d,%d,%d"%(PL[0], PL[1], PL[2]), max(GQ), QUAL

def cal_CIPOS(std, num):
	pos = int(1.96 * std / num ** 0.5)
	return "-%d,%d"%(pos,pos)

def threshold_ref_count(num):
	if num <= 2:
		return 10*num
	elif 3 <= num <= 5:
		return 5*num 
	elif 6 <= num <= 15:
		return 4*num
	else:
		return 3*num

def count_coverage(chr, s, e, f, read_count, up_bound, itround):
	status = 0
	iteration = 0
	primary_num = 0
	for i in f.fetch(chr, s, e):
		iteration += 1
		if i.flag not in [0,16]:
			continue
		primary_num += 1
		if i.reference_start < s and i.reference_end > e:
			read_count.add(i.query_name)
			if len(read_count) >= up_bound:
				status = 1
				break
		if iteration >= itround:
			if float(primary_num/iteration) <= 0.2:
				status = 1
			else:
				status = -1
			break

	return status

def generate_output(args, semi_result, contigINFO, argv):
	
	'''
	Generation of VCF format file.
	VCF version: 4.2
	'''

	# genotype_trigger = TriggerGT[args.genotype]

	svid = dict()
	svid["INS"] = 0
	svid["DEL"] = 0
	svid["BND"] = 0
	svid["DUP"] = 0
	svid["INV"] = 0

	file = open(args.output, 'w')
	Generation_VCF_header(file, contigINFO, args.sample, argv)

	for i in semi_result:
		if i[1] in ["DEL", "INS"]:
			if i[1] == "INS":
				cal_end = int(i[2]) + 1
			else:
				cal_end = int(i[2]) + abs(int(float(i[3])))
			info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};RE={RE};RNAMES={RNAMES}".format(
				PRECISION = "IMPRECISE" if i[8] == "0/0" else "PRECISE", 
				SVTYPE = i[1], 
				SVLEN = i[3], 
				END = str(cal_end), 
				CIPOS = i[5], 
				CILEN = i[6], 
				RE = i[4],
				RNAMES = i[12])
			if i[1] =="DEL":
				info_list += ";STRAND=+-"
			if i[11] == ".":
				filter_lable = "PASS"
			else:
				filter_lable = "PASS" if float(i[11]) >= 5.0 else "q5"
			file.write("{CHR}\t{POS}\t{ID}\tN\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
				CHR = i[0], 
				POS = i[2], 
				ID = "cuteSV.%s.%d"%(i[1], svid[i[1]]),
				ALT = "<%s>"%(i[1]), 
				INFO = info_list, 
				FORMAT = "GT:DR:DV:PL:GQ", 
				GT = i[8],
				DR = i[7],
				RE = i[4],
				PL = i[9],
				GQ = i[10],
				QUAL = i[11],
				PASS = filter_lable))
			svid[i[1]] += 1
		elif i[1] == "DUP":
			cal_end = int(i[2]) + abs(int(float(i[3])))
			info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};STRAND=-+;RNAMES={RNAMES}".format(
				PRECISION = "IMPRECISE" if i[6] == "0/0" else "PRECISE", 
				SVTYPE = i[1], 
				SVLEN = i[3], 
				END = str(cal_end), 
				RE = i[4],
				RNAMES = i[10])
			if i[9] == ".":
				filter_lable = "PASS"
			else:
				filter_lable = "PASS" if float(i[9]) >= 5.0 else "q5"
			file.write("{CHR}\t{POS}\t{ID}\tN\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
				CHR = i[0], 
				POS = i[2], 
				ID = "cuteSV.%s.%d"%(i[1], svid[i[1]]),
				ALT = "<%s>"%(i[1]), 
				INFO = info_list, 
				FORMAT = "GT:DR:DV:PL:GQ", 
				GT = i[6],
				DR = i[5],
				RE = i[4],
				PL = i[7],
				GQ = i[8],
				QUAL = i[9],
				PASS = filter_lable))
			svid[i[1]] += 1
		elif i[1] == "INV":
			cal_end = int(i[2]) + abs(int(float(i[3])))
			info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};STRAND={STRAND};RNAMES={RNAMES}".format(
				PRECISION = "IMPRECISE" if i[6] == "0/0" else "PRECISE", 
				SVTYPE = i[1], 
				SVLEN = i[3], 
				END = str(cal_end), 
				RE = i[4],
				STRAND = i[7],
				RNAMES = i[11])
			if i[10] == ".":
				filter_lable = "PASS"
			else:
				filter_lable = "PASS" if float(i[10]) >= 5.0 else "q5"
			file.write("{CHR}\t{POS}\t{ID}\tN\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
				CHR = i[0], 
				POS = i[2], 
				ID = "cuteSV.%s.%d"%(i[1], svid[i[1]]),
				ALT = "<%s>"%(i[1]), 
				INFO = info_list, 
				FORMAT = "GT:DR:DV:PL:GQ", 
				GT = i[6],
				DR = i[5],
				RE = i[4],
				PL = i[8],
				GQ = i[9],
				QUAL = i[10],
				PASS = filter_lable))
			svid[i[1]] += 1
		else:
			# BND
			info_list = "{PRECISION};SVTYPE={SVTYPE};CHR2={CHR2};END={END};RE={RE};RNAMES={RNAMES}".format(
				PRECISION = "IMPRECISE" if i[7] == "0/0" else "PRECISE", 
				SVTYPE = "BND", 
				CHR2 = i[3], 
				END = i[4], 
				RE = i[5],
				RNAMES = i[11])
			if i[10] == ".":
				filter_lable = "PASS"
			else:
				filter_lable = "PASS" if float(i[10]) >= 5.0 else "q5"
			file.write("{CHR}\t{POS}\t{ID}\tN\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
				CHR = i[0], 
				POS = i[2], 
				ID = "cuteSV.%s.%d"%("BND", svid["BND"]), 
				ALT = i[1], 
				INFO = info_list, 
				FORMAT = "GT:DR:DV:PL:GQ", 
				GT = i[7],
				DR = i[6],
				RE = i[5],
				PL = i[8],
				GQ = i[9],
				QUAL = i[10],
				PASS = filter_lable))
			svid["BND"] += 1
