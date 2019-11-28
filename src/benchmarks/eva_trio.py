import sys
import argparse
import logging
import time

def load_high_conf_bed(path):
	confbed = dict()
	file = open(path,'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		istart = int(seq[1])
		iend = int(seq[2])
		if chr not in confbed:
			confbed[chr] = list()
		confbed[chr].append([istart, iend])
	file.close()
	return confbed

def judge_bed(a, b, clist):
	for i in clist:
		if a <= i[0] and b >= i[0]:
			return 1
		if i[0] >= a and a <= i[1]:
			return 1
	return 0

def parse_BND(tag):
	if len(tag.split('[')) == 1:
		return tag.split(']')[1].split(':')[0], int(tag.split(']')[1].split(':')[1])
	else:
		return tag.split('[')[1].split(':')[0], int(tag.split('[')[1].split(':')[1])

def pase_info(seq):
	info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0, "CHR2": ''}
	for i in seq.split(';'):
		if i.split('=')[0] in ["SVLEN", "END", "RE"]:
			try:
				info[i.split('=')[0]] = abs(int(i.split('=')[1]))
			except:
				pass
		if i.split('=')[0] in ["SVTYPE", "CHR2"]:
			info[i.split('=')[0]] = i.split('=')[1]
	return info

def pase_info_2(seq, seq2):
	info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "SUPPORT": 0, "CHR2": ''}
	for i in seq.split(';'):
		if i.split('=')[0] in ["SVLEN", "END", "SUPPORT"]:
			try:
				info[i.split('=')[0]] = abs(int(i.split('=')[1]))
			except:
				pass
		if i.split('=')[0] in ["SVTYPE"]:
			info[i.split('=')[0]] = i.split('=')[1]
			if i.split('=')[1] == 'BND':
				if seq2[0] == 'N':
					info['CHR2'] = seq[2].split(':')[0][2:]
				else:
					info['CHR2'] = seq[2].split(':')[0][1:]
	return info

def load_callset_cuteSV(path, filter, confbed):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		ALT = seq[2][7:10]
		if ALT == "DUP":
			ALT = "INS"
		if ALT not in ["INS", "INV", "DEL", "DUP", "BND"]:
			continue

		info = pase_info(seq[7])
		if len(confbed) > 0:
			if chr not in confbed:
				continue
			if judge_bed(pos, info["END"], confbed[chr]) == 0:
				continue

		if ALT not in callset:
			callset[ALT] = dict()
		if chr not in callset[ALT]:
			if ALT == 'BND':
				callset[ALT][chr] = dict()
			else:
				callset[ALT][chr] = list()

		if ALT == "BND":
			if info["CHR2"] not in callset[ALT][chr]:
				callset[ALT][chr][info["CHR2"]] = list()
			if info["RE"] >= filter:
				callset[ALT][chr][info["CHR2"]].append([pos, info["END"], 0])
		else:
			if info["SVLEN"] >= 50 and info["RE"] >= filter:
				callset[ALT][chr].append([pos, info["SVLEN"], info["END"], 0])
	file.close()
	return callset


def load_callset_sniffles(path, filter, confbed):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		info = pase_info(seq[7])

		if info["SVTYPE"] == "DUP":
			info["SVTYPE"] = "INS"
		svtype = info["SVTYPE"]

		if svtype == "BND":
			chr_2, pos_2 = parse_BND(seq[4])
			if len(confbed) > 0:
				if chr not in confbed:
					continue
				if judge_bed(pos, pos_2, confbed[chr]) == 0:
					continue

			if svtype not in callset:
				callset[svtype] = dict()
			if chr not in callset[svtype]:
				callset[svtype][chr] = dict()
			if chr_2 not in callset[svtype][chr]:
				callset[svtype][chr][chr_2] = list()

			if info["RE"] >= filter:
				callset[svtype][chr][chr_2].append([pos, pos_2, 0])

		else:
			END = info["END"]
			if len(confbed) > 0:
				if chr not in confbed:
					continue
				if judge_bed(pos, END, confbed[chr]) == 0:
					continue

			if svtype not in callset:
				callset[svtype] = dict()
			if chr not in callset[svtype]:
				callset[svtype][chr] = list()
			if info["SVLEN"] >= 50 and info["RE"] >= filter:
				callset[svtype][chr].append([pos, info["SVLEN"], END, 0])

	file.close()
	return callset

def load_callset_svim(path, filter, confbed):
	base_call = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split("\t")
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		info = pase_info_2(seq[7], seq[4])

		if len(confbed) > 0:
			if chr not in confbed:
				continue
			if judge_bed(pos, info["END"], confbed[chr]) == 0:
				continue

		svtype = info["SVTYPE"]
		if svtype == "BND":
			chr_2, pos_2 = parse_BND(seq[4])
			if len(confbed) > 0:
				if chr not in confbed:
					continue
				if judge_bed(pos, pos_2, confbed[chr]) == 0:
					continue

			if svtype not in base_call:
				base_call[svtype] = dict()
			if chr not in base_call[svtype]:
				base_call[svtype][chr] = dict()
			if chr_2 not in base_call[svtype][chr]:
				base_call[svtype][chr][chr_2] = list()

			if info["SUPPORT"] >= filter:
				base_call[svtype][chr][chr_2].append([pos, pos_2, 0])

		elif info["SVTYPE"] == "INV":
			if info["SVTYPE"] not in base_call:
				base_call[info["SVTYPE"]] = dict()
			if chr not in base_call[info["SVTYPE"]]:
				base_call[info["SVTYPE"]][chr] = list()
			if info["END"] - pos + 1 >= 50:
				base_call[info["SVTYPE"]][chr].append([pos, info["END"] - pos + 1, info["END"], 0])
		else:
			if info["SVTYPE"] not in base_call:
				base_call[info["SVTYPE"]] = dict()
			if chr not in base_call[info["SVTYPE"]]:
				base_call[info["SVTYPE"]][chr] = list()
			if info["SVLEN"] >= 50:
				base_call[info["SVTYPE"]][chr].append([pos, info["SVLEN"], info["END"], 0])
	file.close()
	return base_call

def load_callset_pbsv(path, filter, confbed):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue

		chr = seq[0]
		pos = int(seq[1])
		info = pase_info(seq[7])

		if info["SVTYPE"] not in ["INS", "INV", "DEL", "DUP", "BND"]:
			continue

		try:
			readcount = int(seq[-1].split(':')[1].split(',')[1])
		except:
			continue
		# if readcount < filter:
		# 	continue

		if info["SVTYPE"] == 'BND' and readcount < filter:
			continue

		if info["SVTYPE"] == "DUP":
			info["SVTYPE"] = "INS"

		if info["SVTYPE"] == "BND":
			chr_2, pos_2 = parse_BND(seq[4])
			if len(confbed) > 0:
				if chr not in confbed:
					continue
				if judge_bed(pos, pos_2, confbed[chr]) == 0:
					continue
			if info["SVTYPE"] not in callset:
				callset[info["SVTYPE"]] = dict()
			if chr not in callset[info["SVTYPE"]]:
				callset[info["SVTYPE"]][chr] = dict()
			if chr_2 not in callset[info["SVTYPE"]][chr]:
				callset[info["SVTYPE"]][chr][chr_2] = list()

			callset[info["SVTYPE"]][chr][chr_2].append([pos, pos_2, 0])

		else:
			if info["SVTYPE"] not in callset:
				callset[info["SVTYPE"]] = dict()
			if chr not in callset[info["SVTYPE"]]:
				callset[info["SVTYPE"]][chr] = list()
			if info["SVTYPE"] == "INV":
				info["SVLEN"] = info["END"] - pos + 1
			if len(confbed) > 0:
				if chr not in confbed:
					continue
				if judge_bed(pos, info["END"], confbed[chr]) == 0:
					continue
			if info["SVLEN"] >= 50:
				callset[info["SVTYPE"]][chr].append([pos, info["SVLEN"], info["END"], 0])

	file.close()
	return callset


def eva_record(call_A, call_B, bias, offect):
	for svtype in call_A:
		if svtype not in call_B:
			continue
		else:
			if svtype == "BND":
				for chr1 in call_A[svtype]:
					if chr1 not in call_B[svtype]:
						continue
					else:
						for chr2 in call_A[svtype][chr1]:
							if chr2 not in call_B[svtype][chr1]:
								continue
							else:
								for i in call_A[svtype][chr1][chr2]:
									for j in call_B[svtype][chr1][chr2]:
										if abs(i[0] - j[0]) <= offect and abs(i[1] - j[1]) <= offect:
											i[2] = 1
											j[2] = 1
										else:
											pass
			else:
				for chr in call_A[svtype]:
					if chr not in call_B[svtype]:
						continue
					else:
						for i in call_A[svtype][chr]:
							for j in call_B[svtype][chr]:
								# if min(i[2], j[2]) >= max(i[0], j[0]):
								if i[0] - offect <= j[0] <= i[2] + offect or i[0] - offect <= j[2] <= i[2] + offect or j[0] - offect <= i[0] <= j[2] + offect:
									if min(i[1], j[1])*1.0/max(i[1], j[1]) >= bias:
										i[3] = 1
										j[3] = 1
									else:
										pass


def statistics_true_possitive(callset, SVTYPE="ALL"):
	record = 0
	true_record = 0
	if SVTYPE == "ALL":
		for svtype in callset:
			if svtype == "BND":
				for chr1 in callset[svtype]:
					for chr2 in callset[svtype][chr1]:
						for res in callset[svtype][chr1][chr2]:
							record += 1
							if res[2] == 1:
								true_record += 1
			else:
				for chr in callset[svtype]:
					for res in callset[svtype][chr]:
						record += 1
						if res[3] == 1:
							true_record += 1
	else:
		if SVTYPE == "BND":
			for chr1 in callset[SVTYPE]:
				for chr2 in callset[SVTYPE][chr1]:
					for res in callset[SVTYPE][chr1][chr2]:
						record += 1
						if res[2] == 1:
							true_record += 1
		else:
			for chr in callset[SVTYPE]:
				for res in callset[SVTYPE][chr]:
					record += 1
					if res[3] == 1:
						true_record += 1
	return record, true_record


def main_ctrl(args):	
	if args.include == "NULL":
		confbed = dict()
	else:
		logging.info("Load high confidence region.")
		confbed = load_high_conf_bed(args.include)
	logging.info("Load SV callset of selected caller.")
	if args.CALLER == "cuteSV":
		call_child = load_callset_cuteSV(args.F1, int(args.F1_c), confbed)
		call_father = load_callset_cuteSV(args.MP, int(args.MP_c), confbed)
		call_mother = load_callset_cuteSV(args.FP, int(args.FP_c), confbed)
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "INS", "INV", "BND"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			accuracy = 1.0*child_tr/child_r
			logging.info("%s\tChild:%d\tFather:%d\tMother:%d."%(i, child_r, father_r, mother_r))
			logging.info("%s: ADI is %.2f."%(i, 100 - 100.0 * accuracy))
	elif args.CALLER == "sniffles":
		call_child = load_callset_sniffles(args.F1, int(args.F1_c), confbed)
		# for svtype in call_child:
		# 	if svtype == "BND":
		# 		for chr in call_child[svtype]:
		# 			for chr_2 in call_child[svtype][chr]:
		# 				for i in call_child[svtype][chr][chr_2]:
		# 					print(svtype, chr, chr_2, i)
		# 	else:
		# 		for chr in call_child[svtype]:
		# 			for i in call_child[svtype][chr]:
		# 				print(svtype, chr, i)
		call_father = load_callset_sniffles(args.MP, int(args.MP_c), confbed)
		call_mother = load_callset_sniffles(args.FP, int(args.FP_c), confbed)
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "INS", "INV", "BND"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			accuracy = 1.0*child_tr/child_r
			logging.info("%s\tChild:%d\tFather:%d\tMother:%d."%(i, child_r, father_r, mother_r))
			logging.info("%s: ADI is %.2f."%(i, 100 - 100.0 * accuracy))
	elif args.CALLER == "pbsv":
		call_child = load_callset_pbsv(args.F1, int(args.F1_c), confbed)
		call_father = load_callset_pbsv(args.MP, int(args.MP_c), confbed)
		call_mother = load_callset_pbsv(args.FP, int(args.FP_c), confbed)
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "INS", "INV", "BND"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			accuracy = 1.0*child_tr/child_r
			logging.info("%s\tChild:%d\tFather:%d\tMother:%d."%(i, child_r, father_r, mother_r))
			logging.info("%s: ADI is %.2f."%(i, 100 - 100.0 * accuracy))
	elif args.CALLER == 'svim':
		call_child = load_callset_svim(args.F1, int(args.F1_c), confbed)
		call_father = load_callset_svim(args.MP, int(args.MP_c), confbed)
		call_mother = load_callset_svim(args.FP, int(args.FP_c), confbed)
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "INS", "INV", "BND"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			accuracy = 1.0*child_tr/child_r
			logging.info("%s\tChild:%d\tFather:%d\tMother:%d."%(i, child_r, father_r, mother_r))
			logging.info("%s: ADI is %.2f."%(i, 100 - 100.0 * accuracy))
	child_r, child_tr = statistics_true_possitive(call_child)
	father_r, father_tr = statistics_true_possitive(call_father)
	mother_r, mother_tr = statistics_true_possitive(call_mother)
	accuracy = 1.0*child_tr/child_r
	logging.info("Child:%d\tFather:%d\tMother:%d."%(child_r, father_r, mother_r))
	logging.info("ADI is %.2f."%(100 - 100.0 * accuracy))

def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Evaluate SV callset generated by cuteSV
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="Trio_eval", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("CALLER", type=str, help="Choose a caller")
	parser.add_argument("MP", type=str, help="Male parent callsets")
	parser.add_argument('FP', type=str, help = "Female parent callsets")
	parser.add_argument('F1', type=str, help = "Offspring callsets")
	parser.add_argument("MP_c", type=int, help="Coverage of male parent callsets")
	parser.add_argument('FP_c', type=int, help = "Coverage of female parent callsets")
	parser.add_argument('F1_c', type=int, help = "Coverage of offspring callsets")
	parser.add_argument('-b', '--bias', help = "Bias of overlaping.[%(default)s]", default = 0.7, type = float)
	parser.add_argument('-o', '--offect', help = "Offect of translocation overlaping.[%(default)s]", default = 1000, type = int)
	parser.add_argument('-i', '--include', help = "Include high confident region (bed format).",default = "NULL", type = str)
	# parser.add_argument('-s', '--min_support', help = "Minimum number of reads that support a SV to be reported.[%(default)s]", default = 10, type = int)
	args = parser.parse_args(argv)
	return args

def setupLogging(debug=False):
	logLevel = logging.DEBUG if debug else logging.INFO
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
	logging.info("Running %s" % " ".join(sys.argv))

if __name__ == '__main__':
	main(sys.argv[1:])
