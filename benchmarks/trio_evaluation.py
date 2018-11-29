import sys
import argparse
import logging
import time

def load_callset_cuteSV(path, filter=10):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0] not in callset:
			callset[seq[0]] = dict()
		if seq[1] == 'TRA':
			if seq[1] not in callset[seq[0]]:
				callset[seq[0]][seq[1]] = dict()
			if seq[3] not in callset[seq[0]][seq[1]]:
				callset[seq[0]][seq[1]][seq[3]] = list()
			if int(seq[5]) >= filter:
				callset[seq[0]][seq[1]][seq[3]].append([int(seq[2]), int(seq[4]), 0])
		else:
			if seq[1] not in callset[seq[0]]:
				callset[seq[0]][seq[1]] = list()
			if int(seq[4]) >= filter:
				callset[seq[0]][seq[1]].append([int(seq[2]), int(seq[2])+int(seq[3]), 0])
	file.close()
	return callset

def load_callset_sniffles(path, filter=10):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		if seq[0] not in callset:
			callset[seq[0]] = dict()
		pos_2 = int(seq[7].split(';')[3].split('=')[1])
		if int(seq[-1].split(':')[-1]) < filter:
			continue
		if seq[4][1:-1] == 'TRA':
			if seq[4][1:-1] not in callset[seq[0]]:
				callset[seq[0]][seq[4][1:-1]] = dict()
				chr_2 = seq[7].split(';')[2].split('=')[1]
			if chr_2 not in callset[seq[0]][seq[4][1:-1]]:
				callset[seq[0]][seq[4][1:-1]][chr_2] = list()
			callset[seq[0]][seq[4][1:-1]][chr_2].append([int(seq[1]), pos_2, 0])
		else:
			if seq[4][1:-1] not in callset[seq[0]]:
				callset[seq[0]][seq[4][1:-1]] = list()
			callset[seq[0]][seq[4][1:-1]].append([int(seq[1]), pos_2, 0])
	file.close()
	return callset

def load_callset_pbsv(path, filter=10):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		if seq[0] not in callset:
			callset[seq[0]] = dict()

		if seq[-1].split(":")[0] == "0/1":
			if int(seq[-1].split(",")[1].split(":")[0]) < filter:
				continue
		else:
			if int(seq[-1].split(":")[-1]) < filter:
				continue

		svtype = seq[2].split('.')[1]
		if svtype == 'BND':
			svtype = "TRA"
			if svtype not in callset[seq[0]]:
				callset[seq[0]][svtype] = dict()
				if len(seq[4].split(']')) > 1:
					chr_2 = seq[4].split(']')[1].split(':')[0]
					pos_2 = int(seq[4].split(']')[1].split(':')[1])
				else:
					chr_2 = seq[4].split('[')[1].split(':')[0]
					pos_2 = int(seq[4].split('[')[1].split(':')[1])
			if chr_2 not in callset[seq[0]][svtype]:
				callset[seq[0]][svtype][chr_2] = list()
			callset[seq[0]][svtype][chr_2].append([int(seq[1]), pos_2, 0])
		else:
			if svtype not in callset[seq[0]]:
				callset[seq[0]][svtype] = list()
			# print seq[7]
			try:
				svlen = int(seq[7].split(';')[2].split('=')[1])
				callset[seq[0]][svtype].append([int(seq[1]), int(seq[1])+abs(svlen), 0])
			except:
				pos_2 = int(seq[7].split(';')[1].split('=')[1])
				callset[seq[0]][svtype].append([int(seq[1]), pos_2, 0])
	file.close()
	return callset

def cal_overlap(a_l, a_r, b_l, b_r, bias):
	if a_l < b_l:
		if b_l >= a_r:
			return 0
		else:
			if b_r <= a_r:
				if (a_r - a_l)*bias <= b_r - b_l:
					return 1
				else:
					return 0
			else:
				overlap = a_r - b_l
				if (a_r - a_l)*bias <= overlap and (b_r - b_l)*bias <= overlap:
					return 1
				else:
					return 0
	else:
		if a_l >= b_r:
			return 0
		else:
			if a_r <= b_r:
				if (b_r - b_l)*bias <= a_r - a_l:
					return 1
				else:
					return 0
			else:
				overlap = b_r - a_l
				if (a_r - a_l)*bias <= overlap and (b_r - b_l)*bias <= overlap:
					return 1
				else:
					return 0

def cal_overlap_tra(a_l, a_r, b_l, b_r, bias):
	if abs(a_l-b_l) <= bias and abs(a_r-b_r) <= bias:
		return 1
	else:
		return 0

def eva_record(call_A, call_B, bias, offect):
	for chr in call_A:
		if chr not in call_B:
			continue
		for svtype in call_A[chr]:
			if svtype not in call_B[chr]:
				continue
			if svtype == "TRA":
				for chr_2 in call_A[chr][svtype]:
					if chr_2 not in call_B[chr][svtype]:
						continue
					for i in xrange(len(call_A[chr][svtype][chr_2])):
						for j in xrange(len(call_B[chr][svtype][chr_2])):
							judgement = cal_overlap_tra(call_A[chr][svtype][chr_2][i][0], call_A[chr][svtype][chr_2][i][1], call_B[chr][svtype][chr_2][j][0], call_B[chr][svtype][chr_2][j][1], offect)
							if judgement == 1:
								call_A[chr][svtype][chr_2][i][2] = 1
								call_B[chr][svtype][chr_2][j][2] = 1
								# exist a bug
			elif svtype == "INS":
				for i in xrange(len(call_A[chr][svtype])):
					for j in xrange(len(call_B[chr][svtype])):
						if abs(call_A[chr][svtype][i][0]-call_B[chr][svtype][j][0]) <= offect:
							call_A[chr][svtype][i][2] = 1
							call_B[chr][svtype][j][2] = 1
			else:
				for i in xrange(len(call_A[chr][svtype])):
					for j in xrange(len(call_B[chr][svtype])):
						judgement = cal_overlap(call_A[chr][svtype][i][0], call_A[chr][svtype][i][1], call_B[chr][svtype][j][0], call_B[chr][svtype][j][1], bias)
						if judgement == 1:
							call_A[chr][svtype][i][2] = 1
							call_B[chr][svtype][j][2] = 1

def statistics_true_possitive(callset, SVTYPE="ALL"):
	record = 0
	true_record = 0
	if SVTYPE == "ALL":
		for chr in callset:
			for svtype in callset[chr]:
				if svtype == "TRA":
					for chr_2 in callset[chr][svtype]:
						for res in callset[chr][svtype][chr_2]:
							record += 1
							if res[2] == 1:
								true_record += 1
				else:
					for res in callset[chr][svtype]:
						record += 1
						if res[2] == 1:
							true_record += 1
	else:
		for chr in callset:
			if SVTYPE in callset[chr]:
				if SVTYPE == "TRA":
					for chr_2 in callset[chr][SVTYPE]:
						for res in callset[chr][SVTYPE][chr_2]:
							record += 1
							if res[2] == 1:
								true_record += 1
				else:
					for res in callset[chr][SVTYPE]:
						record += 1
						if res[2] == 1:
							true_record += 1
	return record, true_record

def statistics_true_possitive_print(callset, SVTYPE="ALL"):
	record = 0
	true_record = 0
	if SVTYPE == "ALL":
		for chr in callset:
			for svtype in callset[chr]:
				if svtype == "TRA":
					for chr_2 in callset[chr][svtype]:
						for res in callset[chr][svtype][chr_2]:
							record += 1
							if res[2] == 1:
								true_record += 1
				else:
					for res in callset[chr][svtype]:
						record += 1
						if res[2] == 1:
							true_record += 1
	else:
		for chr in callset:
			if SVTYPE in callset[chr]:
				if SVTYPE == "TRA":
					for chr_2 in callset[chr][SVTYPE]:
						for res in callset[chr][SVTYPE][chr_2]:
							record += 1
							if res[2] == 1:
								true_record += 1
				else:
					for res in callset[chr][SVTYPE]:
						record += 1
						if res[2] == 1:
							true_record += 1
						else:
							print("%s\t%s\t%d\t%d"%(chr, SVTYPE, res[0], res[1]))
	return record, true_record

def main_ctrl(args):
	logging.info("Load SV callset of cuteSV.")
	if args.CALLER == "cuteSV":
		# '''
		call_child = load_callset_cuteSV(args.F1, int(args.F1_c))
		call_father = load_callset_cuteSV(args.MP, int(args.MP_c))
		call_mother = load_callset_cuteSV(args.FP, int(args.FP_c))
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "DUP", "INS", "INV", "TRA"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			sensitivity = 1.0*(mother_tr+father_tr)/(mother_r+father_r)
			accuracy = 1.0*child_tr/child_r
			logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, child_tr, child_r, accuracy))
			logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, mother_tr+father_tr, mother_r+father_r, sensitivity))
			# logging.info("%s: F1-score is %.5f."%(i, 2*sensitivity*accuracy/(sensitivity+accuracy)))
			# logging.info("%s: G-score is %.5f."%(i, (sensitivity*accuracy)**0.05))
			# logging.info("%s: ADI is %.5f."%(i, 1-accuracy))
		# statistics_true_possitive_print(call_child, "DUP")
		# '''
		# call_father = load_callset_cuteSV(args.MP)
		# call_mother = load_callset_cuteSV(args.FP)
		# for i in xrange(7):
		# 	for j in xrange(9):
		# 		call_child = load_callset_cuteSV(args.F1, 10+2*i)
		# 		eva_record(call_child, call_father, 0.1*(j+1), args.offect)
		# 		eva_record(call_child, call_mother, 0.1*(j+1), args.offect)
		# 		child_r, child_tr = statistics_true_possitive(call_child, "DEL")
		# 		ADI = 1.0-1.0*child_tr/child_r
		# 		logging.info("ADI rate is %.5f for coverage %d and bias %.2f."%(ADI, 10+2*i, 0.1*(j+1)))

	elif args.CALLER == "sniffles":
		call_child = load_callset_sniffles(args.F1, int(args.F1_c))
		call_father = load_callset_sniffles(args.MP, int(args.MP_c))
		call_mother = load_callset_sniffles(args.FP, int(args.FP_c))
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "DUP", "INS", "INV", "TRA"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			sensitivity = 1.0*(mother_tr+father_tr)/(mother_r+father_r)
			accuracy = 1.0*child_tr/child_r
			logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, child_tr, child_r, accuracy))
			logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, mother_tr+father_tr, mother_r+father_r, sensitivity))
			logging.info("%s: F1-score is %.5f."%(i, 2*sensitivity*accuracy/(sensitivity+accuracy)))
			logging.info("%s: G-score is %.5f."%(i, (sensitivity*accuracy)**0.05))
	elif args.CALLER == "pbsv":
		call_child = load_callset_pbsv(args.F1, int(args.F1_c))
		call_father = load_callset_pbsv(args.MP, int(args.MP_c))
		call_mother = load_callset_pbsv(args.FP, int(args.FP_c))
		logging.info("Evaluate accuracy and sensitivity.")
		eva_record(call_child, call_father, args.bias, args.offect)
		eva_record(call_child, call_mother, args.bias, args.offect)
		svtype = ["DEL", "INS", "INV", "TRA"]
		for i in svtype:
			child_r, child_tr = statistics_true_possitive(call_child, i)
			father_r, father_tr = statistics_true_possitive(call_father, i)
			mother_r, mother_tr = statistics_true_possitive(call_mother, i)
			try:
				logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, child_tr, child_r, 1.0*child_tr/child_r))
			except:
				logging.info("%s: Accuracy rate is %d/%d."%(i, child_tr, child_r))
			try:
				logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, mother_tr+father_tr, mother_r+father_r, 1.0*(mother_tr+father_tr)/(mother_r+father_r)))
			except:
				logging.info("%s: Sensitivity rate is %d/%d."%(i, mother_tr+father_tr, mother_r+father_r))
	# logging.info("Evaluate accuracy and sensitivity.")
	# eva_record(call_child, call_father, args.bias, args.offect)
	# eva_record(call_child, call_mother, args.bias, args.offect)
	child_r, child_tr = statistics_true_possitive(call_child)
	father_r, father_tr = statistics_true_possitive(call_father)
	mother_r, mother_tr = statistics_true_possitive(call_mother)
	sensitivity = 1.0*(mother_tr+father_tr)/(mother_r+father_r)
	accuracy = 1.0*child_tr/child_r
	logging.info("Accuracy rate is %d/%d(%.5f)."%(child_tr, child_r, accuracy))
	logging.info("Sensitivity rate is %d/%d(%.5f)."%(mother_tr+father_tr, mother_r+father_r, sensitivity))
	logging.info("F1-score is %.5f."%(2*sensitivity*accuracy/(sensitivity+accuracy)))
	# logging.info("G-score is %.5f."%((sensitivity*accuracy)**0.05))

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
	parser.add_argument("MP_c", type=str, help="Coverage of male parent callsets")
	parser.add_argument('FP_c', type=str, help = "Coverage of female parent callsets")
	parser.add_argument('F1_c', type=str, help = "Coverage of offspring callsets")
	parser.add_argument('-b', '--bias', help = "Bias of overlaping.[%(default)s]", default = 0.7, type = float)
	parser.add_argument('-o', '--offect', help = "Offect of translocation overlaping.[%(default)s]", default = 50, type = int)
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
