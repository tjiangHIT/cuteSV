import sys
import argparse
import logging
import time

def transfer_chrID(seq):
	if len(seq.split('_')) > 1:
		if seq.split('_')[1][:2] == 'gl':
			return "GL%s.1"%(seq.split('_')[1][2:])
	else:
		return seq[3:]


def load_golden_answer(path):
	Ans = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		if seq[6] == "lt3":
			continue
		chr = transfer_chrID(seq[0])
		if chr not in Ans:
			Ans[chr] = dict()
		pos = int(seq[1])
		svtype = seq[4][1:-1]
		if svtype not in Ans[chr]:
			Ans[chr][svtype] = list()
		if seq[7].split(';')[3][0] == "S":
			try:
				svlen = int(seq[7].split(';')[3].split('=')[1])
			except:
				svlen = int(seq[7].split(';')[2].split('=')[1]) - pos
		else:
			svlen = int(seq[7].split(';')[2].split('=')[1]) - pos

		if svtype == "INS":
			Ans[chr][svtype].append([pos, abs(svlen), 0])
		else:
			Ans[chr][svtype].append([pos, pos+abs(svlen), 0])
	file.close()
	return Ans

def load_callset_cuteSV(path, filter=10):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		chr = seq[0]
		if chr not in callset:
			callset[chr] = dict()
		svtype = seq[1]
		if svtype not in ["DEL", "INS"]:
			continue
		if svtype not in callset[chr]:
			callset[chr][svtype] = list()
		pos = int(seq[2])
		svlen = int(seq[3])
		if int(seq[-1]) < filter:
			continue
		if svtype == "INS":
			callset[chr][svtype].append([pos, svlen, 0])
		else:
			callset[chr][svtype].append([pos, pos+svlen, 0])
	file.close()
	return callset

def load_callset_sniffles(path, filter=10):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		if chr not in callset:
			callset[chr] = dict()
		pos = int(seq[1])
		svtype = seq[4][1:-1]
		if svtype not in ["DEL", "INS"]:
			continue
		if svtype not in callset[chr]:
			callset[chr][svtype] = list()

		svlen = abs(int(seq[7].split(";")[10].split('=')[1]))
		# print chr, pos, svtype, svlen
		if int(seq[-1].split(':')[-1]) < filter:
			continue

		if svtype == "INS":
			callset[chr][svtype].append([pos, svlen, 0])
		else:
			callset[chr][svtype].append([pos, pos+abs(svlen), 0])
	file.close()
	return callset

def load_callset_pbsv(path, filter=10):
	callset = dict()
	file = open(path, 'r')
	for line in file:
		seq = line.strip('\n').split('\t')
		if seq[0][0] == '#':
			continue
		chr = seq[0]
		if chr not in callset:
			callset[chr] = dict()

		if seq[-1].split(":")[0] == "0/1":
			if int(seq[-1].split(",")[1].split(":")[0]) < filter:
				continue
		else:
			if int(seq[-1].split(":")[-1]) < filter:
				continue

		svtype = seq[2].split('.')[1]

		if svtype not in ["INS", "DEL"]:
			continue
		if svtype not in callset[chr]:
			callset[chr][svtype] = list()
		pos = int(seq[1])
		svlen = int(seq[7].split(';')[2].split('=')[1])
		if svtype == "INS":
			callset[chr][svtype].append([pos, svlen, 0])
		else:
			callset[chr][svtype].append([pos, pos+abs(svlen), 0])
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

def eve_record(Ans, Callset, bias, offect):
	for chr in Callset:
		if chr not in Ans:
			continue
		else:
			if "DEL" in Callset[chr] and "DEL" in Ans[chr]:
				for i in xrange(len(Callset[chr]["DEL"])):
					for j in xrange(len(Ans[chr]["DEL"])):
						judgement = cal_overlap(Callset[chr]["DEL"][i][0], Callset[chr]["DEL"][i][1], Ans[chr]["DEL"][j][0], Ans[chr]["DEL"][j][1], bias)
						if judgement == 1:
							Callset[chr]["DEL"][i][2] = 1
							Ans[chr]["DEL"][j][2] = 1
			if "INS" in Callset[chr] and "INS" in Ans[chr]:
				for i in xrange(len(Callset[chr]["INS"])):
					for j in xrange(len(Ans[chr]["INS"])):
						if abs(Callset[chr]["INS"][i][0]-Ans[chr]["INS"][j][0]) <= offect:
							Callset[chr]["INS"][i][2] = 1
							Ans[chr]["INS"][j][2] = 1

def statistics_true_possitive(callset, SVTYPE="ALL"):
	record = 0
	true_record = 0
	if SVTYPE == "ALL":
		for chr in callset:
			for svtype in callset[chr]:
				for res in callset[chr][svtype]:
					record += 1
					if res[2] == 1:
						true_record += 1
	else:
		for chr in callset:
			if SVTYPE in callset[chr]:
				for res in callset[chr][SVTYPE]:
					record += 1
					if res[2] == 1:
						true_record += 1
					else:
						print chr, SVTYPE, res[0], res[1]
	return record, true_record

def main_ctrl(args):
	logging.info("Load SV Golden Answers.")
	Ans = load_golden_answer(args.ANS)
	# for chr in Ans:
	# 	for sv in Ans[chr]:
	# 		print chr, sv, len(Ans[chr][sv])
	if args.CALLER == "sniffles":
		Callset = load_callset_sniffles(args.CALL, int(args.CALL_c))
		# for chr in Callset:
		# 	for sv in Callset[chr]:
		# 		print chr, sv, len(Callset[chr][sv])
		eve_record(Ans, Callset, args.bias, args.offect)
		svtype = ["DEL", "INS"]
		for i in svtype:
			caller_r, caller_tr = statistics_true_possitive(Callset, i)
			answer_r, answer_tr = statistics_true_possitive(Ans, i)
			sensitivity = 1.0*answer_tr/answer_r
			accuracy = 1.0*caller_tr/caller_r
			logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, caller_tr, caller_r, accuracy))
			logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, answer_tr, answer_r, sensitivity))
	elif args.CALLER == "cuteSV":
		# pass
		Callset = load_callset_cuteSV(args.CALL, int(args.CALL_c))
		eve_record(Ans, Callset, args.bias, args.offect)
		eve_record(Ans, Callset, args.bias, args.offect)
		svtype = ["DEL", "INS"]
		for i in svtype:
			caller_r, caller_tr = statistics_true_possitive(Callset, i)
			answer_r, answer_tr = statistics_true_possitive(Ans, i)
			sensitivity = 1.0*answer_tr/answer_r
			accuracy = 1.0*caller_tr/caller_r
			logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, caller_tr, caller_r, accuracy))
			logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, answer_tr, answer_r, sensitivity))
	elif args.CALLER == "pbsv":
		# pass
		Callset = load_callset_pbsv(args.CALL, int(args.CALL_c))
		eve_record(Ans, Callset, args.bias, args.offect)
		eve_record(Ans, Callset, args.bias, args.offect)
		svtype = ["DEL", "INS"]
		for i in svtype:
			caller_r, caller_tr = statistics_true_possitive(Callset, i)
			answer_r, answer_tr = statistics_true_possitive(Ans, i)
			sensitivity = 1.0*answer_tr/answer_r
			accuracy = 1.0*caller_tr/caller_r
			logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, caller_tr, caller_r, accuracy))
			logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, answer_tr, answer_r, sensitivity))
	caller_r, caller_tr = statistics_true_possitive(Callset)
	answer_r, answer_tr = statistics_true_possitive(Ans)
	sensitivity = 1.0*answer_tr/answer_r
	accuracy = 1.0*caller_tr/caller_r
	logging.info("%s: Accuracy rate is %d/%d(%.5f)."%(i, caller_tr, caller_r, accuracy))
	logging.info("%s: Sensitivity rate is %d/%d(%.5f)."%(i, answer_tr, answer_r, sensitivity))


def main(argv):
	args = parseArgs(argv)
	setupLogging(False)
	# print args
	starttime = time.time()
	main_ctrl(args)
	logging.info("Finished in %0.2f seconds."%(time.time() - starttime))

USAGE="""\
	Evaluate SV callset generated from NA12878
"""

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="NA12878_eval", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("CALLER", type=str, help="Choose a caller.")
	parser.add_argument("ANS", type=str, help="Golden answers.")
	parser.add_argument('CALL', type=str, help = "Callset produced by chosen caller.")
	parser.add_argument("CALL_c", type=str, help="Coverage of chosen callsets")
	parser.add_argument('-b', '--bias', help = "Bias of overlaping.[%(default)s]", default = 0.7, type = float)
	parser.add_argument('-o', '--offect', help = "Offect of INS breakpoint overlaping.[%(default)s]", default = 100, type = int)
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
